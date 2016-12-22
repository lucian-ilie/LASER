############################################################################
# Copyright (c) 2015-2016 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
import copy

from quast_libs import qconfig
from quast_libs.ca_utils.misc import is_same_reference, get_ref_by_chromosome

from quast_libs.log import get_logger
logger = get_logger(qconfig.LOGGER_DEFAULT_NAME)
from quast_libs.qutils import correct_name
from quast_libs.qutils import  icarus_report_str, _start, _end, short_str


class Misassembly:
    LOCAL = 0
    INVERSION = 1
    RELOCATION = 2
    TRANSLOCATION = 3
    INTERSPECTRANSLOCATION = 4  #for --meta, if translocation occurs between chromosomes of different references
    SCAFFOLD_GAP = 5
    FRAGMENTED = 6


class StructuralVariations(object):
    def __init__(self):
        self.inversions = []
        self.relocations = []
        self.translocations = []

    def get_count(self):
        return len(self.inversions) + len(self.relocations) + len(self.translocations)


#class Mapping(object):
#    def __init__(self, s1, e1, s2, e2, len1, len2, idy, ref, contig, ns_pos=None):
#        self[0], self[1], self[2], self[3], self[4], self[5], self[6], self[7], self[8] = s1, e1, s2, e2, len1, len2, idy, ref, contig
#        self.ns_pos = ns_pos

#    @classmethod
#    def from_line(cls, line):
        # line from coords file,e.g.
        # 4324128  4496883  |   112426   285180  |   172755   172756  |  99.9900  | gi|48994873|gb|U00096.2|	NODE_333_length_285180_cov_221082
#        line = line.split()
#        assert line[2] == line[5] == line[8] == line[10] == '|', line
#        contig = line[12]
#        ref = line[11]
#        s1, e1, s2, e2, len1, len2 = [int(line[i]) for i in [0, 1, 3, 4, 6, 7]]
#        idy = float(line[9])
#        return Mapping(s1, e1, s2, e2, len1, len2, idy, ref, contig)
#
#    def __str__(self):
#        return ' '.join(str(x) for x in [self[0], self[1], '|', self[2], self[3], '|', self[4], self[5], '|', self[6], '|', self[7], self[8]])
#
#    def short_str(self):
#        return ' '.join(str(x) for x in [self[0], self[1], '|', self[2], self[3], '|', self[4], self[5]])
#
#    def icarus_report_str(self, ambiguity='', is_best='True'):
#        return '\t'.join(str(x) for x in [self[0], self[1], self[2], self[3], self[7], self[8], self[6], ambiguity, is_best])
#
#    def clone(self):
#        return Mapping.from_line(str(self))

#    def start(self):
#        """Return start on contig (always <= end)"""
#        return min(self[2], self[3])

#    def end(self):
        """Return end on contig (always >= start)"""
#        return max(self[2], self[3])


class IndelsInfo(object):
    def __init__(self):
        self.mismatches = 0
        self.insertions = 0
        self.deletions = 0
        self.indels_list = []

    def __add__(self, other):
        self.mismatches += other.mismatches
        self.insertions += other.insertions
        self.deletions += other.deletions
        self.indels_list += other.indels_list
        return self


def distance_between_alignments(align1, align2, pos_strand1, pos_strand2, cyclic_ref_len=None):
    # returns distance (in reference) between two alignments
    distance_align1_align2 = align2[0] - align1[1] - 1
    distance_align2_align1 = align1[0] - align2[1] - 1
    if pos_strand1 and pos_strand2:            # alignment 1 should be earlier in reference
        distance = distance_align1_align2
    elif not pos_strand1 and not pos_strand2:  # alignment 2 should be earlier in reference
        distance = distance_align2_align1
    else:
        if align2[0] > align1[0]:
            distance = distance_align1_align2
        else:
            distance = distance_align2_align1

    cyclic_moment = False
    if cyclic_ref_len is not None:
        cyclic_distance = distance
        if align1[1] < align2[1] and (cyclic_ref_len + distance_align2_align1) < qconfig.extensive_misassembly_threshold:
            cyclic_distance = cyclic_ref_len + distance_align2_align1
        elif align1[1] >= align2[1] and (cyclic_ref_len + distance_align1_align2) < qconfig.extensive_misassembly_threshold:
            cyclic_distance = cyclic_ref_len + distance_align1_align2
        if abs(cyclic_distance) < abs(distance):
            distance = cyclic_distance
            cyclic_moment = True
    return distance, cyclic_moment


def __get_border_gaps(align1, align2, ref_lens):
    return [min(abs(align[1] - ref_lens[align[7]]), abs(align[0] - 1)) for align in [align1, align2]]


def is_fragmented_ref_fake_translocation(align1, align2, ref_lens):
    # Check whether translocation is caused by fragemented reference and thus should be marked Fake misassembly
    # Return inconsistency value if translocation is fake or None if translocation is real
    # !!! it is assumed that align1[7] != align2[7]
    assert align1[7] != align2[7], "Internal QUAST bug: is_fragmented_ref_fake_translocation() " \
                                     "should be called only if align1[7] != align2[7]"

    if qconfig.check_for_fragmented_ref:
        if qconfig.is_combined_ref and not is_same_reference(align1[7], align2[7]):
            return False
        if all([d <= qconfig.MAX_INDEL_LENGTH for d in __get_border_gaps(align1, align2, ref_lens)]):
            return True
    return False


def is_misassembly(align1, align2, contig_seq, ref_lens, is_cyclic=False, region_struct_variations=None):
    #Calculate inconsistency between distances on the reference and on the contig
#    distance_on_contig = align2.start() - align1.end() - 1
    distance_on_contig = _start(align2) - _end(align1) - 1
    cyclic_ref_lens = ref_lens if is_cyclic else None
    if cyclic_ref_lens is not None and align1[7] == align2[7]:
        distance_on_reference, cyclic_moment = distance_between_alignments(align1, align2, align1[2] < align1[3],
            align2[2] < align2[3], cyclic_ref_lens[align1[7]])
    else:
        distance_on_reference, cyclic_moment = distance_between_alignments(align1, align2, align1[2] < align1[3],
                                                                           align2[2] < align2[3])

    misassembly_internal_overlap = 0
    if distance_on_contig < 0:
        if distance_on_reference >= 0:
            misassembly_internal_overlap = (-distance_on_contig)
        elif (-distance_on_reference) < (-distance_on_contig):
            misassembly_internal_overlap = (distance_on_reference - distance_on_contig)

    strand1 = (align1[2] <= align1[3])
    strand2 = (align2[2] <= align2[3])
    inconsistency = distance_on_reference - distance_on_contig
    aux_data = {"inconsistency": inconsistency, "distance_on_contig": distance_on_contig,
                "misassembly_internal_overlap": misassembly_internal_overlap, "cyclic_moment": cyclic_moment,
                "is_sv": False, "is_scaffold_gap": False}

    if qconfig.scaffolds and contig_seq and check_is_scaffold_gap(inconsistency, contig_seq, align1, align2):
        aux_data["is_scaffold_gap"] = True
        return False, aux_data
    if region_struct_variations and check_sv(align1, align2, inconsistency, region_struct_variations):
        aux_data['is_sv'] = True
        return False, aux_data
    if align1[7] != align2[7] and is_fragmented_ref_fake_translocation(align1, align2, ref_lens):
        aux_data["inconsistency"] = sum(__get_border_gaps(align1, align2, ref_lens))
        return False, aux_data
    if align1[7] != align2[7] or abs(inconsistency) > qconfig.extensive_misassembly_threshold or strand1 != strand2:
        return True, aux_data
    return False, aux_data  # regular local misassembly


def check_sv(align1, align2, inconsistency, region_struct_variations):
    max_error = 100 # qconfig.smgap / 4  # min(2 * qconfig.smgap, max(qconfig.smgap, inconsistency * 0.05))
    max_gap = qconfig.extensive_misassembly_threshold / 4

    def __match_ci(pos, sv):  # check whether pos matches confidence interval of sv
        return sv[0] - max_error <= pos <= sv[1] + max_error

    if align2[0] < align1[0]:
        align1, align2 = align2, align1
    if align1[7] != align2[7]:  # translocation
        for sv in region_struct_variations.translocations:
            if sv[0][7] == align1[7] and sv[1][7] == align2[7] and \
                    __match_ci(align1[1], sv[0]) and __match_ci(align2[0], sv[1]):
                return True
            if sv[0][7] == align2[7] and sv[1][7] == align1[7] and \
                    __match_ci(align2[1], sv[0]) and __match_ci(align1[0], sv[1]):
                return True
    elif (align1[2] < align1[3]) != (align2[2] < align2[3]) and abs(inconsistency) < qconfig.extensive_misassembly_threshold:
        for sv in region_struct_variations.inversions:
            if align1[7] == sv[0][7] and \
                    (__match_ci(align1[0], sv[0]) and __match_ci(align2[0], sv[1])) or \
                    (__match_ci(align1[1], sv[0]) and __match_ci(align2[1], sv[1])):
                return True
    else:
        variations = region_struct_variations.relocations
        for index, sv in enumerate(variations):
            if sv[0][7] == align1[7] and __match_ci(align1[1], sv[0]):
                if __match_ci(align2[0], sv[1]):
                    return True
                # unite large deletion (relocations only)
                prev_end = sv[1][1]
                index_variation = index + 1
                while index_variation < len(variations) and \
                                        variations[index_variation][0][0] - prev_end <= max_gap and \
                                        variations[index_variation][0][7] == align1[7]:
                    sv = variations[index_variation]
                    if __match_ci(align2[0], sv[1]):
                        return True
                    prev_end = sv[1][1]
                    index_variation += 1
    return False


def find_all_sv(bed_fpath):
    if not bed_fpath:
        return None
    region_struct_variations = StructuralVariations()
    f = open(bed_fpath)
    for line in f:
        l = line.split('\t')
        if len(l) > 6 and not line.startswith('#'):
            try:
                align1 = Mapping(s1=int(l[1]), e1=int(l[2]), ref=correct_name(l[0]), s2=None, e2=None, len1=None, len2=None, idy=None, contig=None)
                align2 = Mapping(s1=int(l[4]), e1=int(l[5]),  ref=correct_name(l[3]), s2=None, e2=None, len1=None, len2=None, idy=None, contig=None)
                if align1[7] != align2[7]:
                    region_struct_variations.translocations.append((align1, align2))
                elif 'INV' in l[6]:
                    region_struct_variations.inversions.append((align1, align2))
                elif 'DEL' in l[6]:
                    region_struct_variations.relocations.append((align1, align2))
                else:
                    pass # not supported yet
            except ValueError:
                pass  # incorrect line format
    return region_struct_variations


def check_is_scaffold_gap(inconsistency, contig_seq, align1, align2):
    if abs(inconsistency) <= qconfig.scaffolds_gap_threshold and align1[7] == align2[7] and \
            is_gap_filled_ns(contig_seq, align1, align2) and (align1[2] < align1[3]) == (align2[2] < align2[3]):
        return True
    return False


def exclude_internal_overlaps(align1, align2, i=None, ca_output=None):
    # returns size of align1[5] decrease (or 0 if not changed). It is important for cur_aligned_len calculation

    def __shift_start(align, new_start, indent=''):
        if ca_output is not None:
            print >> ca_output.stdout_f, indent + '%s' % short_str(align),
#            print >> ca_output.stdout_f, indent + '%s' % align.short_str(),
        align = list(align)
        if align[2] < align[3]:
            align[0] += (new_start - align[2])
            align[2] = new_start
            align[5] = align[3] - align[2] + 1
        else:
            align[1] -= (new_start - align[3])
            align[3] = new_start
            align[5] = align[2] - align[3] + 1
        align[4] = align[1] - align[0] + 1
        align = tuple(align)
        if ca_output is not None:
            print >> ca_output.stdout_f, '--> %s' % short_str(align)
#            print >> ca_output.stdout_f, '--> %s' % align.short_str()

    def __shift_end(align, new_end, indent=''):
        if ca_output is not None:
            print >> ca_output.stdout_f, indent + '%s' % short_str(align),
#            print >> ca_output.stdout_f, indent + '%s' % align.short_str(),
        align = list(align)
        if align[2] < align[3]:
            align[1] -= (align[3] - new_end)
            align[3] = new_end
            align[5] = align[3] - align[2] + 1
        else:
            align[0] += (align[2] - new_end)
            align[2] = new_end
            align[5] = align[2] - align[3] + 1
        align[4] = align[1] - align[0] + 1
        align = tuple(align)
        if ca_output is not None:
            print >> ca_output.stdout_f, '--> %s' % short_str(align)
#            print >> ca_output.stdout_f, '--> %s' % align.short_str()

    if qconfig.ambiguity_usage == 'all':
        return 0

#    distance_on_contig = align2.start() - align1.end() - 1
    distance_on_contig = _start(align2) - _end(align1) - 1
    if distance_on_contig >= 0:  # no overlap
        return 0
    prev_len2 = align1[5]
    if ca_output is not None:
        print >> ca_output.stdout_f, '\t\t\tExcluding internal overlap of size %d between Alignment %d and %d: ' \
                                     % (-distance_on_contig, i+1, i+2),
    if qconfig.ambiguity_usage == 'one':  # left only one of two copies (remove overlap from shorter alignment)
        if align1[5] >= align2[5]:
#            __shift_start(align2, align1.end() + 1)
            __shift_start(align2, _end(align1) + 1)
        else:
#            __shift_end(align1, align2.start() - 1)
            __shift_end(align1, _start(align2) - 1)
    elif qconfig.ambiguity_usage == 'none':  # removing both copies
        if ca_output is not None:
            print >> ca_output.stdout_f
#        new_end = align2.start() - 1
#        __shift_start(align2, align1.end() + 1, '\t\t\t  ')
        new_end = _start(align2) - 1
        __shift_start(align2, _end(align1) + 1, '\t\t\t  ')
        __shift_end(align1, new_end, '\t\t\t  ')
    return prev_len2 - align1[5]


def count_ns_and_not_ns_between_aligns(contig_seq, align1, align2):
#    gap_in_contig = contig_seq[align1.end(): align2.start() - 1]
    gap_in_contig = contig_seq[_end(align1): _start(align2) - 1]
    ns_count = gap_in_contig.count('N')
    return ns_count, len(gap_in_contig) - ns_count


def is_gap_filled_ns(contig_seq, align1, align2):
#    gap_in_contig = contig_seq[align1.end(): align2.start() - 1]
    gap_in_contig = contig_seq[_end(align1): _start(align2) - 1]
    if len(gap_in_contig) < qconfig.Ns_break_threshold:
        return False
    return gap_in_contig.count('N')/len(gap_in_contig) > 0.95


def process_misassembled_contig(sorted_aligns, cyclic, aligned_lengths, region_misassemblies, ref_lens, ref_aligns,
                                ref_features, contig_seq, references_misassemblies, region_struct_variations,
                                misassemblies_matched_sv, ca_output, is_ambiguous=False):
    misassembly_internal_overlap = 0
    prev_align = sorted_aligns[0]
    cur_aligned_length = prev_align[5]
    is_misassembled = False
    contig_is_printed = False
    indels_info = IndelsInfo()
    contig_aligned_length = 0  # for internal debugging purposes

    for i in range(len(sorted_aligns) - 1):
        next_align = sorted_aligns[i + 1]

        cur_aligned_length -= exclude_internal_overlaps(prev_align, next_align, i, ca_output)
        is_extensive_misassembly, aux_data = is_misassembly(prev_align, next_align, contig_seq, ref_lens,
                                                            cyclic, region_struct_variations)
        inconsistency = aux_data["inconsistency"]
        distance_on_contig = aux_data["distance_on_contig"]
        misassembly_internal_overlap += aux_data["misassembly_internal_overlap"]
        cyclic_moment = aux_data["cyclic_moment"]
#        print >> ca_output.icarus_out_f, prev_align.icarus_report_str()
        print >> ca_output.icarus_out_f, icarus_report_str(prev_align)
        print >> ca_output.stdout_f, '\t\t\tReal Alignment %d: %s' % (i+1, str(prev_align))

        ref_aligns.setdefault(prev_align[7], []).append(prev_align)
        print >> ca_output.coords_filtered_f, str(prev_align)
        if aux_data["is_sv"]:
            print >> ca_output.stdout_f, '\t\t\t  Not a misassembly (structural variation of the genome) between these two alignments'
            print >> ca_output.icarus_out_f, 'fake: not a misassembly (structural variation of the genome)'
            misassemblies_matched_sv += 1
        elif aux_data["is_scaffold_gap"]:
            print >> ca_output.stdout_f, '\t\t\t  Incorrectly estimated size of scaffold gap between these two alignments:',
            print >> ca_output.stdout_f, 'gap length difference =', inconsistency
            region_misassemblies.append(Misassembly.SCAFFOLD_GAP)
            print >> ca_output.icarus_out_f, 'fake: scaffold gap size wrong estimation'
        elif is_extensive_misassembly:
            is_misassembled = True
            aligned_lengths.append(cur_aligned_length)
            contig_aligned_length += cur_aligned_length
            cur_aligned_length = 0
            if not contig_is_printed:
                print >> ca_output.misassembly_f, prev_align[8]
                contig_is_printed = True
            print >> ca_output.misassembly_f, 'Extensive misassembly (',
            print >> ca_output.stdout_f, '\t\t\t  Extensive misassembly (',
            if prev_align[7] != next_align[7]:  # it is not a Fake translocation, because is_extensive_misassembly is True
                if qconfig.is_combined_ref and \
                        not is_same_reference(prev_align[7], next_align[7]):  # if chromosomes from different references
                        region_misassemblies.append(Misassembly.INTERSPECTRANSLOCATION)
                        ref1, ref2 = get_ref_by_chromosome(prev_align[7]), get_ref_by_chromosome(next_align[7])
                        references_misassemblies[ref1][ref2] += 1
                        references_misassemblies[ref2][ref1] += 1
                        print >> ca_output.stdout_f, 'interspecies translocation',
                        print >> ca_output.misassembly_f, 'interspecies translocation',
                        print >> ca_output.icarus_out_f, 'interspecies translocation'
                else:
                    region_misassemblies.append(Misassembly.TRANSLOCATION)
                    print >> ca_output.stdout_f, 'translocation',
                    print >> ca_output.misassembly_f, 'translocation',
                    print >> ca_output.icarus_out_f, 'translocation'
            elif abs(inconsistency) > qconfig.extensive_misassembly_threshold:
                region_misassemblies.append(Misassembly.RELOCATION)
                print >> ca_output.stdout_f, 'relocation, inconsistency =', inconsistency,
                print >> ca_output.misassembly_f, 'relocation, inconsistency =', inconsistency,
                print >> ca_output.icarus_out_f, 'relocation, inconsistency =', inconsistency
            else: #if strand1 != strand2:
                region_misassemblies.append(Misassembly.INVERSION)
                print >> ca_output.stdout_f, 'inversion',
                print >> ca_output.misassembly_f, 'inversion',
                print >> ca_output.icarus_out_f, 'inversion'
            print >> ca_output.stdout_f, ') between these two alignments'
            print >> ca_output.misassembly_f, ') between %s %s and %s %s' % (prev_align[2], prev_align[3],
                                                                      next_align[2], next_align[3])
            ref_features.setdefault(prev_align[7], {})[prev_align[1]] = 'M'
            ref_features.setdefault(next_align[7], {})[next_align[1]] = 'M'
        else:
            reason_msg = "" + (" (linear representation of circular genome)" if cyclic_moment else "") + \
                         (" (fragmentation of reference genome)" if prev_align[7] != next_align[7] else "")
            if inconsistency == 0 and cyclic_moment:
                print >> ca_output.stdout_f, '\t\t\t  Not a misassembly' + reason_msg + ' between these two alignments'
                print >> ca_output.icarus_out_f, 'fake: not a misassembly' + reason_msg
            elif inconsistency == 0 and prev_align[7] != next_align[7]:  # is_fragmented_ref_fake_translocation is True, because is_extensive_misassembly is False
                print >> ca_output.stdout_f, '\t\t\t  Not a misassembly' + reason_msg + ' between these two alignments'
                region_misassemblies.append(Misassembly.FRAGMENTED)
                print >> ca_output.icarus_out_f, 'fake: not a misassembly' + reason_msg
            elif abs(inconsistency) <= qconfig.MAX_INDEL_LENGTH and \
                    count_ns_and_not_ns_between_aligns(contig_seq, prev_align, next_align)[1] <= qconfig.MAX_INDEL_LENGTH:
                ns_number, not_ns_number = count_ns_and_not_ns_between_aligns(contig_seq, prev_align, next_align)

                if inconsistency == 0:
                    print >> ca_output.stdout_f, ('\t\t\t  Short stretch of %d mismatches and %d Ns between these two alignments' % (not_ns_number, ns_number)) + reason_msg
                    indels_info.mismatches += not_ns_number
                    print >> ca_output.icarus_out_f, 'indel: stretch of mismatches and Ns' + reason_msg
                else:
                    indel_length = abs(inconsistency)
                    indel_class = 'Indel (<= 5bp)' if indel_length <= qconfig.SHORT_INDEL_THRESHOLD else 'Indel (> 5bp)'
                    indel_type = 'insertion' if inconsistency < 0 else 'deletion'
                    mismatches = max(0, not_ns_number - indel_length)
                    print >> ca_output.stdout_f, ('\t\t\t  %s between these two alignments: %s of length %d; %d mismatches'
                                                 % (indel_class, indel_type, indel_length, mismatches)) + reason_msg
                    indels_info.indels_list.append(indel_length)
                    if indel_type == 'insertion':
                        indels_info.insertions += indel_length
                    else:
                        indels_info.deletions += indel_length
                    indels_info.mismatches += mismatches
                    print >> ca_output.icarus_out_f, 'indel: ' + indel_class.lower() + reason_msg
            else:
                if qconfig.strict_NA:
                    aligned_lengths.append(cur_aligned_length)
                    contig_aligned_length += cur_aligned_length
                    cur_aligned_length = 0

                if inconsistency < 0:
                    #There is an overlap between the two alignments, a local misassembly
                    print >> ca_output.stdout_f, '\t\t\t  Overlap between these two alignments (local misassembly).',
                else:
                    #There is a small gap between the two alignments, a local misassembly
                    print >> ca_output.stdout_f, '\t\t\t  Gap between these two alignments (local misassembly).',
                    #print >> plantafile_out, 'Distance on contig =', distance_on_contig, ', distance on reference =', distance_on_reference
                print >> ca_output.stdout_f, 'Inconsistency = ' + str(inconsistency) + reason_msg
                print >> ca_output.icarus_out_f, 'local misassembly'
                region_misassemblies.append(Misassembly.LOCAL)

        prev_align = next_align
        cur_aligned_length += prev_align[5] - (-distance_on_contig if distance_on_contig < 0 else 0)

    #Record the very last alignment
    i = len(sorted_aligns) - 1
    print >> ca_output.stdout_f, '\t\t\tReal Alignment %d: %s' % (i + 1, str(next_align))
#    print >> ca_output.icarus_out_f, next_align.icarus_report_str()
    print >> ca_output.icarus_out_f, icarus_report_str(next_align)
    ref_aligns.setdefault(next_align[7], []).append(next_align)
    print >> ca_output.coords_filtered_f, str(next_align)
    aligned_lengths.append(cur_aligned_length)
    contig_aligned_length += cur_aligned_length

    assert contig_aligned_length <= len(contig_seq), "Internal QUAST bug: contig aligned length is greater than " \
                                                     "contig length (contig: %s, len: %d, aligned: %d)!" % \
                                                     (sorted_aligns[0][8], contig_aligned_length, len(contig_seq))

    return is_misassembled, misassembly_internal_overlap, references_misassemblies, indels_info, misassemblies_matched_sv
