"""
consolidation_circleseq.py
"""

from __future__ import print_function

___author___ = 'Jose Malagon-Lopez'

import argparse
import HTSeq
import pyfaidx
import regex
import os

### The following functions were taken from circleseq ###
""" Get sequences from some reference genome
"""
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return str(seq)

def reverseComplement(seq):
    compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', '.': '.', '-': '-', '_': '_'})
    out_list = [compl[bp] for bp in seq]
    return ''.join(out_list[::-1])

def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
    seq = seq.upper()
    """
    Given a sequence with ambiguous base characters, returns a regex that matches for
    the explicit (unambiguous) base characters
    """
    IUPAC_notation_regex = {'N': '[ATCGN]',
                            'Y': '[CTY]',
                            'R': '[AGR]',
                            'W': '[ATW]',
                            'S': '[CGS]',
                            'A': 'A',
                            'T': 'T',
                            'C': 'C',
                            'G': 'G'}

    pattern = ''

    for c in seq:
        pattern += IUPAC_notation_regex[c]

    if lookahead:
        pattern = '(?b:' + pattern + ')'

    pattern_standard = pattern + '{{s<={0}}}'.format(errors)
    pattern_gap = pattern + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
    return pattern_standard, pattern_gap


"""
Allow for '-' in our search, but do not allow insertions or deletions. 
"""
def extendedPattern(seq, errors):
    IUPAC_notation_regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
    realign_pattern = ''
    for c in seq:
        realign_pattern += IUPAC_notation_regex_extended[c]
    return '(?b:' + realign_pattern + ')' + '{{s<={0}}}'.format(errors)

"""
Recreate A!!! sequence in the window_sequence that matches the conditions given for the fuzzy regex. 
Currently only working for off-targets with at most one bulge !!! 
"""
def realignedSequences(targetsite_sequence, chosen_alignment, errors):
    match_sequence = chosen_alignment.group()
    substitutions, insertions, deletions = chosen_alignment.fuzzy_counts

    # get the .fuzzy_counts associated to the matching sequence after adjusting for indels, where 0 <= INS, DEL <= 1
    realigned_fuzzy = (substitutions, max(0, insertions - 1), max(0, deletions - 1))

    #  recreate the target sequence, with a '-' in the case of an DNA-bulge
    if insertions:
        targetsite_realignments = [targetsite_sequence[:i] + '-' + targetsite_sequence[i:] for i in range(1, len(targetsite_sequence))]
    else:
        targetsite_realignments = [targetsite_sequence]

    realigned_target_sequence, realigned_offtarget_sequence = None, ''  # in case the matching sequence is not founded

    for seq in targetsite_realignments:
        # recreate off-target sequence (with a '-' in the case of an RNA-bulge) and pattern matching the realigned target sequence
        if deletions:
            match_realignments = [match_sequence[:i + 1] + '-' + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
            match_pattern = [match_sequence[:i + 1] + seq[i + 1] + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
        else:
            match_realignments = match_pattern = [match_sequence]

        x = extendedPattern(seq, errors)
        for y_pattern, y_alignment in zip(match_pattern, match_realignments):
            m = regex.search(x, y_pattern, regex.BESTMATCH)
            if m and m.fuzzy_counts == realigned_fuzzy:
                realigned_target_sequence, realigned_offtarget_sequence = seq, y_alignment
    return realigned_target_sequence, realigned_offtarget_sequence


"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match(es).
"""
def alignSequences(targetsite_sequence, window_sequence, max_score=7):

    window_sequence = window_sequence.upper()
    query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_score)

    # Try both strands
    alignments_mm, alignments_bulge = list(), list()
    alignments_mm.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
    alignments_mm.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
    alignments_bulge.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
    alignments_bulge.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

    lowest_distance_score, lowest_mismatch = 100, max_score + 1
    chosen_alignment_b, chosen_alignment_m, chosen_alignment_strand_b, chosen_alignment_strand_m = None, None, '', ''

    # Use regex to find the best match allowing only for mismatches
    for aln_m in alignments_mm:
        strand_m, alignment_type_m, match_m = aln_m
        if match_m != None:
            mismatches, insertions, deletions = match_m.fuzzy_counts
            if mismatches < lowest_mismatch:
                chosen_alignment_m = match_m
                chosen_alignment_strand_m = strand_m
                lowest_mismatch = mismatches

    # Use regex to find the best match allowing for gaps, so that its edit distance is strictly lower than the
    # total number of mismatches of the sequence founded (if any) allowing only for mismatches.
    for aln_b in alignments_bulge:
        strand_b, alignment_type_b, match_b = aln_b
        if match_b != None:
            substitutions, insertions, deletions = match_b.fuzzy_counts
            if insertions or deletions:
                distance_score = substitutions + (insertions + deletions) * 3
                edistance = substitutions + insertions + deletions
                if distance_score < lowest_distance_score and edistance < lowest_mismatch:
                    chosen_alignment_b = match_b
                    chosen_alignment_strand_b = strand_b
                    lowest_distance_score = distance_score

    if chosen_alignment_m:
        offtarget_sequence_no_bulge = chosen_alignment_m.group()
        mismatches = chosen_alignment_m.fuzzy_counts[0]
        start_no_bulge = chosen_alignment_m.start()
        end_no_bulge = chosen_alignment_m.end()
    else:
        offtarget_sequence_no_bulge, mismatches, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '', '', '', '', ''

    bulged_offtarget_sequence, score, length, substitutions, insertions, deletions, bulged_start, bulged_end, realigned_target = \
        '', '', '', '', '', '', '', '', 'none'
    if chosen_alignment_b:
        realigned_target, bulged_offtarget_sequence = realignedSequences(targetsite_sequence, chosen_alignment_b, max_score)
        if bulged_offtarget_sequence:
            length = len(chosen_alignment_b.group())
            substitutions, insertions, deletions = chosen_alignment_b.fuzzy_counts
            score = substitutions + (insertions + deletions) * 3
            bulged_start = chosen_alignment_b.start()
            bulged_end = chosen_alignment_b.end()
        else:
            chosen_alignment_strand_b = ''

    return [offtarget_sequence_no_bulge, mismatches, len(offtarget_sequence_no_bulge), chosen_alignment_strand_m, start_no_bulge, end_no_bulge,
            realigned_target,
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end]


### New functions ###
def arrayOffTargets(identified_file, gap):
    #'Chromosome', 'Start', 'End', 'Name', 'ReadCount', 'Strand'  # 0:5
    #'MappingPositionStart', 'MappingPositionEnd', 'WindowName', 'WindowSequence',  # 6:9
    #'Site_SubstitutionsOnly.Sequence', 'Site_SubstitutionsOnly.NumSubstitutions',  # 10:11
    #'Site_SubstitutionsOnly.Strand', 'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',  # 12:14
    #'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',  # 15:17
    #'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',  # 18:20
    #'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',  #21:23
    #'FileName', 'Cell', 'Targetsite', 'FullName', 'TargetSequence', 'RealignedTargetSequence',  # 24:29
    #'Position.Pvalue', 'Narrow.Pvalue', 'Position.Control.Pvalue', 'Narrow.Control.Pvalue',  # 30:33
  
    g = open(identified_file, 'r')
    sites = g.readlines()[1:]
    g.close()

    ga = HTSeq.GenomicArray("auto", stranded=False, typecode='O')
    sites_dict = {}
    bk = list()
    sites_to_aggregated = {}

    #  collect sites and make genomic intervals for each of them using the Position coordinates
    for line in sites:
        site = line.strip().split('\t')
        Chromosome = site[0]
        PositionStart = int(site[6])
        PositionEnd = int(site[7])
        Name = site[3]

        #  intervals are left-closed and right-open
        ga[HTSeq.GenomicInterval(Chromosome, PositionStart, PositionEnd+1, ".")] = site
        sites_dict[Name] = site

    #  collect sites to be consolidated
    for name in sites_dict:
        Chromosome, Start, End, Name, ReadCount, Strand, PositionStart, PositionEnd, WindowName, WindowSequence, \
        Site_SubstitutionsOnly_Sequence, Site_SubstitutionsOnly_NumSubstitutions, \
        Site_SubstitutionsOnly_Strand, Site_SubstitutionsOnly_Start, Site_SubstitutionsOnly_End, \
        Site_GapsAllowed_Sequence, Site_GapsAllowed_Length, Site_GapsAllowed_Score, \
        Site_GapsAllowed_Substitutions, Site_GapsAllowed_Insertions, Site_GapsAllowed_Deletions, \
        Site_GapsAllowed_Strand, Site_GapsAllowed_Start, Site_GapsAllowed_End, \
        FileName, Cell, Targetsite, FullName, TargetSequence, RealignedTargetSequence, \
        Position_Pvalue, Narrow_Pvalue, Position_Control_Pvalue, Narrow_Control_Pvalue = sites_dict[name]

        if Name not in bk:
            bk += [Name]
            flag = True 
            chromosome = Chromosome
            start = int(PositionStart)
            end = int(PositionEnd)
            to_aggregated = [Name]  # to be used to collect the names of sites to be consolidated with the current site

            while flag:
                new_to_add = []  # to be used to check if during the loop new sites have been added to 'to_aggregated'
                for interval, value in ga[HTSeq.GenomicInterval(chromosome, start-gap, end+gap, ".")].steps():
                    if value is not None:
                        if value[3] not in to_aggregated:
                            new_to_add += [value[3]]
                            bk += [value[3]]
                            to_aggregated += [value[3]]
                            start = min(start, int(value[6]))
                            end = max(end, int(value[7]))

                if not new_to_add:
                    sites_to_aggregated[':'.join([chromosome, str(start)])] = [chromosome, start, end]
                    flag = False
    return ga, sites_to_aggregated


def consolidation(identified_file, max_score, reference, gap, search_radius):

    reference_genome = pyfaidx.Fasta(reference)

    basename = os.path.basename(identified_file)
    output_folder = os.path.dirname(identified_file)
    output_filename = "".join([output_folder, '/', basename[:-4], '_consolidated.txt'])

    output_file = open(output_filename, 'w')
    print('Chromosome', 'Start', 'End', 'Name', 'ReadCount', 'Strand',  # 0:5
          'MappingPositionStart', 'MappingPositionEnd', 'WindowName', 'WindowSequence',  # 6:9
          'Site_SubstitutionsOnly.Sequence', 'Site_SubstitutionsOnly.NumSubstitutions',  # 10:11
          'Site_SubstitutionsOnly.Strand', 'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',  # 12:14
          'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',  # 15:17
          'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',  # 18:20
          'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',  #21:23
          'FileName', 'Cell', 'Targetsite', 'FullName', 'TargetSequence', 'RealignedTargetSequence',  # 24:29
          'Position.Pvalue', 'Narrow.Pvalue', 'Position.Control.Pvalue', 'Narrow.Control.Pvalue',  # 30:33
          'Consolidated.Sites.Count',  # 34
          sep='\t', file=output_file)
    output_file.close()

    ga, sites_to_aggregated = arrayOffTargets(identified_file, gap)

    for tag in sites_to_aggregated:
        Chromosome, PositionStart, PositionEnd = sites_to_aggregated[tag]
        ReadCount, Position_Pvalue, Narrow_Pvalue, Position_Control_Pvalue, Narrow_Control_Pvalue = 0, 1, 1, 1, 1
        total_sites_consolidated = 0

        for interval, value in ga[HTSeq.GenomicInterval(Chromosome, int(PositionStart) - gap, int(PositionEnd) + gap, ".")].steps():
            if value is not None:
                total_sites_consolidated += 1

                PositionStart = min(int(PositionStart), int(value[6]))
                PositionEnd = max(int(PositionEnd), int(value[7]))
                ReadCount += int(value[4])
                FileName, Cell, Targetsite = value[24:27]
                TargetSequence = value[28]
                Position_Pvalue = min(Position_Pvalue, float(value[30]))
                Narrow_Pvalue = min(Narrow_Pvalue, float(value[31]))
                Position_Control_Pvalue = min(Position_Control_Pvalue, float(value[32]))
                Narrow_Control_Pvalue = min(Narrow_Control_Pvalue, float(value[33]))

        window_sequence = get_sequence(reference_genome, Chromosome, PositionStart - search_radius, PositionEnd + search_radius)

        offtarget_sequence_no_bulge, mismatches, len_offtarget_sequence_no_bulge, \
        chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
        realigned_target, \
        bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, \
        chosen_alignment_strand_b, bulged_start, bulged_end = \
            alignSequences(TargetSequence, window_sequence, max_score=max_score)

        # get genomic coordinates of sequences
        mm_start, mm_end, b_start, b_end = '', '', '', ''
        if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '+':
            mm_start = PositionStart - search_radius + int(start_no_bulge)
            mm_end = PositionStart - search_radius + int(end_no_bulge)
        if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '-':
            mm_start = PositionEnd + search_radius - int(end_no_bulge)
            mm_end = PositionEnd + search_radius - int(start_no_bulge)

        if bulged_offtarget_sequence and chosen_alignment_strand_b == '+':
            b_start = PositionStart - search_radius + int(bulged_start)
            b_end = PositionStart - search_radius + int(bulged_end)
        if bulged_offtarget_sequence and chosen_alignment_strand_b == '-':
            b_start = PositionEnd + search_radius - int(bulged_end)
            b_end = PositionEnd + search_radius - int(bulged_start)

        # define overall start and end position, only for annotation
        if offtarget_sequence_no_bulge:
            target_start_absolute = mm_start
            target_end_absolute = mm_end
            target_strand_absolute = chosen_alignment_strand_m
        elif not offtarget_sequence_no_bulge and bulged_offtarget_sequence:
            target_start_absolute = b_start
            target_end_absolute = b_end
            target_strand_absolute = chosen_alignment_strand_b
        else:
            target_start_absolute = iv.start
            target_end_absolute = iv.end
            target_strand_absolute = '*'

        Name = Chromosome + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
        FullName = str(Targetsite) + '_' + str(Cell) + '_' + str(Name) + '_' + str(ReadCount)
        WindowName = Chromosome + ':' + str(PositionStart) + '-' + str(PositionEnd)

        output_line = [Chromosome, target_start_absolute, target_end_absolute, Name, ReadCount, target_strand_absolute,
                       PositionStart, PositionEnd, WindowName,  window_sequence,
                       offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m, mm_start, mm_end,
                       bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
                       chosen_alignment_strand_b, b_start, b_end,
                       FileName, Cell, Targetsite, FullName, TargetSequence, realigned_target,
                       Position_Pvalue, Narrow_Pvalue, Position_Control_Pvalue, Narrow_Control_Pvalue, 
                       total_sites_consolidated]
        with open(output_filename, 'a') as output_file:
            print(*output_line, sep='\t', file=output_file)


def main():
    parser = argparse.ArgumentParser(description='Consolidate sites for very repetitive guides')
    parser.add_argument("--identified_file", help="", required=True)
    parser.add_argument("--max_score", help="", default=7)
    parser.add_argument("--reference", help="", required=True)
    parser.add_argument("--search_radius", help="", default=20)
    parser.add_argument("--gap", help="", default=10)
    args = parser.parse_args()

    print(args)

    consolidation(args.identified_file, int(args.max_score), args.reference, int(args.gap), int(args.search_radius))

if __name__ == "__main__":

    main()
