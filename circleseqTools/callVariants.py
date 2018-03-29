from __future__ import print_function

import subprocess
import sys
import os
import argparse
import regex
import re
import HTSeq
import pyfaidx

"""
Taken from CIRCLEseq
"""
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

def extendedPattern(seq, errors):
    IUPAC_notation_regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
    realign_pattern = ''
    for c in seq:
        realign_pattern += IUPAC_notation_regex_extended[c]
    return '(?b:' + realign_pattern + ')' + '{{s<={0}}}'.format(errors)

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

def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return str(seq)


"""
Run samtools:mpileup and get all identified variants in the window sequences
"""
def snpCall(matched_file, reference, bam_file, out, search_radius):
    basename = os.path.basename(out)
    output_folder = os.path.dirname(out)

    # open matched file
    regions = list()
    with open(matched_file, 'rU') as f:
        f.readline()
        for line in f:
            site = line.strip().split('\t')
            #  chromosome, windowStart, windowEnd, strand, bam, region_basename (=Targetsite_Name)
            regions.append([site[0], int(site[6]) - search_radius, int(site[7]) + search_radius, '*', bam_file, '_'.join([site[26], site[3]])])

    print('Running samtools:mpileup for %s' % basename, file=sys.stderr)
    out_vcf = os.path.join(output_folder, basename + '_mpileup_output')
    if os.path.exists(out_vcf):
        subprocess.check_call('rm -r %s' % out_vcf, shell=True, env=os.environ.copy())
    os.makedirs(out_vcf)
    process_mpileup = open(os.path.join(out_vcf, 'logFile_mpileup'), 'w')

    for item in regions:
        chromosome, windowStart, windowEnd, strand, bam_file, region_basename = item
        region = '%s%s%s%s%s' % (chromosome, ":", int(windowStart), "-", int(windowEnd))
        output = os.path.join(out_vcf, region_basename + '.vcf')

        cl_vcf = 'samtools mpileup -v --region %s --fasta-ref %s %s > %s' % (region, reference, bam_file, output)
        subprocess.check_call(cl_vcf, shell=True, env=os.environ.copy(), stderr=process_mpileup, stdout=process_mpileup)
    process_mpileup.close()

    print('Collecting variants for %s' % basename, file=sys.stderr)
    out_bcf = os.path.join(output_folder, basename + '_output_bcftools')
    if os.path.exists(out_bcf):
        subprocess.check_call('rm -r %s' % out_bcf, shell=True, env=os.environ.copy())
    os.makedirs(out_bcf)
    process_bcftools = open(os.path.join(out_bcf, 'logFile_bcftools'), 'w')

    vcf_files = [f for f in os.listdir(out_vcf) if os.path.isfile(os.path.join(out_vcf, f))]
    for arch in vcf_files:
        if not arch.startswith('.') and arch.endswith('.vcf'):
            name = arch[:-4]
            output = os.path.join(out_bcf, name + '_BCFcall.vcf')

            cl_bcf = 'bcftools call -v -c %s > %s' % (os.path.join(out_vcf, arch), output)
            subprocess.check_call(cl_bcf, shell=True, env=os.environ.copy(), stderr=process_bcftools, stdout=process_bcftools)
    process_bcftools.close()

    print('Collecting significant variant calls for %s' % basename, file=sys.stderr)
    out_svc = os.path.join(output_folder, basename + '_output_svc')
    if os.path.exists(out_svc):
        subprocess.check_call('rm -r %s' % out_svc, shell=True, env=os.environ.copy())
    os.makedirs(out_svc)
    process_svc = open(os.path.join(out_svc, 'logFile_svc'), 'w')

    bcf_files = [f for f in os.listdir(out_bcf) if os.path.isfile(os.path.join(out_bcf, f))]
    for arch in bcf_files:
        if not arch.startswith('.') and arch.endswith('.vcf'):
            name = arch[:-12]
            output = os.path.join(out_svc, name + '_SIGNFcall.txt')

            cl_sed = "sed -n '/##/!p' %s | awk 'FNR>1' > %s" % (os.path.join(out_bcf, arch), output)
            subprocess.check_call(cl_sed, shell=True, env=os.environ.copy(), stderr=process_svc, stdout=process_svc)
    process_svc.close()

    print('Consolidating all the significant variant calls for %s' % basename, file=sys.stderr)
    header = ['targetsite', 'site_name', 'chromosome', 'one_based_position', 'reference', 'variant', 'quality', 'genotype', 'depth', 'PL']
    variants = list()

    svc_files = [f for f in os.listdir(out_svc) if os.path.isfile(os.path.join(out_svc, f))]
    for arch in svc_files:
        if not arch.startswith('.') and arch.endswith('.txt'):
            tag = arch[:-14]
            f = open(os.path.join(out_svc, arch), 'r')
            reads = f.readlines()
            f.close()

            for line in reads:
                item = line.split()
                if 'INDEL' in item[7]:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[3][3:]] + ['_'.join(item[9][4:].split(','))])
                else:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[0][3:]] + ['_'.join(item[9][4:].split(','))])

    out_file = open(out + '_mpileupCall.txt', 'w')
    print(*header, sep='\t', file=out_file)
    for item in variants:
        print(*item, sep='\t', file=out_file)
    out_file.close()

    print('Cleaning up directive for %s' % basename, file=sys.stderr)
    subprocess.check_call('rm -r %s' % out_vcf, shell=True, env=os.environ.copy())
    subprocess.check_call('rm -r %s' % out_bcf, shell=True, env=os.environ.copy())
    subprocess.check_call('rm -r %s' % out_svc, shell=True, env=os.environ.copy())

    print('Done running samtools:mpileup for %s' % basename, file=sys.stderr)
    return variants


"""
Obtain variant off-target sequences
"""
def realignVariantBulge(bulge_sequence, window_sequence_variant, bulge_strand):
    bseq = bulge_sequence.replace('-', '')
    if bulge_strand == '+':
        m_bulge = re.search(bseq, window_sequence_variant, re.I)
    else:
        m_bulge = re.search(bseq, reverseComplement(window_sequence_variant), re.I)
    variant_bseq = m_bulge.group()
    variant_bseq = variant_bseq[:bulge_sequence.find('-')] + '-' + variant_bseq[bulge_sequence.find('-'):]
    return variant_bseq


def SNPreader(snp_file):
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode='O')

    for snp in snp_file:
        basename, snpID, chromosome, one_based_position, reference, variant, quality, genotype, depth, PL = snp
        position = int(one_based_position) - 1
        key = '_'.join([basename, chromosome])
        ga[HTSeq.GenomicInterval(chromosome, position, position + 1, ".")] = [position, reference, variant, genotype, quality, key]
    return ga


def arrayOffTargets(matched_file, search_radius):
    offtargets_dict = {}
    gi_dict = {}

    with open(matched_file, 'r') as g:
        g.readline()
        for line in g:
            site = line.strip().split('\t')

            Chromosome = site[0]
            start = int(site[6]) - search_radius
            end = int(site[7]) + search_radius
            Name = site[3]

            offtargets_dict[Name] = site

            #  create a genomic interval for each window sequence
            gi_dict[Name] = HTSeq.GenomicInterval(Chromosome, start, end, ".")
    return offtargets_dict, gi_dict


def snpAdjustment(matched_file, snp_file, out, mismatch_threshold, search_radius):
    output_file = open(out + '_Variants.txt', 'w')
    print('Chromosome', 'Start', 'End', 'Name', 'ReadCount', 'Strand',
          'Variant_WindowSequence',
          'Variant_Site_SubstitutionsOnly.Sequence', 'Variant_Site_SubstitutionsOnly.NumSubstitutions',
          'Variant_Site_SubstitutionsOnly.Strand',
          'Variant_Site_GapsAllowed.Sequence', 'Variant_Site_GapsAllowed.Length', 
          'Variant_Site_GapsAllowed.Substitutions', 'Variant_Site_GapsAllowed.Insertions', 'Variant_Site_GapsAllowed.Deletions',
          'Variant_Site_GapsAllowed.Strand',
          'Cell', 'Targetsite', 'TargetSequence', 'Variant_RealignedTargetSequence',
          'Reference', 'Variant', 'Genotype', 'Quality',
          sep='\t', file=output_file)
    output_file.close()

    basename = os.path.basename(out)
    offtargets, gi_offtargets = arrayOffTargets(matched_file, search_radius)
    ga_snp = SNPreader(snp_file)

    for name in offtargets:
        variant_flag = False
        site = offtargets[name]
        gi = gi_offtargets[name]

        chromosome = site[0]
        window_sequence = site[9]
        window_sequence = window_sequence.upper()
        cell, targetsite = site[25:27]
        TargetSequence = site[28]
        output01 = site[0:6]
        output03 = [cell, targetsite, TargetSequence]
        ots_nb, ots_bu = site[10], site[15]

        #  obtain variant window sequence
        wkey = '_'.join([basename, chromosome])
        insert_start, insert_end, insert_var, snp_data = list(), list(), list(), {}

        for i, v in ga_snp[gi].steps():
            if v:
                position, reference, variant, genotype, quality, key = v
                if key == wkey:
                    variant = variant.split(',')[0]
                    for n, pos in enumerate(range(gi.start, gi.end)):
                        if pos == int(position):
                            insert_var.append(variant.lower())
                            insert_start.append(n)
                            end_pos = n + len(reference)
                            insert_end.append(end_pos)
                            snp_data[str(position)] = [position, reference, variant, genotype, quality]

        tri = 0
        window_sequence_variant = ''
        for i in range(len(insert_var)):
            variant = insert_var[i]
            pos = insert_start[i]
            window_sequence_variant += window_sequence[tri:pos] + variant.lower()
            tri = insert_end[i]
        window_sequence_variant += window_sequence[tri:]

        #  variant off-target sequences: only proceed if there is a variant in the window sequence
        window_sequence_var = window_sequence_variant.upper()
        if window_sequence_var != window_sequence:
            offtarget_sequence_no_bulge, mismatches, offtarget_sequence_length, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
            realigned_target, \
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end = \
                alignSequences(TargetSequence, window_sequence_var, max_score=mismatch_threshold)

            variant_ots_no_bulge, variant_ots_bulge = '', ''

            #  get variant sequence if the off-target sequences have changed by considering the variant window
            if ots_nb != offtarget_sequence_no_bulge:
                variant_flag = True
                if chosen_alignment_strand_m == '+':
                    m_no_bulge = re.search(offtarget_sequence_no_bulge, window_sequence_variant, re.I)
                else:
                    m_no_bulge = re.search(offtarget_sequence_no_bulge, reverseComplement(window_sequence_variant), re.I)
                variant_ots_no_bulge = m_no_bulge.group()

            if ots_bu != bulged_offtarget_sequence:
                variant_flag = True
                variant_ots_bulge = realignVariantBulge(bulged_offtarget_sequence, window_sequence_variant, chosen_alignment_strand_b)

            # collect and write variant data if we have variant off-target sequence(s)
            if variant_flag:
                total_genotype, total_reference, total_variant, total_quality = '', '', '', ''
                for pos in snp_data:
                    position, reference, variant, genotype, quality = snp_data[pos]
                    if total_genotype != '':
                        total_genotype += ''.join([':', genotype])
                        total_reference += ''.join([':', reference])
                        total_variant += ''.join([':', variant])
                        total_quality += ''.join([':', quality])
                    else:
                        total_genotype += ''.join([genotype])
                        total_reference += ''.join([reference])
                        total_variant += ''.join([variant])
                        total_quality += ''.join([quality])

                output02 = [variant_ots_no_bulge, mismatches, chosen_alignment_strand_m,
                            variant_ots_bulge, length, substitutions, insertions, deletions, chosen_alignment_strand_b]
                output04 = [total_reference, total_variant, total_genotype, total_quality]
                output_line = output01 + [window_sequence_variant] + output02 + output03 + [realigned_target] + output04

                with open(out + '_Variants.txt', 'a') as output_file:
                    print(*output_line, sep='\t', file=output_file)


"""
Main function
"""
def getVariants(matched_file, ref, bam_file, out, search_radius, mismatch_threshold):
    basename = os.path.basename(out)
    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    snp_file = snpCall(matched_file, ref, bam_file, out, search_radius)

    print('Obtaining Variant Off-Target Sequences for %s' % basename, file=sys.stderr)
    snpAdjustment(matched_file, snp_file, out, mismatch_threshold, search_radius)


def main():
    parser = argparse.ArgumentParser(description='Implement samtools:mpileup to identify genomic variants and adjust the off-target sequence when required.')
    parser.add_argument('--matched_file', help="full_path_to/matched file in 'identified' folder", required=True)
    parser.add_argument('--ref', help="Reference Genome Fasta", required=True)
    parser.add_argument('--bam', help="Sorted BAM file", required=True)
    parser.add_argument('--search_radius', help="Search radius around the position window", default=20, type=int)
    parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=7, type=int)
    parser.add_argument('--out', help="Output file basename, with full path", required=True)
    args = parser.parse_args()

    getVariants(args.matched_file, args.ref, args.bam, args.out, args.search_radius, args.mismatch_threshold)

if __name__ == "__main__":
    main()
