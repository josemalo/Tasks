"""
visualization_A4_summary.py
Input file has header.
Input file columns: seqID; Off-Target_Sequence_Mismatch_Only; Off-Target_Sequence_Bulge; Realigned_Off-Target_Sequence; Read_CountS.
"""

from __future__ import print_function
import collections
import svgwrite
import os
import argparse
import logging

logger = logging.getLogger('root')
logger.propagate = False

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', '-': '#B3B3B3'}

def parseSites(infile):
    offtarget_dict = collections.defaultdict(list)

    with open(infile, 'r') as f:
        header = f.readline().rstrip().split('\t')
        sample_names = header[4:]

        for line in f:
            read = line.rstrip('\n').split('\t')
            offtarget_dict[read[0]] = {'seqID': read[0], 'seq': read[1], 'seq_bulged': read[2], 'realigned_target_seq': read[3]}

            for i in range(len(sample_names)):
                sample_name = sample_names[i]
                offtarget_dict[read[0]][sample_name] = int(read[i+4])

    offtargets = [offtarget_dict[index] for index in offtarget_dict.keys()]
    offtargets = sorted(offtargets, key=lambda x: (x[sample_names[0]], x[sample_names[1]]), reverse=True)

    return offtargets, sample_names

def visualizeOfftargets(offtargets, sample_names, outfile, target_seq, title):
    # Initiate canvas
    dwg = svgwrite.Drawing(outfile + '.svg', profile='full', size=(u'100%', 100 + len(offtargets) * (box_size + 1)))

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20

    # Draw ticks
    tick_locations = [1, len(target_seq)]  # limits
    if target_seq.index('N') > len(target_seq)/2:  # PAM on the right end
        tick_locations += range(len(target_seq) + 1)[::10][1:]  # intermediate values
        tick_locations += range(len(target_seq) + 1)[len(target_seq) - 2: len(target_seq)]  # complementing PAM
        tick_locations.sort()
        tick_legend = [str(x) for x in tick_locations[:-3][::-1]] + ['P', 'A', 'M']
    else:
        tick_locations += [range(3, len(target_seq) + 1)[::10][1]]
        tick_locations += range(2, 5)
        tick_locations.sort()
        tick_legend = ['P', 'A', 'M'] + [str(x) for x in [str(x-3) for x in tick_locations[3:]]]

    for x,y in zip(tick_locations, tick_legend):
        dwg.add(dwg.text(y, insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    # Draw target sequence row
    for i, c in enumerate(target_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))

    for i in range(len(sample_names)):
        sample_name = sample_names[i]
        dwg.add(dwg.text(sample_name, insert=(x_offset + box_size * len(target_seq) + 16 + 70*i, y_offset + box_size - 3),
                         style="font-size:15px; font-family:Courier"))

    # Draw aligned sequence rows
    y_offset += 10  # leave some extra space after the reference row
    line_number = 0  # keep track of plotted sequences
    for j, seq in enumerate(offtargets):
        realigned_target_seq = offtargets[j]['realigned_target_seq']
        no_bulge_offtarget_sequence = offtargets[j]['seq']
        bulge_offtarget_sequence = offtargets[j]['seq_bulged']

        k = 0
        y = y_offset + j * box_size

        if no_bulge_offtarget_sequence != '':
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(no_bulge_offtarget_sequence, target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1

        if bulge_offtarget_sequence != '':
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(bulge_offtarget_sequence, realigned_target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(realigned_target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1

        # single off-target sequence
        if no_bulge_offtarget_sequence == '' or bulge_offtarget_sequence == '':
            for i in range(len(sample_names)):
                sample_name = sample_names[i]
                if seq[sample_name] == 0:
                    reads_sample = '-'
                else:
                    reads_sample = str(seq[sample_name])

                reads_text_sample = dwg.text(reads_sample,
                                             insert=(box_size * (len(target_seq) + 1) + 20 + i * 70, y_offset + box_size * (line_number + 2) - 2),
                                             fill='black', style="font-size:15px; font-family:Courier")
                dwg.add(reads_text_sample)
        else:
            for i in range(len(sample_names)):
                sample_name = sample_names[i]
                if seq[sample_name] == 0:
                    reads_sample = '-'
                else:
                    reads_sample = str(seq[sample_name])
                reads_text = dwg.text(reads_sample,
                                      insert=(box_size * (len(target_seq) + 1) + 20 + i * 70, y_offset + box_size * (line_number + 1) + 5),
                                      fill='black', style="font-size:15px; font-family:Courier")
                dwg.add(reads_text)
            reads_text02 = dwg.text(u"\u007D", insert=(
            box_size * (len(target_seq) + 1) + 7, y_offset + box_size * (line_number + 1) + 5),
                                    fill='black', style="font-size:23px; font-family:Courier")
            dwg.add(reads_text02)
    dwg.save()

def visualizeOffTargetSet(infile, outfile, target_seq, title, total_sites):
    # Get offtargets array from file
    offtargets, sample_names = parseSites(infile)

    counter = 0
    bk = 0
    sites = list()
    # Running 'visualizeOfftargets' for every site of 127 sites
    for site in offtargets:
        bk += 1
        if bk % int(total_sites) != 0 and bk <= len(offtargets):
            sites.append(site)
        if bk % int(total_sites) == 0:
            counter += 1
            visualizeOfftargets(sites, sample_names, '_'.join([outfile, str(counter)]), target_seq, title)
            sites = list()
        if bk == len(offtargets):
            counter += 1
            visualizeOfftargets(sites, sample_names, '_'.join([outfile, str(counter)]), target_seq, title)
            print('Last site has been reached')

def main():
    parser = argparse.ArgumentParser(description='Visualization plots with multiple field values')
    parser.add_argument('--infile', help='Input file with header: ots_mm, ots_bulge, ts, id, numeric_fields', required=True)
    parser.add_argument('--outfile', help='Output', required=True)
    parser.add_argument('--ts', help='Target sequence', required=True)
    parser.add_argument('--title', help='guide RNA', required=True)
    parser.add_argument('--total_sites', help='maximum number of sites per column', default=int(127))
    args = parser.parse_args()

    output_folder = os.path.dirname(args.outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    visualizeOffTargetSet(args.infile, args.outfile, args.ts, args.title, args.total_sites)


if __name__ == "__main__":
    main()
