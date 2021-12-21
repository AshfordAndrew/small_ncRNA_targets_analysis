## Andrew Ashford, Podrabsky Lab, Thesis Project; mitosRNA/miRNA transcriptome target analysis
## This Python script will take in the miRanda and RNAhybrid predicted targets for the selected stress-response miRs, in
## addition to the mitochondrial sequences that are enriched in anoxia.

########## Import modules ##########
import csv
import matplotlib.pyplot as plt
import numpy as np

########## Specify file locations ##########
# miRanda and RNAhybrid miRNA predicted targets for A. limnaeus transcriptome:
miranda_mirna_targets_file = './sequence_files/miRanda_A_limnaeus_stress_response_miRNA_transcriptome_targets.txt'
rnahybrid_mirna_targets_file = './sequence_files/RNAhybrid_A_limnaeus_stress_response_miRNA_transcriptome_targets.txt'

# miRanda and RNAhybrid mitosRNA predicted targets for A. limnaeus transcriptome:
miranda_mitosrna_targets_file = './sequence_files/miRanda_A_limnaeus_mitosRNA_transcriptome_targets.txt'
rnahybrid_mitosrna_targets_file = './sequence_files/RNAhybrid_A_limnaeus_mitosRNA_transcriptome_targets.txt'

# Specify the location of the transcript sequence feature locations (5' UTR, CDS, 3' UTR etc.):
transcript_feature_locs_file = './A_limnaeus_mRNA_feature_locations.tsv'

# Specify the location of the "longest transcript per gene" file (for looking more on the gene level).
longest_transcript_per_gene_file = './sequence_files/A_limnaeus_longest_transcript_per_gene.tsv'

########## Read files into/create variables ##########
# Read in the sequence files using the csv module:
# miRanda miRNA transcriptome targets to dictionary (sorted by small RNA name, then transcript accession):
miranda_mirna = {}
with open(miranda_mirna_targets_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        if row and len(row) > 2:
            if '>' in row[0][0] and '>' not in row[0][1]:
                current_srna = row[0].replace('>', '')
                current_target_accession = row[1]
                if current_srna not in miranda_mirna.keys():
                    miranda_mirna[current_srna] = {}
                if current_target_accession not in miranda_mirna[current_srna].keys():
                    miranda_mirna[current_srna][current_target_accession] = row
#print(miranda_mirna)

# RNAhybrid miRNA transcriptome targets to dictionary (sorted by small RNA name, then transcript accession):
rnahybrid_mirna = {}
with open(rnahybrid_mirna_targets_file) as fd:
    rd = csv.reader(fd, delimiter=":", quotechar='"')
    for row in rd:
        current_srna = row[2]
        current_target_accession = row[0]
        if current_srna not in rnahybrid_mirna.keys():
            rnahybrid_mirna[current_srna] = {}
        if current_target_accession not in rnahybrid_mirna[current_srna].keys():
            rnahybrid_mirna[current_srna][current_target_accession] = row
#print(rnahybrid_mirna)

# miRanda mitosRNA transcriptome targets to dictionary (sorted by small RNA name, then transcript accession):
miranda_mitosrna = {}
with open(miranda_mitosrna_targets_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        if row and len(row) > 2:
            if '>' in row[0][0] and '>' not in row[0][1]:
                current_srna = row[0].replace('>', '')
                current_target_accession = row[1]
                if current_srna not in miranda_mitosrna.keys():
                    miranda_mitosrna[current_srna] = {}
                if current_target_accession not in miranda_mitosrna[current_srna].keys():
                    miranda_mitosrna[current_srna][current_target_accession] = row
#print(miranda_mitosrna)

# RNAhybrid mitosRNA transcriptome targets to dictionary (sorted by small RNA name, then transcript accession):
rnahybrid_mitosrna = {}
with open(rnahybrid_mitosrna_targets_file) as fd:
    rd = csv.reader(fd, delimiter=":", quotechar='"')
    for row in rd:
        current_srna = row[2]
        current_target_accession = row[0]
        if current_srna not in rnahybrid_mitosrna.keys():
            rnahybrid_mitosrna[current_srna] = {}
        if current_target_accession not in rnahybrid_mitosrna[current_srna].keys():
            rnahybrid_mitosrna[current_srna][current_target_accession] = row
#print(rnahybrid_mitosrna)

# Read feature locations file into a dictionary (with the transcript accession as the key) using csv module:
transcript_features = {}
with open(transcript_feature_locs_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    next(rd)
    for row in rd:
        if row[0] not in transcript_features.keys():
            transcript_features[row[0]] = row[1:]
#print(transcript_features)

# Read the transcripts from the "longest transcript per gene" file into a list using the csv module:
longest_transcript_per_gene_list = []
with open(longest_transcript_per_gene_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    next(rd)
    for row in rd:
        longest_transcript_per_gene_list.append(row[1])
#print(longest_transcript_per_gene_list)

########## Functions ##########
def predicted_target_transcript_locations(miranda_dict, rnahybrid_dict, transcript_features_dict, longest_transcript_per_gene_list, one_transcript_per_gene, per_mb, delta_g_cutoff):
    '''
    This function will take in a list (generated above) of the miRanda targets and the RNAhybrid targets. It will output
    statistics of the transcript locations of the targets given the cutoff.
    :param miranda_list:
    :param rnahybrid_list:
    :param delta_g_cutoff:
    :return:
    '''
    all_transcripts_dict = {}
    transcript_region_counts = {}
    # Iterate through the miRanda predicted targets and make sure they pass the delta_g_cutoff, additionally, make sure
    # that particular smallRNA has a match in the RNAhybrid target dictionary that ALSO passes the delta_g_cutoff.
    # Iterate through the initial keys in the dictionary, which represent the small RNAs:
    for smallrna in miranda_dict.keys():
        # Iterate through the next level of keys in the dictionary, which represent the transcript:
        for transcript in miranda_dict[smallrna].keys():
            # If "one transcript per gene" is True, make sure the current transcript is in the list:
            if one_transcript_per_gene and transcript not in longest_transcript_per_gene_list:
                continue
            # Check to make sure the delta G of the target are below the threshold set by "delta_g_cutoff" AND make sure
            # that the transcript is not actually a small RNA (XM_ vs XR_):
            if float(miranda_dict[smallrna][transcript][3]) <= delta_g_cutoff and 'XR_' not in transcript:
                # Check to make sure RNAhybrid also found a target for this small RNA/transcript combo:
                if smallrna in rnahybrid_dict.keys():
                    if transcript in rnahybrid_dict[smallrna].keys():
                        # Make sure that the small RNA/transcript combo predicted by RNAhybrid meets the maximum delta g
                        # value set by "delta_g_cutoff":
                        if float(rnahybrid_dict[smallrna][transcript][4]) <= delta_g_cutoff:
                            '''
                            print(smallrna)
                            print(transcript)
                            print(miranda_dict[smallrna][transcript])
                            print(rnahybrid_dict[smallrna][transcript])
                            '''
                            # Add transcript accession to the total transcript dictionary:
                            if smallrna not in all_transcripts_dict.keys():
                                all_transcripts_dict[smallrna] = []
                            if transcript not in all_transcripts_dict[smallrna]:
                                all_transcripts_dict[smallrna].append(transcript)

                            # Retrieve the sequence coordinates from the miRanda predicted target:
                            target_transcript_coordinates = miranda_dict[smallrna][transcript][6].split(' ')

                            # Convert list of strings from miRanda coordinates to list of ints:
                            target_transcript_coordinates = [int(i) for i in target_transcript_coordinates]
                            #print(target_transcript_coordinates)

                            # If the current small RNA is not in the output dictionary, create the output dictionary
                            # with the regions as keys:
                            if smallrna not in transcript_region_counts.keys():
                                transcript_region_counts[smallrna] = {'five_prime_utr': 0, 'coding_region': 0, 'three_prime_utr': 0}

                            # Get the corresponding features from the "transcript_features_dict" for the current small
                            # RNA and convert it from a list of strings to a list of ints:
                            current_transcript_feature_coords = [int(i) for i in transcript_features_dict[transcript]]
                            #print(current_transcript_feature_coords)

                            # Check the transcript target coordinates versus the feature coordinates and add the target
                            # region to the "transcript_region_counts" dictionary:
                            five_prime_target = 0
                            coding_sequence_target = 0
                            three_prime_target = 0
                            for i in range(0, len(target_transcript_coordinates)):
                                if target_transcript_coordinates[i] >= current_transcript_feature_coords[0] and target_transcript_coordinates[i] <= current_transcript_feature_coords[1]:
                                    five_prime_target += 1
                                elif target_transcript_coordinates[i] > current_transcript_feature_coords[2] and target_transcript_coordinates[i] <= current_transcript_feature_coords[3]:
                                    coding_sequence_target += 1
                                elif target_transcript_coordinates[i] > current_transcript_feature_coords[4] and target_transcript_coordinates[i] <= current_transcript_feature_coords[5]:
                                    three_prime_target += 1
                            if five_prime_target == 2:
                                transcript_region_counts[smallrna]['five_prime_utr'] += 1
                            elif coding_sequence_target == 2:
                                transcript_region_counts[smallrna]['coding_region'] += 1
                            elif three_prime_target == 2:
                                transcript_region_counts[smallrna]['three_prime_utr'] += 1
                            else:
                                if five_prime_target == 1:
                                    transcript_region_counts[smallrna]['five_prime_utr'] += 1
                                if coding_sequence_target == 1:
                                    transcript_region_counts[smallrna]['coding_region'] += 1
                                if three_prime_target == 1:
                                    transcript_region_counts[smallrna]['three_prime_utr'] += 1
    if per_mb:
        transcript_region_total_bases = {'five_prime_utr': 0, 'coding_region': 0, 'three_prime_utr': 0}
        if one_transcript_per_gene:
            for i in range(0, len(longest_transcript_per_gene_list)):
                # Iterate through the transcript features and calculate the total bases in the 5' UTR, the coding region, and
                # the 3' UTR:
                #print(transcript)
                #print(transcript_features_dict[transcript])
                transcript_region_total_bases['five_prime_utr'] += int(transcript_features_dict[longest_transcript_per_gene_list[i]][1]) - int(transcript_features_dict[longest_transcript_per_gene_list[i]][0])
                transcript_region_total_bases['coding_region'] += int(transcript_features_dict[longest_transcript_per_gene_list[i]][3]) - int(transcript_features_dict[longest_transcript_per_gene_list[i]][2])
                transcript_region_total_bases['three_prime_utr'] += int(transcript_features_dict[longest_transcript_per_gene_list[i]][5]) - int(transcript_features_dict[longest_transcript_per_gene_list[i]][4])
        else:
            for transcript in transcript_features_dict.keys():
                transcript_region_total_bases['five_prime_utr'] += int(transcript_features_dict[transcript][1]) - int(transcript_features_dict[transcript][0])
                transcript_region_total_bases['coding_region'] += int(transcript_features_dict[transcript][3]) - int(transcript_features_dict[transcript][2])
                transcript_region_total_bases['three_prime_utr'] += int(transcript_features_dict[transcript][5]) - int(transcript_features_dict[transcript][4])

        # Go through target counts and calculate the targets per kilobase (x/x*1000):
        for smallrna in transcript_region_counts.keys():
            transcript_region_counts[smallrna]['five_prime_utr'] = transcript_region_counts[smallrna]['five_prime_utr'] / (transcript_region_total_bases['five_prime_utr']/1000000)
            transcript_region_counts[smallrna]['coding_region'] = transcript_region_counts[smallrna]['coding_region'] / (transcript_region_total_bases['coding_region']/1000000)
            transcript_region_counts[smallrna]['three_prime_utr'] = transcript_region_counts[smallrna]['three_prime_utr'] / (transcript_region_total_bases['three_prime_utr']/1000000)

    return transcript_region_counts, all_transcripts_dict

########## Call functions ##########
transcript_region_counts_dict, all_transcripts = predicted_target_transcript_locations(miranda_mirna, rnahybrid_mirna, transcript_features, longest_transcript_per_gene_list, one_transcript_per_gene=False, per_mb=False, delta_g_cutoff=-25.0)
print(all_transcripts)
'''
for smallrna in transcript_region_counts_dict.keys():
    print(smallrna)
    print(transcript_region_counts_dict[smallrna])
for smallrna in all_transcripts.keys():
    print(smallrna)
    print(all_transcripts[smallrna])
'''

########## Make graphs ##########
# Make bar graphs of target regions:
x = []
y = []
z = []
k = []
# Populate y, z, and k from the dictionary results:
for smallrna in transcript_region_counts_dict.keys():
    if smallrna not in x:
        x.append(smallrna)
    y.append(transcript_region_counts_dict[smallrna]['five_prime_utr'])
    z.append(transcript_region_counts_dict[smallrna]['coding_region'])
    k.append(transcript_region_counts_dict[smallrna]['three_prime_utr'])

def subcategorybar(X, vals, width=0.8):
    n = len(vals)
    _X = np.arange(len(X))
    for i in range(n):
        plt.bar(_X - width / 2. + i / float(n) * width, vals[i], width=width / float(n), align="edge")
    plt.xticks(_X, X)

subcategorybar(x, [y, z, k])

plt.show()