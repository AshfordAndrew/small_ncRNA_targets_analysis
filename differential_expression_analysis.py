## Andrew Ashford, Podrabsky Lab, script to perform differential expression analysis

##### Import modules #####
import csv
import pandas as pd
import matplotlib.pyplot as plt

##### Specify input and output file locations #####
# Specify the location of the predicted targets file:
predicted_targets_utr_file = './sequence_files/RNAhybrid_A_limnaeus_miRNA_5UTR_targets.txt'

# Specify the location of the long and short-term proteomics universes:
lt_proteome_file = './proteomics_files/long-term_universe.xlsx'
st_proteome_file = './proteomics_files/short-term_universe.xlsx'

# Specify the location of the 'mRNA to proteins' conversion file from DAVID:
mrna_to_proteins_file = './proteomics_files/DAVID_A_limnaeus_mRNA_to_proteins.txt'

##### Create variables and read in files #####
# Read the predicted targets into a dictionary of genes where the outermost key is the small RNA name, the inner key is
# the transcript name, at this value there will be a list with the target parameters:
predicted_targets_dict = {}
with open(predicted_targets_utr_file) as fd:
    rd = csv.reader(fd, delimiter=":", quotechar='"')
    for row in rd:
        if row[2] not in predicted_targets_dict.keys():
            predicted_targets_dict[row[2]] = {}
        if row[0] not in predicted_targets_dict[row[2]].keys():
            predicted_targets_dict[row[2]][row[0]] = []
        predicted_targets_dict[row[2]][row[0]].append(row[4:])

# Read in long and short term Excel datasets using the Pandas module:
lt_df = pd.read_excel(lt_proteome_file)
st_df = pd.read_excel(st_proteome_file)

# Alter the "Accession Number" column of the short-term Pandas dataframe to only include the "XP" accession number:
st_df['Accession Number'] = st_df['Accession Number'].str.split('|').str.get(3)
st_df['Quantitative Profile'] = st_df['Quantitative Profile'].replace('[]', 'no change')

# Read the transcript to protein file into a dictionary with the transcript as the keys and the protein as the values:
mrna_to_proteins_dict = {}
with open(mrna_to_proteins_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        if row[0] not in mrna_to_proteins_dict.keys():
            mrna_to_proteins_dict[row[0]] = []
        mrna_to_proteins_dict[row[0]].append(row[1])

# Convert rows of pandas dataframe into searchable dictionary using only the columns containing protein names as the
# values and the columns containing the quantitative profile as the values:
temp_lt_dict = dict(zip(lt_df['Accession Number'], lt_df['Quantitative Profile']))
st_dict = dict(zip(st_df['Accession Number'], st_df['Quantitative Profile']))
lt_dict = {}

# Fix the Quantitative Profiles of the temp_lt_dict and put the fixed dict into the lt_dict variable:
for accession in temp_lt_dict.keys():
    if temp_lt_dict[accession] != '[]':
        current_quant_profile = temp_lt_dict[accession].split(', ')[0] + ', ' + temp_lt_dict[accession].split(', ')[2] + ', ' + temp_lt_dict[accession].split(', ')[3] + ', ' + temp_lt_dict[accession].split(', ')[1]
        lt_dict[accession] = current_quant_profile
    else:
        lt_dict[accession] = 'no change'

##### Functions #####
def create_differential_expression_dict(predicted_targets_dict, mrna_to_proteins_dict, which_small_rna,
                                        which_proteomics, delta_g_cutoff, targets_cutoff):
    quantitative_profile_dict = {}
    for transcript in predicted_targets_dict[which_small_rna].keys():
        current_transcript_targets = predicted_targets_dict[which_small_rna][transcript]
        number_targets_in_threshold = 0
        for i in range(0, len(current_transcript_targets)):
            if float(current_transcript_targets[i][0]) <= delta_g_cutoff:
                number_targets_in_threshold += 1
        if number_targets_in_threshold >= targets_cutoff:
            current_transcript_proteins = mrna_to_proteins_dict[transcript.split('.')[0]]
            for i in range(0, len(current_transcript_proteins)):
                current_protein = current_transcript_proteins[i] + '.1'
                if current_protein in which_proteomics.keys():
                    if which_proteomics[current_protein] not in quantitative_profile_dict.keys():
                        quantitative_profile_dict[which_proteomics[current_protein]] = []
                    quantitative_profile_dict[which_proteomics[current_protein]].append(current_protein)
                    break
    return quantitative_profile_dict

def plot_quantitative_profile(quantitative_profile):
    data = []
    labels = []
    for profile in quantitative_profile.keys():
        labels.append(profile)
        data.append(len(quantitative_profile[profile]))
    plt.xlabel('quantitative profile')
    plt.ylabel('predicted targets in dataset')
    plt.xticks(fontsize=10, rotation=45)
    plt.bar(labels, data)
    plt.show()

##### Call functions #####
which_small_rna = "gene-tools-scramble"
delta_g_cutoff = -24.0
which_proteomics = lt_dict
targets_cutoff = 1
quantitative_profile = create_differential_expression_dict(predicted_targets_dict, mrna_to_proteins_dict,
                                                           which_small_rna, which_proteomics, delta_g_cutoff,
                                                           targets_cutoff)
plot_quantitative_profile(quantitative_profile)
