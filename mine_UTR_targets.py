## Andrew Ashford, Podrabsky Lab, script to analyze the 3' and 5' UTR predicted RNAhybrid targets
## This script will take in a list of 1:1 longest transcript per gene, iterate through that list, and, given a delta G
## cutoff and a number of targets, get the resulting transcripts/figures with overlapping transcripts per small RNA.

##### Import modules #####
import csv
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

##### Specify input and output file locations #####
# Specify the location of the predicted targets file:
predicted_targets_utr_file = './sequence_files/RNAhybrid_A_limnaeus_miRNA_5UTR_targets.txt'
# Specify the location of the longest transcript per gene file:
longest_transcript_file = './sequence_files/A_limnaeus_longest_transcript_per_gene.tsv'
# Specify the location of the long and short-term proteomics universes:
lt_proteome_file = './proteomics_files/long-term_universe.xlsx'
st_proteome_file = './proteomics_files/short-term_universe.xlsx'
# Specify the location of the 'mRNA to proteins' conversion file from DAVID:
mrna_to_proteins_file = './proteomics_files/DAVID_A_limnaeus_mRNA_to_proteins.txt'
# Specify the location of the 'mRNA to gene symbol; conversion file:
mrna_to_gene_symbols_file = './sequence_files/A_limnaeus_RefSeq_mRNA_Accession_to_Gene_Symbol.txt'
# Specify the location of the 5' and 3' UTR sequence fasta files (for converting the transcript dict output to a "per kb
# version):
five_prime_utr_fasta_file = './sequence_files/A_limnaeus_transcriptome_five_prime_UTRs.fasta'
three_prime_utr_fasta_file = './sequence_files/A_limnaeus_transcriptome_three_prime_UTRs.fasta'
# Specify the location of the 'A. limnaeus to D. rerio accessions and GO IDs' file:
alim_to_drerio_accessions_file = './sequence_files/A_limnaeus_to_D_rerio_GO_IDs.tsv'

##### Create variables and read in files #####
# Read the longest transcripts into a list using the CSV module:
longest_transcript_list = []
with open(longest_transcript_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        longest_transcript_list.append(row[1])
#print(longest_transcript_list)

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

# Read the transcript to protein file into a dictionary with the transcript as the keys and the protein as the values:
mrna_to_proteins_dict = {}
with open(mrna_to_proteins_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        if row[0] not in mrna_to_proteins_dict.keys():
            mrna_to_proteins_dict[row[0]] = []
        mrna_to_proteins_dict[row[0]].append(row[1])
#print(mrna_to_proteins_dict)

# Just like the above dictionary only where the protein is the key and the transcript is the value:
proteins_to_mrna_dict = {}
with open(mrna_to_proteins_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        if row[1] not in proteins_to_mrna_dict.keys():
            proteins_to_mrna_dict[row[1]] = []
        proteins_to_mrna_dict[row[1]].append(row[0])
#print(proteins_to_mrna_dict)

# Read the transcript to gene symbol file where the key is the transcript ID and the value is the gene symbol:
mrna_to_gene_symbols_dict = {}
with open(mrna_to_gene_symbols_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        mrna_to_gene_symbols_dict[row[0] + '.1'] = row[1]
#print(mrna_to_gene_symbols_dict)

# Read the 'A. limnaeus to D. rerio accessions and GO IDs to a dataframe:
alim_to_drerio_accessions = pd.read_csv(alim_to_drerio_accessions_file, sep='\t')
#print(alim_to_drerio_accessions)

# Convert rows of pandas dataframe into searchable dictionary using only the columns containing protein names as the
# values and the columns containing the quantitative profile as the values:
lt_dict = dict(zip(lt_df['Accession Number'], lt_df['Quantitative Profile']))
st_dict = dict(zip(st_df['Accession Number'], st_df['Quantitative Profile']))

# Read in the 5' and 3' UTR fasta files:
five_prime_fasta_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(five_prime_utr_fasta_file, "fasta")}
three_prime_fasta_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(three_prime_utr_fasta_file, "fasta")}

##### Functions #####
def find_transcript_targets_within_threshold(longest_transcript_list, predicted_targets_dict, mrna_to_proteins_dict,
                                             delta_g_below, targets_per_transcript, dataset):
    '''
    This function will take in the list of the longest transcripts per gene, iterate through it, and then look through
    the predicted targets to determine whether that target is represented for each small RNA under the given threshold.
    '''
    output_dict = {}
    for small_rna in predicted_targets_dict.keys():
        print("Currently iterating through predicted targets for: " + small_rna)
        for transcript in predicted_targets_dict[small_rna].keys():
            if transcript in longest_transcript_list:
                # Variable to store the number of targets on a transcript that are above the given threshold:
                number_below_threshold = 0
                for i in range(0, len(predicted_targets_dict[small_rna][transcript])):
                    if float(predicted_targets_dict[small_rna][transcript][i][0]) <= delta_g_below:
                        number_below_threshold += 1
                # Check to see if the number of transcript targets above the targets was above a given threshold:
                if number_below_threshold >= targets_per_transcript:
                    # This occurs when the dataset is specified as "None", as in we just want to search the overall
                    # predicted targets outside of our proteome.
                    if not dataset:
                        if small_rna not in output_dict.keys():
                            output_dict[small_rna] = []
                        output_dict[small_rna].append(transcript)
                    # Since we want to look at either the short term or the long term dataset, we need to see if the
                    # corresponding transcripts' proteins are within them. We do this by converting the transcript to a
                    # list of its corresponding proteins, then searching the specified protein dataset's proteome. If
                    # there is a match in the proteome, we will add the transcript to the output.
                    else:
                        proteins_for_current_transcript = mrna_to_proteins_dict[transcript.split('.')[0]]
                        for protein in dataset:
                            if str(protein) != 'nan':
                                for i in range(0, len(proteins_for_current_transcript)):
                                    if proteins_for_current_transcript[i] in protein:
                                        if small_rna not in output_dict.keys():
                                            output_dict[small_rna] = []
                                        output_dict[small_rna].append(transcript)
    return output_dict

def convert_transcript_dict_to_protein_dict(transcript_dict, mrna_to_proteins_dict):
    protein_dict = {}
    for small_rna in transcript_targets_dict.keys():
        if small_rna not in protein_dict.keys():
            protein_dict[small_rna] = []
        for i in range(0, len(transcript_dict[small_rna])):
            protein_dict[small_rna].append(mrna_to_proteins_dict[transcript_dict[small_rna][i].replace('.1', '')][0])
    return protein_dict

def convert_transcript_dict_to_zebrafish_protein_dict(transcript_dict, mrna_to_gene_symbols_dict,
                                                          alim_to_drerio_accessions_df):
    output_dict = {}
    for small_rna in transcript_dict.keys():
        #print('Converting A. limnaeus transcripts to D. rerio proteins for: ' + str(small_rna))
        if small_rna not in output_dict.keys():
            output_dict[small_rna] = []
        current_transcript = transcript_dict[small_rna]
        for i in range(0, len(current_transcript)):
            print('Converting ' + str(i) + ' of ' + str(len(current_transcript)) + ' transcripts to D. rerio proteins'
                                                                                   ' for ' + str(small_rna) + '.')
            if current_transcript[i] in mrna_to_gene_symbols_dict.keys():
                current_gene_symbol = mrna_to_gene_symbols_dict[current_transcript[i]]
                for index, row in alim_to_drerio_accessions_df.iterrows():
                    if str(row['A_limnaeus_gene_name(s)']) != 'nan' and str(row['D_rerio_accession(s)']) != 'nan':
                        if current_gene_symbol in row['A_limnaeus_gene_name(s)'].split(', '):
                            output_dict[small_rna].append(row['D_rerio_accession(s)'].split(', ')[0])
                            break
    return output_dict

def plot_per_kb_distribution(transcript_targets_dict, utr_fasta_dict, proteins_to_mrna_dict, longest_transcript_list,
                             dataset):
    data = []
    labels = []
    # Get the total bases for the current dataset:
    total_bases_for_current_dataset = 0
    if dataset:
        for protein in dataset:
            if str(protein) != 'nan':
                if '|' not in protein and '#' not in protein:
                    current_protein = protein.split('.')[0]
                    for i in range(0, len(proteins_to_mrna_dict[current_protein])):
                        if proteins_to_mrna_dict[current_protein][i] + '.1' in longest_transcript_list and \
                                proteins_to_mrna_dict[current_protein][i] + '.1' in utr_fasta_dict.keys():
                            total_bases_for_current_dataset += len(utr_fasta_dict[(proteins_to_mrna_dict[current_protein][i]
                                                                              ) + '.1'])
                            break
    else:
        for i in range(1, len(longest_transcript_list)):
            if longest_transcript_list[i] in utr_fasta_dict.keys():
                total_bases_for_current_dataset += len(utr_fasta_dict[longest_transcript_list[i]])
    for small_rna in transcript_targets_dict.keys():
        print(small_rna)
        print(len(transcript_targets_dict[small_rna])/total_bases_for_current_dataset*1000)
        labels.append(small_rna)
        data.append(len(transcript_targets_dict[small_rna])/total_bases_for_current_dataset*1000)
    plt.xlabel('small RNA')
    plt.ylabel('number of targets per kbp')
    #plt.xticks(fontsize=10, rotation=45)
    plt.bar(labels, data)
    plt.show()

def plot_total_target_distribution(transcript_targets_dict):
    data = []
    labels = []
    for small_rna in transcript_targets_dict.keys():
        labels.append(small_rna)
        data.append(len(transcript_targets_dict[small_rna]))
    plt.xlabel('small RNA')
    plt.ylabel('number of targets')
    #plt.xticks(fontsize=10, rotation=45)
    plt.bar(labels, data)
    plt.show()

##### Call functions #####
# Specify dataset you are using (lt_dict, st_dict, or None), this is necessary for both generating the transcript
# targets dict, and for plotting:
dataset_to_use = None
# This variable should correspond with the current UTRs you're looking at (three_prime_fasta_dict is using 3' UTRs or
# five_prime_fasta_dict if using 5' UTRs):
utr_fasta = five_prime_fasta_dict
delta_g_cutoff = -24
targets_or_more_per_transcript = 1
transcript_targets_dict = find_transcript_targets_within_threshold(longest_transcript_list, predicted_targets_dict,
                                                                   mrna_to_proteins_dict, delta_g_below=delta_g_cutoff,
                                                                   targets_per_transcript=targets_or_more_per_transcript
                                                                   , dataset=dataset_to_use)
#protein_targets_dict = convert_transcript_dict_to_protein_dict(transcript_targets_dict, mrna_to_proteins_dict)
drer_protein_targets_dict = convert_transcript_dict_to_zebrafish_protein_dict(transcript_targets_dict,
                                                                              mrna_to_gene_symbols_dict,
                                                                              alim_to_drerio_accessions)

for small_rna in drer_protein_targets_dict.keys():
    print(small_rna)
    #print(drer_protein_targets_dict[small_rna])
    for i in range(0, len(drer_protein_targets_dict[small_rna])):
        print(drer_protein_targets_dict[small_rna][i])
    print(len(drer_protein_targets_dict[small_rna]))


'''
for small_rna in transcript_targets_dict.keys():
    print(small_rna)
    print(transcript_targets_dict[small_rna])
    print(len(transcript_targets_dict[small_rna]))
#print(transcript_targets_dict.values())
'''

##### Make plots #####
#plot_total_target_distribution(transcript_targets_dict)
#plot_per_kb_distribution(transcript_targets_dict, utr_fasta, proteins_to_mrna_dict,
#                         longest_transcript_list, dataset_to_use)

##### Write to files #####
