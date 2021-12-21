## Andrew Ashford, Podrabsky Lab, Pull sequences for 5' and 3' UTRs
## This script will take in the feature location tsv file and the A. limnaeus transcriptome file, and pull the
## sequence at the coordinates for the 3' UTR and the 5' UTR per transcript.

##### Import modules #####
import csv
from Bio import SeqIO

##### Specify file locations #####
# Specify the location of the A. limnaeus transcriptome file:
transcriptome_file = './sequence_files/A_limnaeus_transcriptome.fna'
# Specify the location of the transcript feature coordinates file:
features_file = './A_limnaeus_mRNA_feature_locations.tsv'

##### Create variables and read in files #####
# Read the transcriptome file into a dictionary:
transcriptome_dict = SeqIO.to_dict(SeqIO.parse(transcriptome_file, "fasta"))
'''
print(transcriptome_dict['XM_013999454.1'])
print(transcriptome_dict['XM_013999454.1'].id)
print(transcriptome_dict['XM_013999454.1'].description)
print(transcriptome_dict['XM_013999454.1'].seq)
'''

# Read the features file into a dictionary:
transcriptome_features = {}
with open(features_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    next(rd)
    for row in rd:
        transcriptome_features[row[0]] = row[1:]
#print(transcriptome_features)

##### Functions #####
def get_5_prime_utrs(transcriptome_dict, transcriptome_features):
    output_fasta_string = ''
    # Iterate through the transcripts in features dict:
    for transcript in transcriptome_features.keys():
        if transcriptome_features[transcript][1] != '0' and 'XR_' not in transcript:
            current_seq = transcriptome_dict[transcript].seq[int(transcriptome_features[transcript][0]):
                                                             int(transcriptome_features[transcript][1])]
            if output_fasta_string == '':
                output_fasta_string = '>' + str(transcriptome_dict[transcript].description) + '\n' + str(current_seq)
            else:
                output_fasta_string += '\n' + '>' + str(transcriptome_dict[transcript].description) + '\n' + \
                                       str(current_seq)
    return output_fasta_string

def get_3_prime_utrs(transcriptome_dict, transcriptome_features):
    output_fasta_string = ''
    # Iterate through the transcripts in features dict:
    for transcript in transcriptome_features.keys():
        if transcriptome_features[transcript][5] != '0' and 'XR_' not in transcript:
            current_seq = transcriptome_dict[transcript].seq[int(transcriptome_features[transcript][4]):
                                                             int(transcriptome_features[transcript][5])]
            if output_fasta_string == '':
                output_fasta_string = '>' + str(transcriptome_dict[transcript].description) + '\n' + str(current_seq)
            else:
                output_fasta_string += '\n' + '>' + str(transcriptome_dict[transcript].description) + '\n' + \
                                       str(current_seq)
    return output_fasta_string

##### Call functions #####
# Specify output locations for the 5' and 3' transcriptome fasta files:
five_prime_fasta_output = './sequence_files/A_limnaeus_transcriptome_five_prime_UTRs.fasta'
three_prime_fasta_output = './sequence_files/A_limnaeus_transcriptome_three_prime_UTRs.fasta'

# Call functions and save their returns in their own individual variables:
five_prime_utrs = get_5_prime_utrs(transcriptome_dict, transcriptome_features)
print(five_prime_utrs)
three_prime_utrs = get_3_prime_utrs(transcriptome_dict, transcriptome_features)
print(three_prime_utrs)

# Write the fasta files to the above specified locations:
with open(five_prime_fasta_output, 'w') as fh:
    fh.write(five_prime_utrs)
with open(three_prime_fasta_output, 'w') as fh:
    fh.write(three_prime_utrs)