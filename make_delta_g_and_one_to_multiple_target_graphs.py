## Andrew Ashford, Podrabsky Lab, Master's Thesis
## This work will take in the RNAhybrid target results and will generate a graph where the x-axis has decreasing delta-G
## values, and the y-axis has number of targets. Another graph will be generated that shows the effect of increasing the
## number of targets at each delta-G value (for instance, the effects of requiring 2 targets with a dG of ~-19, vs 3
## targets at the same dG value etc).

##### Import modules #####
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

##### Specify file locations #####
# Specify the location of the predicted targets file:
predicted_targets_utr_file = './sequence_files/RNAhybrid_A_limnaeus_mitosRNA_5UTR_targets.txt'
# Specify the location of the longest transcript per gene file:
longest_transcript_file = './sequence_files/A_limnaeus_longest_transcript_per_gene.tsv'

##### Create variables and read in files #####
# Read the longest transcripts into a list using the CSV module:
longest_transcript_set = []
with open(longest_transcript_file) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        longest_transcript_set.append(row[1])
longest_transcript_set = set(longest_transcript_set)
#print(longest_transcript_set)

# Read the RNAhybrid predicted targets into a dictionary (small RNA -> mRNA accession -> list of delta-G values):
predicted_targets_dict = {}
with open(predicted_targets_utr_file) as fd:
    rd = csv.reader(fd, delimiter=":", quotechar='"')
    for row in rd:
        current_small_rna = row[2]
        current_mrna = row[0]
        current_delta_g = float(row[4])
        if current_small_rna not in predicted_targets_dict.keys():
            predicted_targets_dict[current_small_rna] = {}
        if current_mrna not in predicted_targets_dict[current_small_rna].keys():
            predicted_targets_dict[current_small_rna][current_mrna] = []
        predicted_targets_dict[current_small_rna][current_mrna].append(current_delta_g)
#print(predicted_targets_dict)

##### Functions #####
def generate_delta_g_graph_data_total(longest_transcript_set, predicted_targets_dict, delta_g_low, delta_g_high):
    output_graph_dict = {}
    # Specify the range of delta-G values you want on the x-axis of the graph:
    for i in range(delta_g_low, delta_g_high):
        print('Gathering results for the current dG threshold: ' + str(i))
        for small_rna in predicted_targets_dict.keys():
            # Add the current small RNA and delta-G threshold to the output dictionary for the graph:
            if small_rna not in output_graph_dict.keys():
                output_graph_dict[small_rna] = {}
            if i not in output_graph_dict[small_rna].keys():
                output_graph_dict[small_rna][i] = 0
            for mrna in predicted_targets_dict[small_rna].keys():
                if mrna in longest_transcript_set:
                    ###########################
                    total_within_threshold = 0
                    ###########################
                    for l in range(0, len(predicted_targets_dict[small_rna][mrna])):
                        current_target_dg = predicted_targets_dict[small_rna][mrna][l]
                        ##### denotes 2+ targets part of script #####
                        if float(i) >= current_target_dg:
                            total_within_threshold += 1
                            if total_within_threshold >= 2:
                                output_graph_dict[small_rna][i] += 1
                        #############################################
                        '''
                        if float(i) >= current_target_dg:
                            output_graph_dict[small_rna][i] += 1
                        '''
    return output_graph_dict

def generate_delta_g_plot(dg_dict):
    for small_rna in dg_dict.keys():
        current_x_axis = []
        current_y_axis = []
        for delta_g_value in dg_dict[small_rna].keys():
            current_x_axis.append(int(delta_g_value))
            current_y_axis.append(int(dg_dict[small_rna][delta_g_value]))
        plt.plot(current_x_axis, current_y_axis, label=small_rna)
    plt.xlabel('target count cutoff')
    plt.ylabel('number of targets')
    ax = plt.gca()
    #ax.invert_xaxis()
    ax.margins(x=0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    plt.legend(fontsize=10, loc='upper right')
    plt.show()

def generate_target_cutoff_graph_dictionary(longest_transcript_set, predicted_targets_dict, delta_g_cutoff, target_low,
                                            target_high):
    output_graph_dict = {}
    for i in range(target_low, target_high):
        print('Gathering results for the current target count threshold at -24 dG: ' + str(i))
        for small_rna in predicted_targets_dict.keys():
            # Add the current small RNA and delta-G threshold to the output dictionary for the graph:
            if small_rna not in output_graph_dict.keys():
                output_graph_dict[small_rna] = {}
            if i not in output_graph_dict[small_rna].keys():
                output_graph_dict[small_rna][i] = 0
            for mrna in predicted_targets_dict[small_rna].keys():
                if mrna in longest_transcript_set:
                    ###########################
                    total_within_threshold = 0
                    ###########################
                    for l in range(0, len(predicted_targets_dict[small_rna][mrna])):
                        current_target_dg = predicted_targets_dict[small_rna][mrna][l]
                        if delta_g_cutoff >= current_target_dg:
                            total_within_threshold += 1
                    if total_within_threshold >= (i):
                        output_graph_dict[small_rna][i] += total_within_threshold

    return output_graph_dict

##### Call functions #####
'''
dg_high = -15
dg_low = -32
dg_graph_dict = generate_delta_g_graph_data_total(longest_transcript_set, predicted_targets_dict, dg_low, dg_high)
print(dg_graph_dict)
'''

targets_low = 1
targets_high = 20
delta_g_cutoff = -23.0
targets_graph_dict = generate_target_cutoff_graph_dictionary(longest_transcript_set, predicted_targets_dict,
                                                             delta_g_cutoff, targets_low, targets_high)
print(targets_graph_dict)


##### Generate graphs #####
#generate_delta_g_plot(dg_graph_dict)
generate_delta_g_plot(targets_graph_dict)