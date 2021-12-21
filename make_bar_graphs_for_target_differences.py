## Andrew Ashford, Podrabsky Lab, script to make bar graphs of the differences in targets counts between all, long-term,
## and short-term datasets to check for enrichment.

##### Import modules #####
import pandas as pd
import matplotlib.pyplot as plt

##### Specify file locations #####
# Specify the location of the targets per kb spreadsheet file:
targets_per_kb_spreadsheet_file = './all_miRNA_&_mitosRNA_target_counts_per_kB_at_-24.xlsx'

##### Create variables and read in files #####
# Read in the targets per kb excel file using the pandas module:
columns = ["small RNAs", "5' targets ALL", "5' targets LT", "5' targets ST", "5' long-term - all",
           "5' short-term - all", "3' targets ALL", "3' targets LT", "3' targets ST", "3' long-term - all",
           "3' short-term - all"]
targets_per_kb_spreadsheet_df = pd.read_excel(targets_per_kb_spreadsheet_file, names=columns)

##### Functions #####
def generate_five_prime_lt_minus_all_bar_graphs(df):
    df.plot(kind="bar", x="small RNAs", y=["5' long-term - all", "3' long-term - all"], stacked=True)
    plt.ylabel('targets per kb difference')
    plt.show()

def generate_five_prime_st_minus_all_bar_graphs(df):
    df.plot(kind="bar", x="small RNAs", y=["5' short-term - all", "3' short-term - all"], stacked=True)
    plt.ylabel('targets per kb difference')
    plt.show()

##### Create graphs #####
generate_five_prime_lt_minus_all_bar_graphs(targets_per_kb_spreadsheet_df)
generate_five_prime_st_minus_all_bar_graphs(targets_per_kb_spreadsheet_df)