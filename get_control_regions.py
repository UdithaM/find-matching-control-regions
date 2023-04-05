import pandas as pd
import random
import sys
import os
import time
from pandas.io.common import EmptyDataError



# parameters
annotation_filename = sys.argv[1] 
corsiv_filename = sys.argv[2]
output_filename = sys.argv[3]
# get chromosome in the type 'chr4'
chromosome = sys.argv[4]
genic = sys.argv[5]

# # parameters
# annotation_filename = "../../Annotations/chr21_annotated.csv"
# corsiv_filename = 'NEW_DATASET/corsiv_regions_that_has_missing_genic_control.csv'
# output_filename = 'NEW_DATASET/chrom_21_controls_new.csv'
# # get chromosome in the type 'chr4'
# chromosome = 'chr21'
# genic = 'genic'

def find_match(cpg, size, bin_dict):
    """
    Finds a control region in the input list with the input size.
    
    Input:
    cpg - the number of CpGs in the CoRSIV
    size - an integer representing the desired size for the control region
    bin_dict - a dictionary with keys being the start coordinate of the bin and values being the number of CpGs in the bin

    Returns:
    start coord and end coord of the control region
    """
    done = False
    count = 0
    # if after 200 times we still can't find a control, we will loosen up the requirement and try again, can change 200 to whatever number
    while not done and count <= 200:
        # randomly pick a bin
        coord, bin_cpg = random.choice(list(bin_dict.items()))
        # start at a bin with cpg number > 0
        if bin_cpg != 0:
            # sums the cpg number in neighboring bins with a total of input desired size
            cpg_num = 0
            # flag is used to mark if exception is raised
            flag = 0
            for i in range(coord, coord+size, 100):
                # if the start bin we pick happens to end at a corsiv region and results in key error, we try again
                try:
                    cpg_num += bin_dict[i]
                except KeyError:
                    flag = 1
                    break
           
            # return if we find one
            if ((cpg_num == cpg) and (not flag) and (genic == 'genic')):
                # bedtools command to check whether there are any genes overlapping 3kb upstream or downstream region
                # including control region itself
                end_coordinate = (coord + size - 1) 
                data = {'chr': [chromosome],'start': [coord - 3000],'end': [end_coordinate + 3000]}
                # Create DataFrame
                temp_df = pd.DataFrame(data)
                control_region_file_name = 'chrom_21_control/control_region_for_%s_%s_%s.bed' % (chromosome, coord, end_coordinate)
                gene_overlap_file_name = 'chrom_21_control/control_region_for_%s_%s_%s_gene_overlap.bed' % (chromosome, coord, end_coordinate)
                temp_df.to_csv(control_region_file_name , sep='\t', header=False, index=False)
                print("SAVED!")
                time.sleep(3)
                os.system('bedtools intersect -a %s -b 100bp_annotations/every_known_gene.bed -wb -wa > %s' % (control_region_file_name, gene_overlap_file_name))
#                 !bedtools intersect -a {control_region_file_name} -b corsiv_100kb_split1000/every_known_gene.bed -wb -wa > {gene_overlap_file_name}
                time.sleep(3)
                try:
                    gene_overlap_df = pd.read_csv(gene_overlap_file_name, header=None, sep="\t")
                    done = True
                    return coord, coord + size - 1
                except EmptyDataError:
                    continue
            elif (cpg_num == cpg and (not flag) and (not genic)):
                    done = True
                    return coord, coord + size - 1
            # otherwise try again
            else:
                count += 1
                continue
        else:
            continue
    # return 0 indicates a failure to find a control region of the desired size after 200 attempts
    return 0


# reads in the annotation file with columns "Bin Name" and "CpG Count"
annotation = pd.read_csv(annotation_filename, index_col=0)

# this dictionary stores the number of cpg in the 100bp bin as value and uses the start coord of the bin as key
bin_cpg_num = {}

# fills up the dictionary
for idx, row in annotation.iterrows():
    coord = int(row['Bin Name'].split('_')[1]) - 99
    bin_cpg_num[coord] = row['CpG Count']

print(bin_cpg_num)
# reads in CoRSIV file
# columns = [chr, start, end, cpg_count, id, ..] format
# should be 1-based coordinates
corsiv_df = pd.read_csv(corsiv_filename)
corsiv_df = corsiv_df.loc[corsiv_df['chr'] == chromosome]
# get each row corresponds to the first bin in each CoRSIV
# corsiv_df['first_bin_end'] = corsiv_df['start'] + 99
# corsiv_first_bin = corsiv_df[['chr','start', 'first_bin_end', 'cpg_count', 'id']]

# stores in the form of (start coord, end coord, cpg num) for each CoRSIV
corsiv_cpg_num = []

# fills up the list corsiv_cpg_num
for idx, row in corsiv_df.iterrows():
    start = int(row['start'])
    end = int(row['end'])
    cpg_num = int(row['cpg_count'])
    id_ = row['id']
    corsiv_cpg_num.append((start, end, cpg_num, id_))

# delete known corsivs from the dictionary
# for block in corsiv_cpg_num:
#     for num in range(block[0], block[1], 100):
#         del bin_cpg_num[num]


block_id = []
start_coord = []
end_coord = []
exact_match = []

for block in corsiv_cpg_num:
    # attempt 1-200, match = 1 means an exact match and match = 0 means an approximate
    match = 1
    result = find_match(block[2], block[1]-block[0]+1, bin_cpg_num)
    # attempt 201-400
    if not result:
        result = find_match(block[2], block[1]-block[0]+101, bin_cpg_num)
        match = 0
    # attempt 401-600
    if not result:
        result = find_match(block[2], block[1]-block[0]+201, bin_cpg_num)
        match = 0
    print(result)
    block_id.append(block[3])
    start_coord.append(result[0])
    end_coord.append(result[1])
    exact_match.append(match)
# outputs the dataframe
output = pd.DataFrame({"CoRSIV ID":block_id,"Start Coord":start_coord,"End Coord":end_coord,"Exact Match?":exact_match})
output.to_csv(output_filename, header = True, index=0)

