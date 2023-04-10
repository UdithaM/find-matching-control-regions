import pandas as pd
import random
import sys
import os
import time
# from pandas.io.common import EmptyDataError
import argparse
 
parser = argparse.ArgumentParser(description="--------Get Control Regions For CoRSIV Regions--------",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-g", "--genic", action="store_true", help="Whether genic control regions are needed. Make sure that gene file is provided with if this is set to True. Use -f or --gene_file to provide gene file.")
parser.add_argument("-f", "--gene_file",  help="Path to file containing genes in whole genome in BED format. Columns of gene file should be file chr, start, end, ... and should be in 1-based coordinates. Column headers should not be in the bed file.")
parser.add_argument("Annotation_File", help="Path to annotation file (CSV) for a given chromosome. Columns of annotation file should be Bin Name,CpG Count, ... and should be in 1-based coordinates. Values in the column 'Bin Name' should be in 'chrX_<ENDING_COORDINATE_OF_BIN>' format. (ex: chr21_4300)")
parser.add_argument("Annotation_File_Bin_Size", help="Bin Size of the annotation file provided. ( ex: 100, 200)")
parser.add_argument("CoRSIV_File", help="Path to CoRSIV file (CSV) for a given chromosome. Columns of CoRSIV file should be chr, start, end, cpg_count, id, ... and should be in 1-based coordinates")
parser.add_argument("Chromosome", help="Chromosome in the format chrx.")
parser.add_argument("Output_File", help="Path of the output file containing control regions.")
args = parser.parse_args()
config = vars(args)
if (config['genic'] and config['gene_file'] is None):
    print("ERROR: Please enter the path to file containing genes in whole genome in BED format. Columns of gene file should be file chr, start, end, ... and should be in 1-based coordinates. Column headers should not be in the bed file. ")



# parameters
annotation_filename = config['Annotation_File'] 
annotation_bin_size = int(config['Annotation_File_Bin_Size'])
corsiv_filename = config['CoRSIV_File']
output_filename = config['Output_File']
chromosome = config['Chromosome']
genic = config['genic']
gene_filename = config['gene_file']

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
    while not done and count <= 400:
        # randomly pick a bin
        coord, bin_cpg = random.choice(list(bin_dict.items()))
        # start at a bin with cpg number > 0
        if bin_cpg != 0:
            # sums the cpg number in neighboring bins with a total of input desired size
            cpg_num = 0
            # flag is used to mark if exception is raised
            flag = 0
            for i in range(coord, coord+size, annotation_bin_size ):
                # if the start bin we pick happens to end at a corsiv region and results in key error, we try again
                try:
                    cpg_num += bin_dict[i]
                except KeyError:
                    flag = 1
                    break
           
            # return if we find one
            if ((cpg_num == cpg) and (not flag) and genic):
                # bedtools command to check whether there are any genes overlapping 3kb upstream or downstream region
                # including control region itself
                end_coordinate = (coord + size - 1) 
                data = {'chr': [chromosome],'start': [coord - 3000],'end': [end_coordinate + 3000]}
                # Create a dataFrame for bedtools intersect including control region
                temp_df = pd.DataFrame(data)
                control_folder = '%s_controls_temp'  % chromosome
                # Create a folder for bedtools
                os.system('mkdir -p %s' % control_folder)
                time.sleep(2)
                control_region_file_name = '%s/control_region_for_%s_%s_%s.bed' % (control_folder,chromosome, coord, end_coordinate)
                gene_overlap_file_name = '%s/control_region_for_%s_%s_%s_gene_overlap.bed' % (control_folder,chromosome, coord, end_coordinate)
                temp_df.to_csv(control_region_file_name , sep='\t', header=False, index=False)
                time.sleep(3)
                # Check for overlaps with genes within 3kb upstream and downstream
                os.system('bedtools intersect -a %s -b %s -wb -wa > %s' % (control_region_file_name, gene_filename, gene_overlap_file_name))
#                 !bedtools intersect -a {control_region_file_name} -b corsiv_100kb_split1000/every_known_gene.bed -wb -wa > {gene_overlap_file_name}
                time.sleep(3)
                try:
                    # If genic return matching genic control region
                    gene_overlap_df = pd.read_csv(gene_overlap_file_name, header=None, sep="\t")
                    done = True
                    return coord, coord + size - 1
                except pd.errors.EmptyDataError:
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
annotation = pd.read_csv(annotation_filename)

# this dictionary stores the number of cpg in the 100bp bin as value and uses the start coord of the bin as key
bin_cpg_num = {}

# fills up the dictionary
for idx, row in annotation.iterrows():
    coord = int(row['Bin Name'].split('_')[1]) - (annotation_bin_size -1)
    bin_cpg_num[coord] = row['CpG Count']

print(bin_cpg_num)
# reads in CoRSIV file
# columns = [chr, start, end, cpg_count, id, ..] format
# should be 1-based coordinates
corsiv_df = pd.read_csv(corsiv_filename)
corsiv_df = corsiv_df.loc[corsiv_df['chr'] == chromosome]

# stores in the form of (start coord, end coord, cpg num) for each CoRSIV
corsiv_cpg_num = []

# fills up the list corsiv_cpg_num
for idx, row in corsiv_df.iterrows():
    start = int(row['start'])
    end = int(row['end'])
    cpg_num = int(row['cpg_count'])
    id_ = row['id']
    corsiv_cpg_num.append((start, end, cpg_num, id_))


block_id = []
start_coord = []
end_coord = []
exact_match = []

failed_block_id = []

for block in corsiv_cpg_num:
    # attempt 1-400, match = 1 means an exact match and match = 0 means an approximate
    match = 1
    result = find_match(block[2], block[1]-block[0]+1, bin_cpg_num)
    # attempt 401-800
    if not result:
        result = find_match(block[2], block[1]-block[0]+101, bin_cpg_num)
        match = 0
    # attempt 801-1200
    if not result:
        result = find_match(block[2], block[1]-block[0]+201, bin_cpg_num)
        match = 0
    if not result:
        result = find_match(block[2], block[1]-block[0], bin_cpg_num -1)
        match = 0
    if not result:
        result = find_match(block[2], block[1]-block[0], bin_cpg_num +1)
        match = 0
    if not result:
        result = find_match(block[2], block[1]-block[0] + 101, bin_cpg_num +1)
        match = 0
    if not result:
        result = find_match(block[2], block[1]-block[0] + 201, bin_cpg_num +1)
        match = 0
    if not result:
        result = find_match(block[2], block[1]-block[0], bin_cpg_num -2)
        match = 0
    if not result:
        result = find_match(block[2], block[1]-block[0], bin_cpg_num +2)
        match = 0

    print(result)
    try:
        start_coord.append(result[0])
        block_id.append(block[3])
        end_coord.append(result[1])
        exact_match.append(match)
    except TypeError:
        failed_block_id.append(block[3])
        print("Failed to find a matching control region for CoRSIV with id %s !" % (block[3]))
# outputs the dataframe
output = pd.DataFrame({"CoRSIV ID":block_id,"Start Coord":start_coord,"End Coord":end_coord,"Exact Match?":exact_match})
output.to_csv(output_filename, header = True, index=0)
failed_output = pd.DataFrame({"Failed CoRSIV ID":failed_block_id})
failed_output.to_csv("Failed_CoRSIVs_%s.csv" % chromosome, header = True)

