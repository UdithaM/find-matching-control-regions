# find-matching-control-regions
##### Finds matching control regions to a set of given CoRSIV regions, in terms of CpG density, region size and in the vicinity of genes (optional).

&nbsp;
&nbsp;

## Steps to use

#### 1)  Clone this Github repository to your Linux computer/Server.
#### 2) Go into the folder.

```sh
cd find-matching-control-regions
```

### 3) Make sure you get the correct paths and values as input arguments. You can type following command to find more about the input arguments.

```sh
python get_control_regions.py --help
```
#### Output 

```
--------Get Control Regions For CoRSIV Regions--------

positional arguments:
  Annotation_File       Path to annotation file (CSV) for a given chromosome. Columns of annotation file should be Bin Name,CpG Count, ... and should be in 1-based coordinates. Values in the column 'Bin Name' should be in 'chrX_<ENDING_COORDINATE_OF_BIN>' format. (ex: chr21_4300)
  Annotation_File_Bin_Size
                        Bin Size of the annotation file provided. ( ex: 100, 200)
  CoRSIV_File           Path to CoRSIV file (CSV) for a given chromosome. Columns of CoRSIV file should be chr, start, end, cpg_count, id, ... and should be in 1-based coordinates
  Chromosome            Chromosome in the format chrx.
  Output_File           Path of the output file containing control regions.

optional arguments:
  -h, --help            show this help message and exit
  -g, --genic           Whether genic control regions are needed. Make sure that gene file is provided with if this is set to True. Use -f or --gene_file to provide gene file. (default: False)
  -f GENE_FILE, --gene_file GENE_FILE
                        Path to file containing genes in whole genome in BED format. Columns of gene file should be file chr, start, end, ... and should be in 1-based coordinates. Column headers should not be in the bed file.
```

## Note : 

#### In the column "Bin Name" in the annotation file, values should be in the following format.

```sh
chrX_<ENDING_COORDINATE_OF_BIN>
```
#### Example: 
```sh
chr21_4300
```


### 4)  Use the following "sample" command to run the Python script to get control regions (without genic option).
```sh
python get_control_regions.py 100bp_annotations/chr21_annotated.csv 100 corsiv_regions.csv chr21 chrom_21_control_regions.csv
```

### 5)  Use the following "sample" command to run the Python script to get control regions with genic option (Optional).
```sh
python get_control_regions.py -g -f every_known_gene.bed 100bp_annotations/chr21_annotated.csv 100 corsiv_regions.csv  chr21 chrom_21_control_regions.csv
```

## License

MIT

