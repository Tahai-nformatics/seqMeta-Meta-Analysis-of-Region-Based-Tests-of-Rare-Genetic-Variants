## seqMeta Pipeline

An automatic pipeline to run seqMeta from a VCF/GEN file.

### Software Requirements
1. Python 2.7 (or Python 3)
2. Pandas 0.23.4 (or above)
3. DRMAA 0.7.9 (or above)
4. pybgen 0.7.0 (or above)
5. R 3.2.3 (or above)
6. seqMeta 1.6.7

#### Optional (recommended for faster, parallel processing)
1. bsub/LSF

### Quick install
+ Install Pandas Python module:
```
pip install pandas drmaa pybgen --user
```
Eliminate `--user` option for installing in the system directory. 

+ Install seqMeta R package:
```
install.packages(seqMeta)
```

### Features of running the pipeline
1. Start from VCF file (full pipeline)
2. Start from GEN file
3. Resume to (re-)run seqMeta analysis given all per-gene genotype files that have been already generated

### Examples
#### Batch (all chromosomes):
Parallel mode (distributed processing with bsub/LSF). Per-chromosome jobs will be distributed across bsub/LSF hosts

+ Parallel mode example:
```
bash seqmeta_pipeline_bsub.sh /PATH/OF/INPUT/VCF/DIR /PATH/OF/PHENOTYPE/FILE /PATH/OF/SNP/INFO/DIR /PATH/OF/OUTPUT/DIR PROJECT_TAG
```
```
USAGE: seqmeta_pipeline_bsub.sh <vcf_dir> <phenotype_file> <snp_info_dir> <out_dir> <project_tag>
NOTE: <vcf_dir> should contain chr<chrNum>.dose.vcf.gz per-chromosome VCF files
Output will be written to <out_dir>/chr<chrNum> subfolders for each chromosome
IMPORTANT: subject IDs in <phenotype_file> should match subject IDs in VCF files in the same order!
```

Consecutive mode (local multiprocessing will be used). Per-chromosome jobs will be executed consecutively (i.e., chr1, chr2, ...). The individual job will be parallelized using multiprocessing.

+ Sequential mode example:
```
bash seqmeta_pipeline_consecutive.sh /PATH/OF/INPUT/VCF/DIR /PATH/OF/PHENOTYPE/FILE /PATH/OF/SNP/INFO/DIR /PATH/OF/OUTPUT/DIR PROJECT_TAG
```
```
bash seqmeta_pipeline_consecutive.sh
USAGE: seqmeta_pipeline_consecutive.sh <vcf_dir> <phenotype_file> <snp_info_dir> <out_dir> <project_tag>
NOTE: <vcf_dir> should contain chr<chrNum>.dose.vcf.gz per-chromosome VCF files
Output will be written to <out_dir>/chr<chrNum> subfolders for each chromosome
IMPORTANT: subject IDs in <phenotype_file> should match subject IDs in VCF files in the same order!
```



#### Individual chromosome
```
python seqmeta_pipeline.py -p PREFIX -j PROJECT_TAG -w /PATH/OF/OUTPUT/DIR -c CPU# /PATH/OF/INPUT/VCF/FILE /PATH/OF/SNP/INFO/DIR /PATH/OF/PHENOTYPE/FILE CHR#
```
```
Usage:
Start a new analysis: seqmeta_pipeline.py [options] <VCF file/-g GEN file/-b BGEN file> <Directory of SNP info files> <phenotype file> <chr number>
Resume previous work: seqmeta_pipeline.py [options] <-r seqMeta_pipeline.log>

Version: 1.2.0

Options:
  -h, --help            show this help message and exit
  -c CPU, --cpu=CPU     How many jobs submitted simultaneously at one time 
                        [Default: cpu_count()]
  -g, --gen             Specify the input file is a GEN file instead of
                        default VCF file [Default: False]
  -b, --bgen            Specify the input file is a BGEN file instead of
                        default VCF file [Default: False]
  -p PREFIX, --prefix=PREFIX
                        Add a prefix tag in file names of all outputs
                        [Default: '']
  -r RESUME_LOG, --resume=RESUME_LOG
                        Give a seqMeta_pipeline.log file to resume previous
                        work [Default: '']
  -u, --unverify        Skip to verify sample IDs and their order between VCF
                        and phenotype files. (Refer to -v option) [Default:
                        False]
  -v VERIFY_COL_NAME, --verify=VERIFY_COL_NAME
                        Specify a column name of the phenotype file to verify
                        sample IDs and their order between VCF and phenotype
                        files. (e.g., '-v IID' indicates sample IDs in the IID
                        column) [Default: 'IID']
  -w WORKING_DIR, --workdir=WORKING_DIR
                        The working directory for input/output files of
                        seqMeta [Default: the directory of the seqmeta python 
						script]
  -j PROJECT, --project=PROJECT
                        Project job tag [Default: '']
  --rscript_bin=RSCRIPT_BIN_PATH
                        Path of executable file for Rscript [Default:
                        'Rscript']
  --version             Show the version of this script
```

#### Individual gene
```
Rscript run_seqmeta.R GENE_NAME /PATH/OF/DOSAGE/FILE /PATH/OF/SNP/INFO/FILE /PATH/OF/PHENOTYPE/FILE PREFIX /PATH/OF/OUTPUT/DIR
```

### Note
1. All genes with less than 3 variants will not be run (seqMeta restriction)
2. Phenotype file requirements:
	+ IID column: should match VCF subject IDs (same IDs, in the same order)
	+ status column: subjects with -9 status will be filtered out when running prepScores; values expected to be 1 or 2 (seqMeta requires 0<=y<=1)
	+ pc1, pc2, pc3 columns must be present
3. IMPORTANT:
   + Expected seqMeta memory usage can be computed as `n^2 x 8 x 10 x (1.5)`, where `n` is the number of SNPs. E.g., for a machine with 64GB of RAM the maximum number of SNPs will be approx.
