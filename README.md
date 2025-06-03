# Differential-Expression-Analysis
A tutorial guiding the use through differential expression analysis on RNA-Seq data

## Input Data:
- Raw Reads from sequencing, typically in fastq or fastq.gz format
- Reference genome in fasta format with GTF and GFF files
  - Ex: https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6/
  - Easiest way to get this is to navigate to the Genome page on NCBI, click FTP, and then
    copy the link for *_genomic.fna.gz, then run ```wget your_link``` on the server. Repeat
    this process with the *_genomic.gtf.gz and *_genomic.gff.gz files.

## Creating a Conda Environment:
- This step is necessary to install all the tools and scripts for this tutorial. The below code
  contains the steps for the R version of this pipeline. The tutorial assumes you run all scripts
  from this directory.
```
conda create --name rnaseq --yes
conda activate rnaseq
conda install trimmomatic htseq hisat2 bioconductor-deseq2 bioconductor-apeglm r-ggplot2 r-rcolorbrewer r-pheatmap r-base --yes r-data.table --yes
git clone https://github.com/LawrenceGordon/Differential-Expression-Analysis.git
chmod -R 755 Differential-Expression-Analysis/
cd Differential-Expression-Analysis/
```

## Typical RNA-Seq Workflow:

### 1. Trim raw reads with Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic): 

```
mkdir trimmed-reads

trimmomatic PE -threads 24 $forward $reverse\ 
    trimmed-reads/$paired_forward trimmed-reads/$unpaired_forward\
    trimmed-reads/$paired_reverse trimmed-reads/$unpaired_reverse\
    ILLUMINACLIP:$ADAPTERS:2:30:10:8:true\
    LEADING:3 TRAILING:3\
    SLIDINGWINDOW:4:10 MINLEN:36
```
This will perform trimming on one set of paired end reads.
Replace `$forward` with the R1 reads and `$reverse` with the R2 reads for each sample.
Replace `$paired_forward`, `$unpaired_forward`, `$paired_reverse`, and `$unpaired_reverse` with the desired names.
`$ADAPTERS` must be a path to a fasta file with the sequencing adapters used.


### 2. Build index for reference genome with hisat2-build (https://daehwankimlab.github.io/hisat2/manual/):

```
hisat2-build <reference_in> <ht2_base>

hisat2-build GCF_000002985.6_WBcel235_genomic.fna GCF_000002985.6_WBcel235_genomic.fna.idx
```

This will create several index files for your reference.
```
GCF_000002985.6_WBcel235_genomic.fna.idx.1.ht2
GCF_000002985.6_WBcel235_genomic.fna.idx.2.ht2
...
GCF_000002985.6_WBcel235_genomic.fna.idx.8.ht2
```

### 3. Map paired and trimmed reads to reference genome using hisat2:
The max intron length can be determined using the script `intron-length.awk`. This script is called with the gff as follows:
```
./intron-length.awk TYPE=CDS GCF_000002985.6_WBcel235_genomic.gff
```
The script will return the following info to standard out:
`MINIMUM_INTRON_LENGTH MAXIMUM_INTRON_LENGTH MAXIMUM_SUM_INTRON_LENGTHS`

```
hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]

hisat2 -x GCF_000002985.6_WBcel235_genomic.fna.idx -1 $paired_forward -2 $paired_reverse -S $sam_name -p 24 --max-intronlen MAXIMUM_INTRON_LENGTH
```
This will map one set of paired and trimmed reads to the reference genome.
Replace `$paired_forward` with the trimmed paired R1 reads and `$paired_reverse` with the trimmed paired R2 reads for each sample.
`$sam_name` refers to the output name of the alignment file, this file will be in the SAM format (https://samtools.github.io/hts-specs/SAMv1.pdf)

### 4. Counting reads in features with htseq-count:
```
htseq-count [options] <alignment_files> <gff_file>

htseq-count -f sam $sam_name GCF_000002985.6_WBcel235_genomic.gtf -t gene > $sam_name.features.txt
```
This will count the reads for one sample per gene.
Replace `$sam_name` with the SAM file for each sample.

### 5. Build count matrix using bash commands:
The output `sample.sam.features.txt` will be a table with two columns. The first column will be gene names, and the second will be respective counts for the sample. Other aligners can be used as long as the output is tab delimited with the gene names in the first column and the counts in the second column.

```
$ less -S sample.sam.features.txt

ATP6	2
CELE_2L52.1	22
CELE_2L52.2	6
CELE_2RSSE.1	4
CELE_2RSSE.3	15
...
__no_feature	13501027
__ambiguous	70624
__too_low_aQual	370657
__not_aligned	206645
__alignment_not_unique	1111979
```

To combine these individual counts files into a single matrix we can use the bash script `htseq_count_to_counts_matrix.sh`:
```
./htseq_count_to_counts_matrix.sh sample1.sam.features.txt sample2.sam.features.txt ... sampleX.sam.features.txt
```

The resulting file gene.counts.matrix can be used in DESeq2. It should have one column with gene names, and one column per sample with counts.

View this file by typing `less -S gene.counts.matrix`

```
	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8
ATP6	2	0	1	4	5	4	1	3
CELE_2L52.1	1	0	0	3	0	5	0	2
CELE_2L52.2	0	0	0	0	0	0	0	0
CELE_2RSSE.1	0	1	3	1	0	3	9	5
CELE_2RSSE.3	0	0	0	0	0	0	0	0
...
__no_feature	13501027	15563129	15994076	14471195	15563966	16712179	13239092	14673742
__ambiguous	70624	68283	80613	62025	68120	77720	60419	61401
__too_low_aQual	370657	420911	475688	374362	461920	435703	349582	390528
__not_aligned	206645	235034	224522	211690	235723	227077	192822	205528
__alignment_not_unique	1111979	911511	1157079	775521	1041590	847218	709790	715698
```

### 6.1. Running DESeq2 (Python Implementation):
The script `unfiltered-deseq.py` can be used to run DESeq2 in Python. This requires the pydeseq2 library from https://github.com/owkin/PyDESeq2.
For visualization, bioinfokit is used (https://github.com/reneshbedre/bioinfokit)
Both libraries can be installed using pip.

Once all libraries are installed, test the command using `./unfiltered-deseq.py -h`. This will bring up the help menu.
```
usage: unfiltered-deseq.py [-h] [--alpha ALPHA] [--no_cooks] [--no_independent] [--plot] [--save] [--abnormal ABNORMAL [ABNORMAL ...]] [--normal NORMAL [NORMAL ...]]
                           [--out OUT]
                           matrix

Performs DESeq2 in python

positional arguments:
  matrix                Read counts matrix generated by kallisto

options:
  -h, --help            show this help message and exit
  --alpha ALPHA, -p ALPHA
                        P-value and adjusted p-value significance threshold
  --no_cooks            whether to not filter and refit cooks outliers
  --no_independent      whether to not independent filter
  --plot                whether to plot final data
  --save                whether to store intermediates in pkl format
  --abnormal ABNORMAL [ABNORMAL ...]
                        sample prefix for abnormal samples
  --normal NORMAL [NORMAL ...]
                        sample prefix for normal samples
  --out OUT, -o OUT     output prefix for all files
```
Different options can be changed depending on desired alpha value, and whether cooks refitting, independent filtering, and plotting data are desired.
In addition you must specify the abnormal and normal sample prefixes.

If we want to run DESeq2 with alpha of 0.05, no cooks refitting, samples 1-4 as normal and samples 5-8 as abnormal, we can run the following:
```
./unfiltered-deseq.py gene.counts.matrix --alpha 0.05 --no_cooks --plot --normal sample1 sample2 sample3 sample4 --abnormal sample5 sample6 sample7 sample8 --out deseq
```
The results from the DESeq2 run will be stored in multiple files. These correspond to the raw results, the non NaN results, the LFC shrunk results,
and the isolated upregulated and downregulated gene csvs. By default, these upregulated and downregulated genes have an adjusted p-value < alpha
and log2FoldChange > 1 or < -1 respectively.

Of most importance are the files with significant genes, which will be at `{output_prefix}_upregulated_shrunk_non_nan_unfiltered_deseq_cooks_refit.csv` and `{output_prefix}_downregulated_shrunk_non_nan_unfiltered_deseq_cooks_refit.csv`.

Take a look at these output files:
```
$ less -S {output_prefix}_upregulated_shrunk_non_nan_unfiltered_deseq_cooks_refit.csv

,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
CELE_B0564.2,151.88632115624983,2.2893446679057856,0.4500016752259603,5.561048038462878,2.6815936391943935e-08,0.00022406055652288755
CELE_C23H3.4,241.7169722395296,1.9885214216629015,0.5258183740515873,4.40508314073915,1.057433135709663e-05,0.022088456413555224
CELE_C32F10.2,448.69117281224624,1.4029473790901585,0.3655887965173853,4.455525895228331,8.368779311591284e-06,0.019978667296571707
CELE_F32A5.3,48.64588240406012,1.3517904699454375,0.3893730838546,4.140799039298287,3.460980822414412e-05,0.0444895773256671
CELE_H01G02.3,48.12283560868961,1.2947949334902555,0.3589221400094513,4.246205050651833,2.1742166766293018e-05,0.03314852255624289
```

### 6.2. Running DESeq2 (R Alternative):
The script `deseq.r` can be used to run DESeq2 in R. The R version will produce PCA plots and heatmaps for the data.
The multiple test p-value adjustment method in R is different from the Python method, producing slightly different (less strict) results.
To run the script, the libraries DESeq2, pheatmap, RColorBrewer, apeglm, ggplot2, and data.table must be installed.
```
Rscript ./deseq.r <counts tsv> <affected sample prefix> <alpha value>
```
The script uses positional arguments including the path to the `gene.counts.matrix` file, the affected sample condition prefix, and the alpha value. Multiple files are produced including `raw_results.csv`, `deseq_results_non_nan.csv`, `shrunk_deseq_results_non_nan.csv`, `upregulated.csv`, `downregulated.csv`, `shrunk_upregulated.csv`, and `shrunk_downregulated.csv`. The upregulated and downregulated results are filtered on the provided alpha value. In addition, a heatmap with sample cladogram is saved to `heatmap.pdf` while a PCA plot is saved to `PCA.png`.

### 7. Summarizing DESeq2 results:
```
$ ./summarize_deseq_data_with_gff.py -h

usage: summarize_deseq_data_with_gff.py [-h] [--gff GFF] [--counts COUNTS] [--out OUT] deseq

Filters gff file based on filtered fasta file

positional arguments:
  deseq              csv file in standard deseq format

options:
  -h, --help         show this help message and exit
  --gff GFF          gff file containing functional annotations
  --counts COUNTS    counts matrix containing read counts
  --out OUT, -o OUT  prefix for output files
```
Providing a GFF file, `gene.counts.matrix`, and the deseq results to the script `summarize_deseq_data_with_gff.py` will pair the GFF annotation information with the features from the counts matrix to provide context for significant features.
```
./summarize_deseq_data_with_gff.py --gff GCF_000002985.6_WBcel235_genomic.gff --counts gene.counts.matrix --out upregulated_summary deseq_upregulated_shrunk_non_nan_unfiltered_deseq_cooks_refit.csv
```

```
$ less -S upregulated_summary.csv

feature,annotation,EV_1,EV_2,EV_3,EV_4,R10_1,R10_2,R10_3,R10_4,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
CELE_B0564.2,Fe2OG dioxygenase domain-containing protein,42,56,44,47,78,473,284,240,151.88632115624983,2.2893446679057856,0.4500016752259603,5.561048038462878,2.6815936391943935e-08,0.00022406055652288755
CELE_C23H3.4,Serine palmitoyltransferase 1,52,95,102,98,77,724,455,415,241.7169722395296,1.9885214216629015,0.5258183740515873,4.40508314073915,1.057433135709663e-05,0.022088456413555224
CELE_C32F10.2,Retinoblastoma-like protein homolog lin-35,209,209,281,215,258,952,706,857,448.69117281224624,1.4029473790901585,0.3655887965173853,4.455525895228331,8.368779311591284e-06,0.019978667296571707
CELE_F32A5.3,Serine carboxypeptidase ctsa-3.1,38,18,28,14,66,88,59,82,48.64588240406012,1.3517904699454375,0.3893730838546,4.140799039298287,3.460980822414412e-05,0.0444895773256671
CELE_H01G02.3,Uncharacterized protein,27,33,19,23,51,107,77,55,48.12283560868961,1.2947949334902555,0.3589221400094513,4.246205050651833,2.1742166766293018e-05,0.03314852255624289
```

## PALADIN Workflow:

### 1. Run raw reads through PALADIN (https://github.com/ToniWestbrook/paladin): 

```
for forward in raw-reads/*R1*
do
	output=$(basename $forward | cut -f 1 -d'.')
	paladin align -o $output -t 24 -P $proxy $database $forward
done
```

This will perform PALADIN alignment on all forward reads, putting results in the same prefix as the read.
Replace `$proxy` with a proxy address if working on an offline node and `$database` with the path to the desired database (Ex: UniRef90).

### 2. Combine PALADIN .tsv reports into a single counts matrix:
```
$ ./paladin_tsv_to_count_matrix.py -h

usage: paladin_tsv_to_count_matrix.py [-h] [--minqual MINQUAL] [--out OUT] pal [pal ...]

Combines paladin TSVs into count matrix

positional arguments:
  pal                   paladin tsvs

optional arguments:
  -h, --help            show this help message and exit
  --minqual MINQUAL, -q MINQUAL
                        minimum average quality to filter by
  --out OUT, -o OUT     output name for final file
```
This will combine all PALADIN .tsv results into a single count matrix in R DESeq2 format. An optional `--minqual` flag can be used to set a minimum quality for kept counts.
```
./paladin_tsv_to_count_matrix.py *.tsv --out paladin_count_matrix.tsv
```

```
$ less -S paladin_count_matrix.tsv

	sample1	sample2	sample3	sample4	sample5	sample6	
16ABF_BIFLN	0	0	0	0	0	1	
3AHDP_RUMGV	0	0	0	0	0	0	...
3AHD_EGGLE	0	1	2	0	0	0
...
```

### 3. Running DESeq2 (R Alternative):
The script `deseq.r` can be used to run DESeq2 in R. The R version will produce PCA plots and heatmaps for the data.
The multiple test p-value adjustment method in R is different from the Python method, producing slightly different (less strict) results.
To run the script, the libraries DESeq2, pheatmap, RColorBrewer, apeglm, ggplot2, and data.table must be installed.
```
Rscript ./deseq.r <counts tsv> <affected sample prefix> <alpha value>
```
The script uses positional arguments including the path to the `paladin_count_matrix.tsv` file, the affected sample condition prefix, and the alpha value. Multiple files are produced including `raw_results.csv`, `deseq_results_non_nan.csv`, `shrunk_deseq_results_non_nan.csv`, `upregulated.csv`, `downregulated.csv`, `shrunk_upregulated.csv`, and `shrunk_downregulated.csv`. The upregulated and downregulated results are filtered on the provided alpha value. In addition, a heatmap with sample cladogram is saved to `heatmap.pdf` while a PCA plot is saved to `PCA.png`.

### 4. Summarizing DESeq2 results:
```
$ ./summarize_deseq_data_with_paladin.py -h

usage: summarize_deseq_data_with_paladin.py [-h] [--pal PAL [PAL ...]] [--counts COUNTS] [--out OUT] deseq

Filters gff file based on filtered fasta file

positional arguments:
  deseq                csv file in standard deseq format

optional arguments:
  -h, --help           show this help message and exit
  --pal PAL [PAL ...]  paladin files containing functional annotations
  --counts COUNTS      counts matrix containing read counts
  --out OUT, -o OUT    prefix for output files
```
Providing the PALADIN tsv files, `paladin_count_matrix.tsv`, and the deseq results to the script `summarize_deseq_data_with_gff.py` will pair the PALADIN annotation information with the features from the counts matrix to provide context for significant features.
```
./summarize_deseq_data_with_paladin.py deseq_results_non_nan.csv --pal *tsv --counts paladin_count_matrix.tsv --out deseq_gene_summary
```

```
$ less -S deseq_gene_summary.csv

feature,annotation,C054,C055_t2dm,C056_t2dm,C057,C058,C059,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
A0A014AUH4_9GAMM,Acinetobacter sp. 25977_7;Uncharacterized protein;J810_4084,0,0,28,0,0,0,4.30397859596638,0.373512097764134,1.15638025333893,0.323001103387623,0.746694406435262,0.999926231467233
A0A015S3H7_BACFG,Bacteroides fragilis str. 3986 N(B)19;Alanine racemase (EC 5.1.1.1);alr M085_0754,0,7,17,0,0,0,3.32151320511241,0.290512098421384,1.00191425068489,0.28995704794377,0.771849097381583,0.999926231467233
A0A015T1R3_BACFG,Bacteroides fragilis str. 3988T(B)14;Conjugative transposon TraM protein;M124_4633,0,0,3,0,0,0,3.21828923466149,-1.26343332614334,1.33070141327592,-0.949449150304144,0.342392222559533,0.999926231467233

```
