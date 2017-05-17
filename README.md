# Practical Examples

Companion repository for TBD paper.

## Prerequisites

### Software

- [DIAMOND 0.8.36][diamond]
- [MetaPhlAn2][metaphlan2]
- various utilities required by the [Microbiome Helper][microbiome_helper]

### Environment

In addition to having the prerequisite software installed, the steps below
assume a Ubuntu environment (though macOS will likely work as well) and
Python 3.6 installed.

To install the additional Python dependencies, [Anaconda/Miniconda][anaconda]
is recommended. Assuming Miniconda is installed, run the following from the
project root:

```
$ conda env create --file environment.yml
$ source activate examples
```

All commands below assume that the `examples` Python environment is active,
unless otherwise noted.

### Databases

Two [Universal Protein Resource (UniProt)][uniprot] databases are required (see
the [UniProt downloads page][uniprot_downloads]). This project uses the 2017_04
release, though newer releases may give similar results:
  - [UniProt Reference Clusters at 100% identity (UniRef100)][uniref100_ftp]
    - both XML and FASTA versions required
  - [UniProt Knowledgebase (UniProtKB) ID Mapping][uniprotkb_ftp]

To download the latest version of the databases (older releases can be
downloaded from the [UniProt FTP site][previous_rel_ftp]):

```
$ mkdir uniref
$ curl -o uniref/uniref100.fasta.gz \
    ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
$ curl -o uniref/idmapping.dat.gz \
    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz    
```

The SHA1 checksums for the databases are:

```
85ebde5758f2067175a4ac29a8245a495b7ce1f3  uniref/uniref100.fasta.gz
35d862ea431f6271a8f92d70e983684c9d02080d  uniref/uniref100.xml.gz
ba34ed44c27d63b96d1d48f4f914b57d826eec48  uniref/idmapping.dat.gz
```

Once the downloads are complete, a version of the UniRef100 appropriate for use
with DIAMOND must be generated:

```
$ diamond makedb --threads 24 \
    --in uniref/uniref100.fasta.gz \
    --db uniref/uniref100 \
    | tee uniref/makedb.uniref100.log
```

In addition, we generate some new mapping files for downstream analyses. One
to map UniRef100 sequences to UniRef90 and UniRef50:

```
$ python scripts/make_uniref_mapping.py \
    uniref/uniref100.xml.gz uniref/uniref_mapping.tsv
```

And another to map UniProtKB accession numbers to NCBI taxa IDs, 
UniRef Clusters and KEGG Orthologs, Modules and Pathways:

```
$ python scripts/make_uniprotkb_mapping.py \
    uniref/idmapping.dat.gz \
    kegg_id_mapping.tsv \
    --output-mapping uniref/uniprot_to_other.tsv
```

Dataset Processing
------------------

Raw data will be stored in `data/`:

```
$ mkdir data
```

### Download data

HMP FASTQ data was downloaded from the NCBI Sequence Read Archive (SRA) using
their [SRA tools][sra_tools].

To get started quickly, simply run:

```
$ bash scripts/prefetch_cmds.sh
```

The `metadata/` folder contains more details about the samples:
- `metadata/HMASM.csv` contains the raw info from the HMP website for
  all of the FASTQs
- `metadata/HMASM_test_subset.csv` is a subset of 10 stool and 10
  tongue dorsum samples (both randomly selected).
- the commands in `SRS_ids_convert.R` (in the `scripts/` directory in the
  in the project root) were used to convert sample ids (SRS) to run ids (SRR)
  (corresponding to the downloaded FASTQs). This mapping was saved to
  `metadata/SRS_SRR_ids_linked.txt`.
- Prefetch (a SRA tools command) can be used to download all the files
  (automatically put in your user folder, e.g. `/home/user1/ncbi/public/sra)`).
  To create a shell script to run prefetch in fewer lines:

      $ tail -n +2 SRS_SRR_ids_linked.txt | awk '{ print "prefetch "$2 }' \
        >> prefetch_cmds.sh

### Prepare sample sequences and QC

First, convert from `.sra` to FASTQ:

```
$ mkdir -p data/raw_fastqs
$ fastq-dump \
    --gzip \
    --skip-technical \
    --readids \
    --dumpbase \
    --split-files \
    --clip \
    --outdir data/raw_fastqs \
    ~/ncbi/public/sra/*.sra
```

Concatenate reads for each sample:

```
$ mkdir -p data/concatenated_fastq/
$ python scripts/cat_sample_reads.py data/raw_fastqs data/concatenated_fastq
```

Next, filter the reads:

```
$ run_trimmomatic.pl \
    --leading_quality 5 \
    --trailing_quality 5 \
    --required_quality 15 \
    --window_size 4 \
    --min_length 70 \
    --jar /usr/local/prg/Trimmomatic-0.36/trimmomatic-0.36.jar \
    --thread 20 \
    --out_dir data/trimmomatic_filtered \
    data/concatenated_fastq/*.fastq.gz
```

Concatenate the paired ends together:

```
$ concat_paired_end.pl -p 4 -o data/filtered_reads \
    data/trimmomatic_filtered/*_paired*fastq*
```

Convert FASTQ files to FASTA:

```
$ parallel --jobs 2 \
  'zcat data/filtered_reads/{} \
    | fastq_to_fasta -v -n -z -o data/fasta/{/.}.fasta.gz \
    > fasta/{/.}.log' \
  ::: data/filtered_reads/*fastq.gz
$ rename s/.fastq// data/fasta/*
```

After conversion, each sample will have 10s to 100s of millions of reads.
Subsample to ~7M reads to make downstream analysis faster:

```
$ mkdir data/fasta-subsampled-7M
$ parallel --jobs 2 \
  'export SEED=$RANDOM; \
    echo "{} $SEED" >> data/fasta-subsampled-7M/seeds.txt; \
    seqtk sample -s $SEED {} 7000000 > data/fasta-subsampled-7M/{/.}' \
  ::: data/fasta/*fasta*
$ gzip data/fasta-subsampled-7M/*.fasta
```

## Analysis

### Functional composition

Run DIAMOND search against the UniRef100 database. Importantly, we configure
the search to only return the top hit for each query at an identity cutoff of
90%. This is a time consuming step, taking approximately 12 hours on a server
with 48 vCPUs:

```
$ parallel --jobs 2 \
  'diamond blastx \
      --db uniref/uniref100.2017_04.dmnd \
      --query {} \
      --max-target-seqs 1 \
      --id 0.9 \
      --out data/diamond_out/{/.}.txt \
      --block-size 6 \
      --threads 20 \
      --index-chunks 1 \
    > data/diamond_out/{/.}.log' \
    ::: data/fasta-subsampled-7M/*fasta*
```

Now, generate tables of alignment counts and scale down to other databases:

```
# No mapping (UniRef100)
$ python scripts/aggregate_alignments.py \
    --output-file tables/uniref100.spf \
    --summary-file uniref100_summary.tsv \
    diamond_out/hmp/uniref100/SRS0*.txt

# UniRef 90 and 50
$ parallel --jobs 2 --link \
    'python scripts/aggregate_alignments.py \
      --mapping-file uniref/uniref_mapping.tsv \
      --from-type UniRef100 \
      --to-type {1} \
      --output-file tables/uniref100_to_{2}.spf \
      --summary-file uniref100_to_{2}_summary.tsv \
      diamond_out/hmp/uniref100/SRS0*.txt' \
    ::: UniRef90 UniRef50 ::: uniref90 uniref50

# KOs, Modules and Pathways
$ parallel --jobs 3 --link \
    'python scripts/aggregate_alignments.py \
      --mapping-file uniref/uniprotkb_to_other.tsv \
      --from-type UniRef100 \
      --to-type {1} \
      --output-file tables/uniref100_to_{2}.spf \
      --summary-file uniref100_to_{2}_summary.tsv \
      diamond_out/hmp/uniref100/SRS0*.txt' \
    ::: KO "KEGG Modules" "KEGG Pathways" ::: kos modules pathways
```

The previous commands will have generated summary files that we will now
combine for convenience:

```
$ for f in *_summary.tsv; do tail -n+2 $f | grep Alignments \
    >> alignments_summary.tsv; done
$ for f in *_summary.tsv; do tail -n+2 $f | grep Subjects \
    >> subjects_summary.tsv; done
$ (head uniref100_summary.tsv; cat alignments_summary.tsv) > \
    tables/functional_alignments_summary.tsv
$ (head uniref100_summary.tsv; cat subjects_summary.tsv) > \
    tables/functional_subjects_summary.tsv
$ rm *_summary.tsv
```

### Taxonomic composition (Using MetaPhlAn2)

Run MetaPhlAn2:

```
$ run_metaphlan2.pl -p 20 -o metaphlan-taxonomy.txt
$ metaphlan_to_stamp.pl metaphlan-taxonomy.txt > tables/metaphlan-taxonomy.spf
$ rm -rf metaphlan_out; rm metaphlan-taxonomy.txt
```

### Generating figures

To generate the figures, run the following R scripts:

```
$ mkdir figures
$ Rscript scripts/num_taxa_per_func_distribution.R
$ Rscript scripts/function_taxonomy_box_plots.R
```

## Citing and Contributors

Citation TBD

Contributors (to this repository):
  - Gavin Douglas
  - Karl Leuschen


[anaconda]: https://www.continuum.io/anaconda-overview
[previous_rel_ftp]: ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/
[uniprot_downloads]: http://www.uniprot.org/downloads
[uniprotkb_ftp]: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
[uniref100_ftp]: ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/
