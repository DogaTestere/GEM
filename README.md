# iumobg/model_creation


[![GitHub Actions CI Status](https://github.com/iumobg/model_creation/actions/workflows/nf-test.yml/badge.svg)](https://github.com/iumobg/model_creation/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/iumobg/model_creation/actions/workflows/linting.yml/badge.svg)](https://github.com/iumobg/model_creation/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/iumobg/model_creation)

## Introduction

**iumobg/model_creation** is a bioinformatics pipeline that converts raw FASTQ sequencing reads into a genome-scale metabolic model (GEM) and performs constraint-based simulations. The full pipeline conducts quality control, genome assembly(de novo or reference-guided) and annotation(ab-inito or reference-based) then uses chosen database(currently BiGG, KEGG) to construct a GEM.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->

1. Quality Control via [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Genome Assembly
2.a. De Novo Assembly via [`SPAdes`](https://github.com/ablab/spades)
2.b. Reference-guided Assembly via [`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/), [`SAMtools`](https://www.htslib.org/) and [`BCFtools`](https://samtools.github.io/bcftools/)
3. Genome Annotation
3.a.1. Ab inito annotation for prokaryotes via Prodigal
3.a.2. Database guided annotation for prokaryotes via miniprot
3.b.1. Ab inito annotation for eukaryotes via GeneMark-ES
3.b.2. Database guided annotation for eukaryotes via miniprot
4. Database enrichment via python scripts : Currently defaults to KEGG, has options for BiGG and NCBI
5. Model creation via CobraPy

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. 

First, prepare a samplesheet with your input data that looks as follows:
`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,has_ref,ref_file,genome_type,ann_type,db_type
SAMPLE1, ec1_S1_L001_R1_001.fastq.gz, ec1_S1_L001_R2_001.fastq.gz,true, e_coli_ref.fasta, pro, ref_in, KEGG
```
Each row represents one biological sample.
The columns are defined as follows:
- `sample`: Unique sample identifier
- `fastq_1`: FASTQ file for read 1 (File path)
- `fastq_2`: FASTQ file for read 2 (File path)
- `has_ref`: Whether the FASTQ file have a reference to use for assembly (Boolean)
- `ref_file` : Reference file to be used in Reference-guided assembly (File Path)
- `genome_type` : Identifier used for prokaryotes(pro) and eukaroytes(euk) (String)
- `ann_type` : Whether ab initio(ab_in) or database guided annotation(ref_in) should be used (String)
- `pep_accession` : RefSeq accession number that would be used for downloading the protein information from NCBI when guided annotation is used. (String)
- `pep_file` : Protein information file path when guided annotation is used. (File Path)
- `db_type` : Name of the database that should be used to mainly search for information in order to build the GEM. Defaults to KEGG, has 'BiGG' and 'NCBI' options. (String)

Now, you can run the pipeline using:

```bash
nextflow run iumobg/model_creation \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

iumobg/model_creation was originally written by Doga Yasemen Testere, İrem Ay.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use iumobg/model_creation for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
