# lifebit-ai/relate: Usage

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run lifebit-ai/relate --input ..
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull lifebit-ai/relate
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/relate releases page](https://github.com/nf-core/relate/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Required relate pipeline arguments
### `--input`
File with list of full paths to bcf files and their indexes. Bcf files can be compressed but in a readable for bcftools format.

Example argument:
```
--input s3://lifebit-featured-datasets/projects/gel/siteqc/input.csv
```
Example file content:
```
bcf,index
s3://lifebit-featured-datasets/projects/gel/siteqc/test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz,s3://lifebit-featured-datasets/projects/gel/siteqc/test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz.csi
```
### `--inputDir`
Input dir for annotation ".txt" files for each of the regions (or chunks ), comming from the metrics compoment and first aggregate annotation of SiteQC pipeline. 

Example argument:
```
--inputDir s3://lifebit-featured-datasets/projects/gel/siteqc/ARtestFiles/Annotation_newtest/
```
Example file content of txt files inside directory:
```
chr10	52955340	A	G	PASS
chr10	52955649	T	C	PASS
chr10	52955673	T	G	PASS
chr10	52955677	A	G	PASS
chr10	52955680	A	C	PASS
chr10	52956494	C	T	PASS
chr10	52956566	G	A	PASS
chr10	52956631	A	G	PASS
chr10	52956723	T	A	PASS
chr10	52957225	G	A	PASS
chr10	52957263	C	T	PASS
chr10	52957282	A	G	PASS
chr10	52957321	C	T	PASS
chr10	52957381	A	G	PASS
chr10	52957427	G	T	PASS
```
### `--inputMichiganLDfileExclude`
File with regions to be filtered out for improving quality of sites selected.

Example argument:
```
--inputMichiganLDfileExclude s3://lifebit-featured-datasets/projects/gel/siteqc/MichiganLD_liftover_exclude_regions_PARSED.txt 
```
Example file content:
```
chr1 47534328 51534328 1
chr2 133742429 137242430 1
chr2 182135273 189135274 1
chr3 47458510 49962567 1
chr3 83450849 86950850 1
chr5 98664296 101164296 1
chr5 129664307 132664308 1
chr5 136164311 139164311 1
chr6 24999772 35032223 1
chr6 139678863 142178863 1
chr8 7142478 13142491 1
chr8 110987771 113987771 1
chr11 87789108 90766832 1
chr12 109062195 111562196 1
chr20 33412194 35912078 1
```

### `--inputPCsancestryrelated`
File with Principal Components information comming from reference resources of GEL for the inferred ancestries from the 30k dataset.

Example argument:
```
--inputPCsancestryrelated s3://lifebit-featured-datasets/projects/gel/siteqc/aggV2_bedmerge_30KSNPs_labkeyV9_08062020_update_PCsancestryrelated.tsv
```
Example file content:
```
plate_key Pc1 Pc2 Pc3 Pc4 Pc5 Pc6 Pc7 Pc8 Pc9 Pc10 AFR SAS EAS EUR AMR unrelated_set
HG002 -0.00264664 -0.001827457 -0.000782506 -0.00126352 -0.000198434 -0.000456 0.000632326 -0.000606751 -0.00293623 0.00449653 0 0 0 0.99 0.01 1
HG003 -0.00181639 -0.000681688 -0.00173264 -0.000721908 -0.0004302 -0.000634145 0.0000852272 -0.00421031 0.00599535 0.00272445 0 0 0 0.01 0.99 0
HG004 -0.00180975 -0.00269357 -0.00153047 -0.00223267 0.00132476 -0.000542491 -0.0008464999 -0.00379903 -0.0001633 0.00455183 0 0 0 1 0 0
```


### `--inputFinalPlatekeys`
File with a list of platekeys to be included in the analysis , required for "create_final_king_vcf" process.

Example argument:
```
--inputFinalPlatekeys s3://lifebit-featured-datasets/projects/gel/siteqc/sampleList.txt
```
Example file content:
```
HG002
sample_1
sample_2
sample_3
HG003
HG004
```


### `--inputProbs200K`
File with Ancestry assignments from GEL's reference resources.

Example argument:
```
--inputProbs200K s3://lifebit-featured-datasets/projects/gel/siteqc/aggV2_ancestry_assignment_probs_1KGP3_200K.tsv
```
Example file content:
```
Sample AFR SAS EAS EUR AMR ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI
HG002 0 0 0 0.89 0.01 1 0 0 0 0.01 0 0 0 0 0 0.21 0 0 0.01 0 0 0 0 0 0 0 0 0 0 0.01 0
HG003 0 0 0 0.01 0.99 0 0 0 0 0.01 0 0 0 0 0 0.24 0 0 0 0 0 0 0 0 0 0 0 0 0 0.01 0
HG004 0 0 0 1 0 0 0 0 0 0.35 0 0 0 0 0 0.41 0 0 0.03 0 0 0 0 0 0 0 0 0 0 0.16 0
```

### `--inputUNRELATED_1KGP3`
File required for the infer_ancestry process , is a two column tab separated file with platekeys on each column.

Example argument:
```
--inputUNRELATED_1KGP3 s3://lifebit-featured-datasets/projects/gel/siteqc/UNRELATED_1KGP3.samples
```
Example file content:
```
HG00096	HG00096
HG00097	HG00097
HG00099	HG00099
HG00100	HG00100
HG00101	HG00101
HG00102	HG00102
HG00103	HG00103
HG00105	HG00105
HG00106	HG00106
HG00107	HG00107
HG00108	HG00108
HG00109	HG00109
HG00110	HG00110
HG00111	HG00111
```

### `--input1KGP3`
File required for the infer_ancestry process , is a three column tab separated file with Sample(Platekey), Family ID and Population asignments.

Example argument:
```
--input1KGP3 s3://lifebit-featured-datasets/projects/gel/siteqc/1KGP3.sample_table
```
Example file content:
```
Sample	Family ID	Population
HG00096	HG00096	GBR
HG00097	HG00097	GBR
HG00098	HG00098	GBR
HG00099	HG00099	GBR
HG00100	HG00100	GBR
HG00101	HG00101	GBR
HG00102	HG00102	GBR
```

### `--inputSuper_pop_codes`
File required for the infer_ancestry process, its tab separated file with population and super population relations and related information.

Example argument:
```
--inputSuper_pop_codes s3://lifebit-featured-datasets/projects/gel/siteqc/super_pop_codes.tsv
```
Example file content:
```
Population	PopulationD	Super_Population
CHB	Han Chinese in Beijing, China	EAS
JPT	Japanese in Tokyo, Japan	EAS
CHS	Southern Han Chinese	EAS
CDX	Chinese Dai in Xishuangbanna, China	EAS
KHV	Kinh in Ho Chi Minh City, Vietnam	EAS
CEU	Utah Residents (CEPH) with Northern and Western European Ancestry	EUR
```

### `--input05both1K100K_eigenvec`
File required for the infer_ancestry process, its tab separated file with eigenvectors.

Example argument:
```
`--input05both1K100K_eigenvec s3://lifebit-featured-datasets/projects/gel/siteqc/1KGP3_30K_unrel_autosomes.eigenvec
```
Example file content:
```
HG00096 HG00096 -0.0324126 0.000243448 0.00664191 -0.0131188 0.00715375 0.0153836 0.0457208 -0.0213814 0.0183829 -0.0217617 -0.0241451 0.0562236 -0.0437782 0.00829286 0.0866145 -0.0442724 -0.0484038 0.0465044 0.0462053 -0.0131035
HG00097 HG00097 -0.0316075 0.00138343 0.00602363 -0.0169859 0.00928749 -0.00330704 0.052626 -0.0140852 0.0236413 0.00480714 -0.0178737 0.0187442 -0.0834035 -0.000704524 -0.00238823 -0.0899047 0.00142339 0.0164213 0.0256184 -0.0573523
HG00099 HG00099 -0.0320885 -0.000516738 0.0103002 -0.0112502 0.0195807 0.00096675 0.0772091 -0.0222773 0.0188928 0.0118046 0.0155667 -0.0125489 -0.0275191 -2.1464e-05 0.0182104 -0.0459963 -0.0136182 -0.0236272 0.0540468 0.0127228
```
### `--inputGELprojection_proj_eigenvec`
File required for the infer_ancestry process, its tab separated file with eigenvectors and it has "NA" values in the last column.

Example argument:
```
--inputGELprojection_proj_eigenvec s3://lifebit-featured-datasets/projects/gel/siteqc/1KGP3_30K_unrel_autosomes.eigenvec_TMPPROJ
```
Example file content:
```
HG00096 HG00096 -0.0324126 0.000243448 0.00664191 -0.0131188 0.00715375 0.0153836 0.0457208 -0.0213814 0.0183829 -0.0217617 -0.0241451 0.0562236 -0.0437782 0.00829286 0.0866145 -0.0442724 -0.0484038 0.0465044 0.0462053 -0.0131035 NA
HG00097 HG00097 -0.0316075 0.00138343 0.00602363 -0.0169859 0.00928749 -0.00330704 0.052626 -0.0140852 0.0236413 0.00480714 -0.0178737 0.0187442 -0.0834035 -0.000704524 -0.00238823 -0.0899047 0.00142339 0.0164213 0.0256184 -0.0573523 NA
```
## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/relate`](https://hub.docker.com/r/nfcore/relate/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/relate`](https://hub.docker.com/r/nfcore/relate/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
