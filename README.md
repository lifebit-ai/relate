# Ancestry and Relatedness inference workflow
Developed by Genomics England and Lifebit

## Introduction

The Ancestry and Relatedness pipeline is the second pipeline in a row of three Genomics England's pipelines for annotating aggregate genomic data:

[SiteQC](https://github.com/lifebit-ai/siteqc) > **Ancestry and Relatedness** > [Annotate](https://github.com/lifebit-ai/annotate)

The purpose of the current pipeline is to infer relatedness and ancestry of participants based on their germline variants and site quality metrics produced by SiteQC pipeline. The main output of Ancestry and Relatedness pipeline is participant ancestry assignment. This information together with all QC mitrics created by SiteQC pipeline is used in the final Annotate pipeline to add new annotations to raw input aggregate genotypic files.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

### Install dependencies: 
[`nextflow`](https://nf-co.re/usage/installation) and [`docker`](https://docs.docker.com/engine/installation/)

### Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run lifebit-ai/relate -profile test
```

### Start running your own analysis!
Modify the comand with your own input data.
```bash
nextflow run lifebit-ai/relate \
    --input                 's3://lifebit-featured-datasets/projects/gel/siteqc/input.csv' \
    --included_samples      's3://lifebit-featured-datasets/projects/gel/siteqc/sampleList.txt' \
    --siteqc_results_dir    's3://lifebit-featured-datasets/projects/gel/siteqc/ARtestFiles/siteqc_example_res' \
    --michigan_ld_file      's3://lifebit-featured-datasets/projects/gel/siteqc/MichiganLD_liftover_exclude_regions_PARSED.txt' \
    --ancestry_probs        's3://lifebit-featured-datasets/projects/gel/siteqc/aggV2_R9_M30K_1KGP3_ancestry_assignment_probs.tsv' \
    --pcs_ancestry          's3://lifebit-featured-datasets/projects/gel/siteqc/aggV2_bedmerge_30KSNPs_labkeyV9_08062020_update_PCsancestryrelated.tsv' \
    --unrelated_1kgp3       's3://lifebit-featured-datasets/projects/gel/siteqc/UNRELATED_1KGP3.samples' \
    --samplelist_1kgp3      's3://lifebit-featured-datasets/projects/gel/siteqc/1KGP3.sample_table' \
    --example_eigenvec      's3://lifebit-featured-datasets/projects/gel/siteqc/1KGP3_30K_unrel_autosomes.eigenvec' \
    --example_proj_eigenvec 's3://lifebit-featured-datasets/projects/gel/siteqc/1KGP3_30K_unrel_autosomes.eigenvec_TMPPROJ' \
    --super_pop_codes       's3://lifebit-featured-datasets/projects/gel/siteqc/super_pop_codes.tsv'

```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The lifebit-ai/relate pipeline comes with documentation about the pipeline which you can find in the [`docs/` directory](docs).

## Credits

The code for the `relate` workflow was originally written by Daniel Rhodes and the Genomics England Bioinformatics team.

The Nextflow implementation of the pipeline was developed by Christina Chatzipantsiou, Salvador Martinez and Vladyslav Dembrovskyi for use by Genomics England and collaborators.

## Pipeline DAG
Direct Acyclic Graph representation of the Ancestry and Relatedness pipeline.
- Bubbles describe individual steps of the pipeline, that in Nextflow are calledd `processes`.
- Arrows show the data flow between processes and correspond to Nextflow `channels`. Label on the right of each arrow (channel) is the name of the channel (or applied modification) that can help understand the kind of data transferred by the channel.

The presented graph shows the high level of pipeline structure. Each process and channel will be executed as many times as input data will require.

![image](https://user-images.githubusercontent.com/64809705/101158079-d2a20e80-3633-11eb-8374-ca4a463707dc.png)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

> <img src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-square.svg" width="20"/></a> **NOTE**: This pipeline was created using the nf-core template.  For further information or help with nf-core pipelines, you can get in touch with the core developers and community on [Slack](https://nfcore.slack.com/channels/lifebit-ai/relate) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  lifebit-ai/relate for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` pre-print as follows:  
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
