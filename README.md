**Genomics England's workflow for variant quality control and annotation**.

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run lifebit-ai/relate -profile test,<docker/singularity/conda>
```

iv. Start running your own analysis!

<!-- TODO nf-core: Update the default command above used to run the pipeline -->
```bash
nextflow run lifebit-ai/relate -profile <docker/singularity/conda> --input test.csv
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The lifebit-ai/relate pipeline comes with documentation about the pipeline which you can find in the [`docs/` directory](docs).

## Credits

The code for the `relate` workflow was originally written by Daniel Rhodes and the Genomics England Bioinformatics team.

The Nextflow implementation of the pipeline was developed by Christina Chatzipantsiou for use by Genomics England and collaborators.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

> <img src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo-square.svg" width="20"/></a> **NOTE**: This pipeline was created using the nf-core template.  For further information or help with nf-core pipelines, you can get in touch with the core developers and community on [Slack](https://nfcore.slack.com/channels/lifebit-ai/relate) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  lifebit-ai/relate for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` pre-print as follows:  
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
