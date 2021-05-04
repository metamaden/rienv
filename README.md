# rienv

Author: Sean Maden

Programmatic setup for retained intron detection tools with `conda`. This resource currently 
supports environment setup for the following tools:

* [IntEREst](https://bioconductor.org/packages/release/bioc/html/IntEREst.html)
* [SIRFindeR](https://github.com/lbroseus/SIRFindeR/)
* [superintronic](https://github.com/sa-lee/superintronic)
* [iREAD](https://github.com/genemine/iread/)
* [Keep Me Around (KMA)](https://github.com/adamtongji/kma)

# Make environments

Run the shell script `conda.sh`:

```
sh conda.sh
```

Alternatively, set up environments from their `.yml` files, located at `./inst/yml/`:

```
conda env create -f ./inst/yml/env_interest.yml
```

# Tests

Scripts to test environments, especially using any provided vignette examples, are contained at `./inst/test`.

## Cross-tool tests

The test BAM files provided by the `IntEREst` package are contained in the folder `./inst/bam`. These consist of
a small example BAM `.bam` and BAM index `.bai` file aligned with hg19. The example BAM data was used to test 
each of the retained intron tools, where scripts and test results are contained at `./inst/bam/tests/tooltests/`.
