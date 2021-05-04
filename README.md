# rienv

Programmatic setup for retained intron detection tools with `conda`. This resource currently 
supports environment setup for the following tools:

* [IntEREst](https://bioconductor.org/packages/release/bioc/html/IntEREst.html), 
* [SIRFindeR](https://github.com/lbroseus/SIRFindeR/)
* [superintronic](https://github.com/sa-lee/superintronic)
* [iREAD](https://github.com/genemine/iread/)
* [Keep Me Around](https://github.com/adamtongji/kma).

# Make environments

Run the shell script `conda.sh`:

```
sh conda.sh
```

Alternatively, set up environments from their `.yml` files, located at `./inst/yml/`:

```
conda env create -f ./inst/yml/env_interest.yml
```
