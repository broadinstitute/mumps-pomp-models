## Code for SEIR  mumps models

This folder contains R notebooks and scripts to calculate the MLE of the parameters of the mumps epidemic models for Harvard University and Ohio State University.

### How to use

The .Rmd notebooks can be run from within RStudio, using the desired properties files. Otherwise, the .R scripts for Harvard and OSU can be run from the command line as follows:

```
./run_harvard fast.properties
```

```
./run_ohio fast.properties
```

### Publication with resutls

The results from these scripts are presented in the following pre-print:

https://www.medrxiv.org/content/10.1101/2020.07.31.20166348v2

~Note:~ The calculations use the "opt" parameter set, so the reported point estimates of the epidemiological parameters correspond to the maximum of the log-likelihood function.