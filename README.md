# tigaRcycle

The R-package `tigaRcycle` performs integrative detection of circadian signals in time series omics data, both sequencing counts RNAseq and continuous microarray data. The package offers a broad range of functions for visualisation of rithmic signals and comparison between different diets and time points. Also, pathway enrichment analysis and visualisation using different repositories can be performed comparing different conditions.

Note: if you have a choice to use either Windows or Unix/Linux, opt for the latter. `tigaR` runs more efficiently under Unix/Linux than under Windows. NOTE:  when running `tigaR` you may see *** WARNINGS ***  from `INLA` (e.g. on eigenvalues, or on convergence, or even something like 18500 Aborted...). They can currently not be suppressed, because they are produced by C-code. Please ignore them.

![image](https://user-images.githubusercontent.com/22052679/150277203-646d6d85-482a-44ab-8e30-c20d260179fe.png)

# References

Publication related to `tigaRcycle` include:

- Miok, V., Wilting, S.M., Van de Wiel, M.A., Jaspers, A., van Noort, P.I., Brakenhoff, R.H., Snijders, P.J.F., Steenbergen, R.D.M., Van Wieringen, W.N. (2014), "[tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities](https://doi.org/10.1186/1471-2105-15-327)", *BMC Bioinformatics*, 15, 327. 

Please cite the relevant publication if you use `tigaRcycle`.
