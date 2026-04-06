#   IICD Parameter Estimation Workshop: Parameter Inference Approaches for Data-Driven Mathematical Modeling

##  Installation

This library can be installed with

```{r}
devtools::install_github("dinhngockhanh/IICDWorkshopParameterInference")
```

If there are issues with library `abcsmcrf`, it can be installed with

```{r}
devtools::install_github("dinhngockhanh/abcsmcrf")
```

Please check that the following libraries are installed and can be used before the Workshop:

```{r}
library(deSolve)
library(truncnorm)
library(crayon)
library(ggplot2)
library(EasyABC)
library(parallel)
library(abcsmcrf)
library(IICDWorkshopParameterInference)
```

##  References
1.  Dinh KN, Liu C, Xiang Z, Liu Z, Tavaré S. Approximate Bayesian computation sequential Monte Carlo via random forests. Stat Comput 35, 219 (2025). https://doi.org/10.1007/s11222-025-10748-x
2.  Sisson SA, Fan Y, Beaumont M, editors. Handbook of approximate Bayesian computation. CRC press; 2018 Sep 3. https://doi.org/10.1201/9781315117195
3.  Tavaré S, Balding DJ, Griffiths RC, Donnelly P. Inferring coalescence times from DNA sequence data. Genetics. 1997 Feb 1;145(2):505-18. https://doi.org/10.1093/genetics/145.2.505
4.  Beaumont MA, Cornuet JM, Marin JM, Robert CP. Adaptive approximate Bayesian computation. Biometrika. 2009 Dec 1;96(4):983-90. https://doi.org/10.1093/biomet/asp052
