# madrat based MAgPIE Input Data Library

R package **mrMAgPIE**, version **0.6.0**

  

## Purpose and Functionality

Provides functions for MAgPIE country and cellular input data generation.


## Installation

For installation of the most recent package version an additional repository has to be added in R:

```r
options(repos = c(CRAN = "@CRAN@", pik = "https://rse.pik-potsdam.de/r/packages"))
```
The additional repository can be made available permanently by adding the line above to a file called `.Rprofile` stored in the home folder of your system (`Sys.glob("~")` in R returns the home directory).

After that the most recent version of the package can be installed using `install.packages`:

```r 
install.packages("mrMAgPIE")
```

Package updates can be installed using `update.packages` (make sure that the additional repository has been added before running that command):

```r 
update.packages()
```

## Questions / Problems

In case of questions / problems please contact Kristine Karstens <karstens@pik-potsdam.de>.

## Citation

To cite package **mrMAgPIE** in publications use:

Karstens K, Dietrich J, Chen D, Windisch M, Alves M, Beier F, Mishra A (2020). _mrMAgPIE: madrat
based MAgPIE Input Data Library_. R package version 0.6.0.

A BibTeX entry for LaTeX users is

 ```latex
@Manual{,
  title = {mrMAgPIE: madrat based MAgPIE Input Data Library},
  author = {Kristine Karstens and Jan Philipp Dietrich and David Chen and Michael Windisch and Marcos Alves and Felicitas Beier and Abhijeet Mishra},
  year = {2020},
  note = {R package version 0.6.0},
}
```

