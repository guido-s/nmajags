##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s)
##
setwd("..")
document("nmajags")


##
## (3) Build R package and PDF file with help pages
##
build("nmajags")
build_manual("nmajags")


##
## (4) Install R package
##
install("nmajags")


##
## (5) Check R package
##
check("nmajags")


##
## (6) Check examples
##
run_examples("nmajags", fresh = TRUE,
             run_dontrun = TRUE, run_donttest = TRUE)
warnings()
