Under the supervision of Dr. Osvaldo Espin-Garcia, I ported an existing genetic algorithm routine from native R to C using Rcpp and RcppArmadillo in order to decrease runtime. 

The files _kofnGA.cpp_ and _kofnGArcpp.R_ were written by myself (with some code utilized from external sources). _optimJTC_v2_rcpp.R_ was written by Dr. Osvaldo Espin-Garcia prior to our collaboration - I only added 2 blocks of code (810 - 820, 826 - 831) to it in order for the file to comply with the code written in _kofnGA.cpp_ and _kofnGArcpp.R_. 
