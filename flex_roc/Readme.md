# Code for flexible ROC area calculation.
*flex_roc.cc*

	USAGE: ./fex_roc [datafile] [width] [opt: roc output file]


* A positive prediction +/- width registers as success.
* Note: setting width=0 is standard ROC calculation


	* dependency: gsl (for numerical integration)
	* compile command: g++ -O3 flex_roc.cc -o flex_roc -lgsl


**Copyright 2011 Ishanu Chattopadhyay <ishanu@uchicago.edu>**

