This file contains data and scripts to estimate parameters
for models for the data from the paper:

    Rzhetsky, A., H. Shatkay, and W. J. Wilbur.  How to get the most
    out of your curation effort.  PLoS Computational
    Biology. 5(5). 2009.  

    doi: 10.1371/journal.pcbi.1000391

The data is a munged from of the data distributed with the paper as
"Dataset S1".  Specifically, the data was read into MATLAB and then
printed out in human-readable comma-separated-value (CSV) format.

After installing pyanno, you can run the data using our
sample MLE and MAP configurations, in either Windows or
Unix, using:

Maximum Likelihood Estimator (MLE):

     % cd $PYANNO_HOME
     % python examples/rzhetsky_2009/mle.py examples/rzhetsky_2009/data/p85.csv

Maximum a Posteriori Estimator (MAP):

     % cd $PYANNO_HOME
     % python examples/rzhetsky_2009/map.py examples/rzhetsky_2009/data/p85.csv

To change the priors for the MAP estimates, edit the map.py file.
Setting the priors all to 1.0 provides the same estimates as the MLE.

To run the other two annotations, replace "p85.csv" with "e85.csv" and
"f85.csv".  These have different shaped data, but the scripts are
generic enough to handle them all.

