I checked it in rather confusingly into two R files:
 
   pyanno/examples/sim.R

   pyanno/examples/fit.R

The input format is exactly the same for the fit
program, and it's what produced by the simulator.  

First, install the MCMCpack package (either with this
on the command line or however you do it in R studio):

> install.packages("MCMCpack");

Then set the working directory to the examples:

> setwd('/Users/carp/github/pyanno/examples')

Then run the simulator:

> source('sim.R')
Loading required package: coda
Loading required package: lattice
Loading required package: MASS
##
## Markov Chain Monte Carlo Package (MCMCpack)
## Copyright (C) 2003-2013 Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park
##
## Support provided by the U.S. National Science Foundation
## (Grants SES-0350646 and SES-0350613)
##

You may need to install.packages("coda") if you don't get it
for free with MCMCpack.

Then you can look at the data it created, which starts like this:

1	12	1
1	13	3
2	1	1
2	5	3
...


It is tab separated data with rows consisting of:

   item-id    annotator-id    label-id

All numbering is now from 1, I'm afraid.   This'll make
life easier downstream when we use Stan, too.

Then, you can run EM to fit the Dawid-Skene model with:

> source('fit.R')
[1] "K= 3 ; N= 5944 ; J= 20 ;I= 1000"
[1] "epoch= 1  log likelihood= -3890.67263325177"
[1] "epoch= 2  log likelihood= -3876.80855989868  relative_diff= 0.00356341297764285"
[1] "epoch= 3  log likelihood= -3875.9051763058  relative_diff= 0.00023302249232276"
[1] "epoch= 4  log likelihood= -3875.64751549489  relative_diff= 6.64775837333088e-05"
[1] "epoch= 5  log likelihood= -3875.52832556922  relative_diff= 3.07535515537382e-05"
[1] "epoch= 6  log likelihood= -3875.46482599982  relative_diff= 1.63847517205452e-05"
[1] "epoch= 7  log likelihood= -3875.42897802664  relative_diff= 9.24998027095082e-06"
[1] "epoch= 8  log likelihood= -3875.40810194752  relative_diff= 5.38677891939462e-06"
[1] "epoch= 9  log likelihood= -3875.39570203434  relative_diff= 3.19964061892114e-06"
[1] "epoch= 10  log likelihood= -3875.38823455795  relative_diff= 1.9268939125221e-06"
[1] "epoch= 11  log likelihood= -3875.38369183004  relative_diff= 1.17219943905468e-06"
[1] "epoch= 12  log likelihood= -3875.380907116  relative_diff= 7.18564731456025e-07"
[1] "epoch= 13  log likelihood= -3875.37918992909  relative_diff= 4.4310145321872e-07"
[1] "epoch= 14  log likelihood= -3875.37812606752  relative_diff= 2.74518056893307e-07"
[1] "epoch= 15  log likelihood= -3875.37746449237  relative_diff= 1.70712413497694e-07"
[1] "epoch= 16  log likelihood= -3875.37705183265  relative_diff= 1.06482458266986e-07"
[1] "epoch= 17  log likelihood= -3875.37679379397  relative_diff= 6.65841474909895e-08"
[1] "epoch= 18  log likelihood= -3875.37663210938  relative_diff= 4.17209996550323e-08"
[1] "epoch= 19  log likelihood= -3875.37653062658  relative_diff= 2.6186563028152e-08"
[1] "epoch= 20  log likelihood= -3875.3764668394  relative_diff= 1.64596090525843e-08"
[1] "epoch= 21  log likelihood= -3875.37642669812  relative_diff= 1.03580350701242e-08"
[1] "epoch= 22  log likelihood= -3875.37640141196  relative_diff= 6.52482670431205e-09"
[1] "FINISHED."

You'll find that this writes three files, 

  pi_hat.tsv     : prevalence estimate

  theta_hat.tsv  : annotator response estimates

  y_hat.tsv      : true category estimates

They all have headers which should make it clear.  I put them all into
DB style.

pi_hat.tsv
------------
category	prob
1	0.135940187127211
2	0.421016790400775
3	0.443043022472014

theta_hat.tsv
-------------
annotator	reference	response	prob
1	1	1	0.642676289401797
1	1	2	0.357323710598203
1	1	3	3.38043285559312e-48
1	2	1	0.117765826491172
1	2	2	0.83435279149599
1	2	3	0.0478813820128381
1	3	1	0.141760703967389
1	3	2	0.0592905873841705
1	3	3	0.798948708648441
2	1	1	0.924339004923274
2	1	2	3.91032986384067e-42

z_hat.tsv
---------
item	category	prob
1	1	0.497784635325564
1	2	0.0272926868690016
1	3	0.474922677805434
2	1	5.93583169362603e-34
2	2	0.00802662850296925
2	3	0.991973371497031
3	1	3.28892862843673e-86

