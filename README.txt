README: pyanno 1.0
============================================================
pyanno 1.0 is a suite of Python libraries for modeling the data coding
or diagnostic testing process to support inferences required for data
curation.  Specifically, pyanno implements statistical models for
inferring from multiply coded categorical data 

    * annotator accuracies and biases,
    * gold standard categories of items,
    * prevalence of categories in population, and
    * population distribution of annotator accuracies and biases.

The models include Dawid and Skene's (1979) multinomial model with a
maximum likelihood estimator implemented with EM.  A generalization of
Dawid and Skene's model with Dirichlet priors on prevalence and
estimator accuracy is also implemented.  This model is also estimated
with EM.  Finally, Rzhetsky et al.'s (2009) version of the Dawid and
Skene model with tied accuracies and uniform errors.


LICENSING
------------------------------------------------------------
pyAnno is licensed under the Apache License, Version 2.0.
For more information, see

     License.txt 



INSTALLATION
------------------------------------------------------------
pyAnno comes with a standard Python installer, which
is described in:

     Installation.txt



USE
------------------------------------------------------------
pyAnno's user guide is at:

     UserGuide.txt



UNDERSTANDING
------------------------------------------------------------
The statistical models on which pyAnno is based are
described in:

     Models.txt



DEVELOPING
------------------------------------------------------------
Links to the version control repository and explanations
of coding conventions and unit testing are in:

     DeveloperGuide.txt



CONTRIBUTORS
------------------------------------------------------------
Bob Carpenter    (Columbia University,   Statistics)
Andrey Rzhetsky  (University of Chicago, Medicine)
James Evans      (University of Chicago, Sociology)


