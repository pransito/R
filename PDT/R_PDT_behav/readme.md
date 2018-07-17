R_PDT
-----

Code for publication of "Cue-induced effects on decision-making help detect subjects with gambling disorder"
by: Alexander Genauck (2018)


How to use
----------

Download the whole zip and put it into some working directory. The working directory you will have to indicate
at the beginning of two scripts:

1)
The machine learning part is started with the script "severity_pred_loop_v6.R"; put the working directory (root_wd)

2)
The univariate testing part / hierarchical regression (lme4) modeling part is done with the
/univariate_testing/glmer_accRate_la_cat.R script; also here set the working directory in the beginning

Careful: Both scripts are intense and take hours to run. You might want to just try out parts of it. The scripts
also install and load many R packages (see in the beginning of those scripts). They are all useful and should not
hurt your installation of R or bother your other code. However, revise the list before you run the code and decide
if you would like to continue.