# Phylogeographic Concordance Factors (PCFs)


___
**PCFs.py** will calculate phylogeographic concordance factors for a set of
co-distributed species. The method is described in Satler and Carstens
(2016).


___
## Instructions for use

Python script (written in version 2.7.6) uses mbsum and bucky (Ané et al. 2007;
Larget et al. 2010). You will want to download those two programs and place them in your
$PATH. Script was tested with bucky v1.4.3 but should work on the latest version.

Script uses species tree distributions as input; this has been tested with trees from
both *BEAST and SNAPP. Important things before running script:

1. OTUs need to be the **same** across all species, but order of the OTUs is not important.
2. Tree distributions should contain the same number of trees across each species.
If necessary, thin distributions until they are matching in length.

Can then place species tree distributions in a folder (separate from script), set the
tree folder as your current working directory, then call the script. This will run the
analysis, provide folders for the mbsum and bucky output, plus write a table summarizing
all bucky runs with their average PCF value. This can then be imported into excel and
sorted by *K* and PCF values.

___
example usage:

**python ../script/PCFs.py *.trees**
___

If you have any questions or comments, please email me at jordansatler@gmail.com.

___
References:

Ané, C., B. Larget, D.A. Baum, S.D. Smith, and A. Rokas. 2007. Bayesian estimation of
concordance among gene trees. *Molecular Biology and Evolution*, **24**:412–426.

Larget, B., S.K. Kotha, C.N. Dewey, and C. Ané. 2010. BUCKy: Gene tree / species tree
reconciliation with the Bayesian concordance analysis. *Bioinformatics*, **26**:2910–2911.

Satler, J.D., and B.C. Carstens. 2016. Phylogeographic concordance factors quantify
phylogeographic congruence among co-distributed species in the *Sarracenia alata*
pitcher plant system. *Evolution*, **70**:1105–1119..
