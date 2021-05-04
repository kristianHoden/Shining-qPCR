# Shining-qPCR
Shiny app that can analyze the qPCR data and contains multiple barplot plotting options.

Run it in R with the following command:

runGitHub( "Shining-qPCR", "KristianHoden")

<a href="https://zenodo.org/badge/latestdoi/351428022"><img src="https://zenodo.org/badge/351428022.svg" alt="DOI"></a>

## Features
ddCt qPCR analysis  
Outlier analysis and removal  
Label renaming  
Ordering of samples  
Automatic statistics (Both ANOVA and different t-tests)  
Easy plotting and exporting  

## Preparations
### Before running the qPCR machine
Please label the wells Sample-Target-ReplicateNumber. E.g. Samp1-Targ1-1 or Condition1-Gene1-1. See the first column in the example excel document (exampleSheet) for further details.

### After running the qPCR machine but before running the Shining-qPCR app
Collect the wells column, Cq column and (if desired) Melt Temperature column in a separate sheet. See the example sheet (exampleSheet) for further details. The __Shining-qPCR__ app automatically sets the last sheet of the excel document. No further manual processing is needed.   

## When running the Shining-qPCR app

The app is quite intuitive. Start with choosing what excel document to analyze and the different analysis options will appear. Please remember to push the __ddCt__ button before going to the Results tab.  

Don't hesitate if you have any questions.
