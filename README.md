# Information Browser for CFTR Regulatory Regions

In the process of researching the regulation of CFTR, information from literature review becomes difficult to sift through.
This project includes a very light "database" for organizing basic info as well as an RShiny app for browsing. The goal is 
to be able to use this to efficiently ask questions about enhancer specificity, chromatin interactions, and transcription
factor targeting while keeping track of sources. 

## Basic Info

If downloading for the first time the setupDB.R file should be run to make the sqlite file with necessary tables. After that,
there should be no issue running the Shiny app. Note that while data entry is allowed, data removal is not implemented through
the Shiny app just yet. That is a TODO item that will be noted later.

### Package Installation

Only R packages are needed to run this project. The renv lockfile should have everything needed.

## Data Structuring Rationale

The point of keeping track of experiment type is to gauge how reliable/robust the result is. Negative result are also desired
here because they can tell us a bit about cell type specificity. 

## Future Plans

Add ability to remove data, particularly if there is a need to remove all data associated with specific publications. 

## Authors

Carter Nakagawa
