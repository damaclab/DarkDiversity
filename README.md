# Dark Diversity

This repository contains code and instructions for computing dark diversity before performing association rule mining on a species occurrence dataset to extract meaningful information from the absent data of species occurrence using R.

## Table of Contents

- [Introduction](#Introduction)
An enormous load of growing biodiversity data needs algorithmic care for accurate data management, and therefore the term computational biodiversity comes. Instead of relying purely on presence data, the probabilistic forecast of member distribution including the regions of not occurrence can neutralize biodiversity loss by restoring potential ecosystems. This work is aiming at revealing the perspective of computational biodiversity as a counteract for biodiversity loss by correlating the concept of dark diversity. The computation of the dark
diversity is accompanied by a data mining algorithm for establishing rules with more nobility to manage the depletion of biodiversity.
We generate a dataset for the Indian estuarine ecosystem and show the use of our approach by ending up with rules worthwhile for the ecologists.
These would step up reinforcing biological diversity via introducing or rehabilitating specific faunal groups to an estuary under survey.


- [Prerequisites](#Prerequisites)
  
- Before using this code, ensure that the following prerequisites are needed to be installed:

•	R: It can be downloaded and installed from the official website: https://www.r-project.org/

•	RStudio: RStudio is an integrated development environment (IDE) for R. It provides a user-friendly interface and additional tools for R programming. RStudio can be downloaded and installed from the official website: https://www.rstudio.com/.




- [Installation](#Installation)

- To use this code, follow these steps:
1.	Clone this repository to a local machine or download the ZIP file and extract it.
2.	Open the R environment or RStudio.
3.	Set the working directory to the location where it is cloned or extracted from   the repository.
4.	Install the required R packages by running the following command in the R console: For e.g.

|install.packages(c("arules", "arulesViz"))|
| :- |
- [Data Preparation and computing dark diversity](#DataPreration)
    - We use 75\% of the dataset as training data. Then we predict the probable occurrences of not-found classes at multiple estuaries for the remaining 25\% of 
      the dataset. UNO method uses any of the 3 different methods (minobs(), minpred(), binminpred()) to calculate threshold value for predicting abundances. 
      After having tested with all the 3 different methods, it has been found that minobs() and minpreds() are suitable for our case. minpred() uses the 
      smallest predicted value for cases where a class observed and minobs() uses the smallest positive value observed for a class. Data is first binarized in 
      the case of binminpred().Results obtained from the UNO function will be preprocessed.

  - Min-max normalization is performed on the result of UNO Function.
    
  - Discretization of all the obtained values and apply textual leveling to make the dataset comprehensible.
 
  - Generating a set of association rules. 

## ORIGINAL PAPER
* **Moumita Ghosh, Anirban Roy, Kartick Chandra Mondal** - *Determining Dark Diversity of Different Faunal Groups in Indian Estuarine Ecosystem: A New Approach With Computational Biodiversity* - [Link To Paper](https://link.springer.com/chapter/10.1007/978-981-16-4435-1_16)
