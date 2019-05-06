# CircadianTools
A Collection of Tools for Detecting Rhythmic Genes
## Overview
Allows Cosinor Models and Turning Point Analysis to be easily carried out on transcriptomics data using R. 
## Install Guide
From R:
```{r}
install.packages("devtools")
devtools::("nathansam/CircadianTools")
```
## Cosinor Plotting
```{r}
cosinorplot("comp102333_c0_seq21", Laurasmappings)
```
![](cosinorex.png)


## Turnpoint Plotting
```{r}
turningplot("comp101252_c0_seq2", Laurasmappings)
```
![](turnpointex.png)


