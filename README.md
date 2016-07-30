# redBlocks

redBlocks: Filtering differential expression results for common co-expression pairs

========

### Summary

![summary](https://github.com/sarbal/redBlocks/blob/master/imgs/freq_tally.png "Method summary")

#################################################################################################
### 1. Setting up the environment
#################################################################################################

You first need to label the directory you have saved redBlocks in: eg. ``` masterdir="C:/redBlocks/" ```
And label your experiment: eg. ``` name = "my_experiment" ```
All output will be saved in: ``` out = paste(masterdir,name,sep="/") ```
You will also need to download the frequency tally network data from this link: ```  ```.

#################################################################################################
### 2. Expression data
#################################################################################################

#### a. To run redBlocks, you need a read counts dataset.
Currently, we have not implemented any pre-processing steps, so please make sure that the data is
set up as a matrix, with columns as your individual samples, and rows as genes.
The row names should be labelled by their gene entrez IDs.
The columns should be labelled by their sample IDs.
The expression dataset should be placed in the variable: ``` counts ``` 

#### b. Alternatively, you can use your DE output from either edgeR or DESeq2. 
Other methods can be used. 

#################################################################################################
### 3. Running redBlocks
#################################################################################################

Once the environment variables and the expression data is loaded, you can run the
script ``` run_redBlocks.R ``` :
In an R session: ``` source("run_redBlocks_on_data.R") ``` 


#################################################################################################
### 4. Output summary
#################################################################################################

![summary](https://github.com/sarbal/redBlocks/blob/master/imgs/output.png "Method summary")

#################################################################################################
### 5. Intepreting results
#################################################################################################


#################################################################################################
### 6. Extras
#################################################################################################
