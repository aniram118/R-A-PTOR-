### Poly(A) tail predictor 

### ( Please contact admin to download code ) 
 Note : Current version supports sam files from human sample sequencing only 
## Prerequisites ##

Python3 or higher 

### Recommended ###

    Create Virtual environment using Anaconda3 or virtualenv

### Required Modules 

    Numpy
    itertools
    Matplotlib
    Seaborn
    Pandas

## Usage 

Extract file Raptor.tar.gz. This will generate a folder with the required scripts and files 

To run Raptor on test files 

    python raptor.py -i *file*.sam

## Output

Two csv files will be generated : 

1. ****tail_info.csv**** --->  Containing 3' terminal UMR region information 
2. ****hexamers.csv**** ---> Containing hexamer signal information for each read 




