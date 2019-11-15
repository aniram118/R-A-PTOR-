# RAPTOR: A tool for systematic identification of poly A tails and 3' unmapped regions from single molecule direct RNA sequencing

Raptor provides a comprehensive report of unmapped region (UMR)lengths, sequencences & composition, conserved poly A hexamer regions, nucleotide base composition,annotation information and a range of analysis plots.

The motivation behind this tool is to fill a gap in analysis tools as there is no current tool for analyzing 3' unmapped regions from long read sequencing data. Our hopes with this project are to continue to maintain this tool, implement new useful features, and respond to feedback from the community.

**Requirements**
  
- Python 2.7 or >= 3.0
- Matplotlib 
- seaborn 
- pandas 
- numpy 
- argparse 
- bedtools
- pybedtools

**Change Log**




**Installation**

git clone https://github.iu.edu/anigovin/RAPTOR.git


**Flags**


 Usage : raptor.py [-i < sam file >] [-b] [-hex] [-hum] [-u]
 - `-i`,`--input` Input Sam file
 - `-b`,`--bases` Generate a csv file of nucleotide composition of UMR sequences
 - `-hex`,`--hexamer` Generate a csv file and plots of Hexamer frequencies in the UMR sequences
 - `-u`,`--umrplots` Generate analsysis plots for UMR sequences
 - `-hum`,`--human` Generate a csv file containing comprehensive information of UMR sequences for Human sample
 - `-mo`,`--mouse` Generate a csv file of comprehensive information of UMR sequences for Mouse sample
 - `-f`,`--fly` Generate a csv file of comprehensive information of UMR sequences for Fly sample
 - `-y`,`--yeast` Generate a csv file of comprehensive information of UMR sequences for Saccharomyces Cerivisiae
 - `-h`,`--help` Show help message and exit


**Testing**

Once the repository has been downloaded, look for the file test_file.zip and unzip it which produces a test sam file.
Run as follows :
      
      python3 raptor.py -i test_file.sam -b -hex -u -hum

**Description of Output files** 

- 3_UMR_seqs.txt - A text file containing 3' UMR sequences 
- 3_UMR_tails.txt - A text file containing 3' UMR sequences and information for each read 
- 3_UMR.fastq - 3' UMR sequences in fastq file format 
- 3_UMR.bed - 3' UMR sequences in bed file format
- Base_composition.csv - A csv file contating individual nucleotide composition for UMRs in each read
- Hexamer_composition.csv - A csv file containing all possible hexamers and their percentage
- Top_Hexamers.pdf - A distribution plot of top 15 hexamers from the list of Hexamers detected
- UMR_analysis.csv - A csv file containing information on UMR sequences and associated features(read coordinates, length, etc)  and annotation information
- UMR_analysis.pdf - Analysis plots of UMR sequences 



**Credits**

@anigovin
@rkadumu

**License**
