# RAPTOR: A tool for systematic identification of poly A tails and 3' unmapped regions from single molecule direct RNA sequencing

Raptor provides a comprehensive report of unmapped region (UMR)lengths, sequence composition, conserved poly A hexamer regions, nucleotide base composition, and annotation information.

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


 Usage : raptor.py [-i < sam file >] [-b] [-hex] [-hum]
 - `-i`,`--input` Input Sam file
 - `-b`,`--bases` Generate a csv file of nucleotide composition of UMR sequences
 - `-hex`,`--hexamer` Generate a csv file and plots of Hexamer frequencies in the UMR sequences
 - `-hum`,`--human` Generate a csv file containing comprehensive information of UMR sequences for Human sample
 - `-mo`,`--mouse` Generate a csv file of comprehensive information of UMR sequences for Mouse sample
 - `-f`,`--fly` Generate a csv file of comprehensive information of UMR sequences for Fly sample
 - `-y`,`--yeast` Generate a csv file of comprehensive information of UMR sequences for Saccharomyces Cerivisiae
 - `-h`,`--help` Show help message and exit

Usage : transcript_plots.py [-i < fastq file >]
- `-i` Input fastq file of transcripts

Usage : umr_plots.py [-i < fastq file >]
- `-i` Input fastq file of UMRs from raptor.py script

**Testing**

Once the repository has been downloaded, look for the file test_file.zip and unzip it which produces a test sam file.
Run as follows :
      
      python3 raptor.py -i test_file.sam -b -hex -hum
      
      python3 transcript_plots.py transcript.fastq
      
      python3 umr_plots.py 3_UMR.fastq



**Credits**

@anigovin
@rkadumu

**License**
