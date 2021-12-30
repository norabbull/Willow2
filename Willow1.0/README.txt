

################################################

Willow1.0 contain all program code necessary to reproduce results of the
master project: The GDR. A novel approach to detect large-scale
genomic sequence patterns, 2021.

Most code is written in Python version 3.8.8, and some additional code
is written in R version 4.0.3. 

################################################


### Overview:
  
- How to run Willow1.0
- Input file descriptions
  1. Configuration file
  2. Phylogenetic distance matrix files
  3. Group information file


### How to run: ################################

1. Terminal: 
$ python treeMain.py main_config.yml

2. In IDE: 
- Specify config file path in bottom of treeMain.py.

Job to run is specified in main_config file.

### Input files description #################### 


1. main_config.py
The configuration file (main_config.yml) in yaml-format contains all optional input parameters to Willow. 
It is the only input-argument given to run a Willow job. 
See `main_config_DESCRIPTION` and `main_config_EXAMPLE` (in documentation) for further details. 


2. Phylogenetic distance matrix files
Distance matrices are created from Newick-format phylogenetic distance files in R, see "R_scripts/newick2phydist.R" for details.
The distance matrices are stored in a folder, where all files in folder are processed in one job-run.
Path to folder is specified in main_config.
Files are in csv-format.

Example of a 3x3 distance matrix file: 

```
,AFR___LWK___NA19308,AFR___LWK___NA19317,EUR___GBR___HG00096
AFR___LWK___NA19308,0.0,0.00236,0.00236
AFR___LWK___NA1931,0.00236,0.00236,0.00118
EUR___GBR___HG00096,0.00236,0.00118,0.00354

```

All row- and column names refer to a sample/leaf node in the phylogenetic 
tree and must be formatted as a string that includes group information for 
all specified group categories. The group names, given by a 3-letter acronym, 
are separated by three underscores, "___". The order of specified groups 
in the string determines which group category the group belong to. 
This information corresponds to how group information is given in the group information file.    


3. Group information file
The group information file contains information of all specified groups for all categories. The file is at TSV-format, with tab-separated columns.
The file contains three columns: `GroupCategory`, `GroupName` and `GroupDescription`, see example below. 

EXAMPLE:
In this example, there are 2 different group categories specified: "SUB" and "SUPER". 
The group names "FIN", "GBR" and "IBS" are contained in the "SUB" category, whereas "EUR" and "EAS" are contained in the "SUPER" category.
The `GroupDescription` column contains optional additional group information (optional).
 

```
GroupCategory   GroupName   GroupDescription
SUB             FIN         Finnish in Finland
SUB             GBR         British in England and Scotland
SUB             IBS         Iberian Population in Spain
SUPER           EUR         European
SUPER           EAS         East Asians

```


########################################################################################
PROJECT DETAILS
########################################################################################

Configuration files to run calculations performed in the master-project is specified in "jobs".

Formatting and plotting of result data are contained in the "treeInspection" folder, including: 

- treeHelpers.py 		- Help-functions to ease data-loading  
    
- dataWrangle.py  		- Format result data for plotting and inspection.
        
- plots.py        		- Produce plots of result data

- single_gene_inspection.py	- 

Relative paths used in project are included but must be changed in case of 
result reproduction. 


#############################################

