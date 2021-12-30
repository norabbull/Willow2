



#############################################

Willow1.0 contain all program code necessary to reproduce results of the
master project: The GDR. A novel approach to detect large-scale
genomic sequence patterns, 2021.

Most code is written in Python version 3.8.8, and some additional code
is written in R version 4.0.3. 

############################################

Willow1.0 program is contained in the "src" folder, including: 

- treeInformation.py
    CONTAINS:   Core-class "treeInfo" to initiate and store information of a tree.

- treeMetrics.py
    CONTAINS:   Class "treeMetrics" to perform calculations of tree. Inherits
                "treeInfo".
- treeRun.py
    CONTAINS:   Class "treeRun" to configure tree instance calculations and 
                save reults for input files.
- treeMain.py   
    CONTAINS:   main function to read configuration and run.
    
See further descriptions in scripts.

############################################

Data wrangling and plots are found in the "treeInspection" folder, including: 


- treeHelpers.py
    CONTAINS: Help-functions to ease data-loading.
    
    
- dataWrangle.py
   CONTAINS: Code to format result data for plotting and inspection.
        - non-zero phydists
        - unique_seq
        - GDR values for super- and sub populations
        - GDR null-distribution values
"""
        
- plots.py
- single_gene_inspection.py

Relative paths used in project are included but must be changed in case of 
result reproduction. 


#############################################

