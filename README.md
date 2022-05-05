# expressyouRcell-examples

Go to https://github.com/gittina/expressyouRcell and follow the instructions for installing the package.

Then, download the example.R script and the ```ensg_dt.RData``` data.table. The latter will be used in the example.R script to associate gene IDs with gene symbols. 

Now, you can set you working directory within the ```setwd("hereyouhavetowriteyourfolderpath")``` function.
Execute now the example.R script. Two folders will be created in the first lines of code. 
```
dir.create("data")
dir.create("movies")
```

The "data" and "movies" folders, for storing the datasets reported in our paper and for storing the animated pictograms produced by expressyouRcell, respectively. 

Download the datasets as reported in the Key Resource Table (see STAR Methods) and save them within the "data" folder previously created.

To execute the map_gene_localization() function you first need to download a GTF file. Please, go to https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.primary_assembly.annotation.gtf.gz and store it within your working directory.

As reported in the example.R script you can now execute:
```
gene_loc_table <- map_gene_localization(gene_set = "gencode.vM22.primary_assembly.annotation.gtf")
```

Otherwise you can just download the gene_loc_table.RData and load it to your R environment.

```
load("gene_loc_table.RData")
```
