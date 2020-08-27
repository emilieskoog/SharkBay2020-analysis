# README

#### This will walk you through how to create a phylogenomic tree using anvi'o (with information overlaid such as MAG completeness, in this case!)

> Note: For more detailed explanation and further help (straight from the source), I highly recommend you pay the anvi'o workflow for phylogenomics tutorial a visit, which can be found [here](http://merenlab.org/2017/06/07/phylogenomics/). If you are interested in how I personally created this figure, then feel free to keep reading. 

## :memo: Phylogenomic tree reconstruction using anvi'o

### Step 1: Creating a conda environment for anvi'o
If you do not already have an anvi'o conda environment created, you can either follow instructions [here](http://merenlab.org/2016/06/26/installation-v2/) (straight from the source!) or use the conda environment that I have already created for you. You will need to follow the below directions for unpacking this environment. 

First, download the anvio5.5.yml file which can be found in this repository [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/PhylogenomicTree/anvio5.5.yml). Navigate to your terminal and make sure your path aligns with where your anvio5.5.yml file is. Execute the following:
```
conda env create -f anvio5.5.yml
```
You can verify that it was created by executing the following:
```
conda info --envs
```
Anvio5.5 should appear in this list. Sidenote - a fantastic resource for managing conda environments can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

> Note: This is version 5.5. Anvi'o has since come out with newer versions.  
> 
***Make sure to activate your newly created environment:***
```
conda activate anvio5.5
```




### Step 2: Generate an anvi’o contigs database for each MAG

You will need to create a directory with all MAGs of interest inside. 'cd' into this directory. The following code will help you with this (again taken straight from the Meren Lab webpage tutorial [here](http://merenlab.org/2017/06/07/phylogenomics/)):
```
cd your/path/to/MAGs

for i in *fa
do
	anvi-script-FASTA-to-contigs-db $i
done
```

MAG data for this project can be found [here](https://zenodo.org/record/3874996#.X0apURNKjRY). This step will result in the creation of new files with a '.db' file extension for every fasta file in the directory.

### Step 3: Get amino acid sequences out of contigs databases and concatenate them 

You will want to create a file that defines each contig database. You can see my example [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/PhylogenomicTree/combined_external_genomes.txt). Download this file into your working directory for use for the following command: 
```
anvi-get-sequences-for-hmm-hits --external-genomes combined_external_genomes.txt \
                                -o concatenated-proteins.fa \
                                --hmm-source Campbell_et_al \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate
```
I encourage you to read more about the ```--gene-names``` flag options on the Meren website tutorial and be intentional about the genes you choose for your own dataset. 

The outputted concatenated-proteins.fa file will be used in the final step.

### Step 4: Generate a newick-formatted tree

```
anvi-gen-phylogenomic-tree -f concatenated-proteins.fa \
                           -o phylogenomic-tree.txt
```

### Step 5: Generate a file with the information you wish to have included in the tree

View my file [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/PhylogenomicTree/info_overlay.txt) (and download to working directory for following step). My file includes information on different classifications of taxonomy as MAG completeness. I find Mike Lee's [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/six-glorious-commands) page quite useful for simple commands for creating, cutting, and pasting together tab-delimited files. 

### Step 6: Visualize your tree
FINALLY! The product of the hard work! Let's take a look.

> Note: If you are using an HPC for this step, you will need to set up an SSH tunnel with local and remote port forwarding which can be done by following [this](http://merenlab.org/2015/11/28/visualizing-from-a-server/) tutorial. If you are running it locally, feel free to run the following:

```
anvi-interactive -p phylogenomic-profile.db \
                 -d info_overlay.txt \
                 -t phylogenomic-tree.txt \
                 --title "Phylogenomic distribution and MAG completeness" \
                 --manual
```
If all went well, you should have something that (after some color additions and formatting in anvi-interactive) looks like this!


![](https://i.imgur.com/GrY3ALx.png)





