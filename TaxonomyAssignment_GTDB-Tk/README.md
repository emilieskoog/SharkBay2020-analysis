# README

#### This file will walk you through taxonomy assignment for all Shark Bay MAGs using GTDB-Tk.

>Note: This analysis was done using v0.2.2 of GTDB-Tk. More recent versions have come out since.


## :memo: Assigning Taxonomy with GTDB-Tk

### Step 1: Download GTDB-Tk

The first step, of course, is to acquire the program. Because the program is large, I recommend installing GTDB-Tk on an HPC. You can find instructions on how to download the program on the GTDB-Tk website [here](https://github.com/Ecogenomics/GTDBTk). I find that the Bioconda installation and creating a new conda environment best. 

```
conda create --name gtdbtk

conda activate gtdbtk

conda install -c bioconda gtdbtk
```
### Step 2: Run it!

The MAG files used can be found [here](https://zenodo.org/record/3874996#.X0gJmBNKhLY). You want to make sure you download this folder into your working directory. Once downloaded, execute this command in your working directory. 

```
gtdbtk classify_wf --cpus 20 --genome_dir MAGS --out_dir outfile_MAGS_gtdbtk --extension fa
```

GTDB-Tk just produced many many files for you. You can view a few of the output files [here](https://github.com/emilieskoog/SharkBay2020-analysis/tree/master/TaxonomyAssignment_GTDB-Tk/outfile_MAGS_gtdbtk). The most essential file for our purposes is that which breaks down the taxonomy for each MAG. You can take a look at this summary file [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/TaxonomyAssignment_GTDB-Tk/outfile_MAGS_gtdbtk/gtdbtk.bac120.classification_pplacer2.txt). 

Now that we have assigned taxonomy using GTDB-Tk, we can make more sense of these MAGs and place the "who" in context of the "what".  