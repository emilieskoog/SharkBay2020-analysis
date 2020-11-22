# README

#### Below, I walk through how I blasted for genes within my MAGs in order to create a presence/absence "heatmap" in R. 

## :memo: Acquiring gene hits from dataset

#### Before we can plot any results, we need to first identify hits for each gene of interest within each MAG. I did this using NCBI BLAST v. 2.6.0 on an HPC. For four of my six genes of interest - sulfotransferase, sulfatase, formylglycine-generating enzyme, and polysaccharide exporter gene - I acquired the PFAM identifiers from the NCBI conserved domain database. However, the glycosyltransferase and glycosyl hydrolase genes each have MANY different families of genes each with a different specificity for transferring sugars to different substrates, so JGI IMG annotations were first used to identify the presence of each unique glycosyltransferase and glycosyl hydrolase PFAM hit present within each MAG. Each of these PFAM hits were then taken from the NCBI conserved domain database and blasted using BLAST 2.6.0+ to identify e-value and percent identity of each glycosyltransferase and glycosyl hydrolase gene (which was not given by JGI IMG). I will break down the following workflow into the two ways these different gene hit files were created, first focusing on the protein sequences taken from NCBI directly and then those first identified in JGI IMG.  

## :one: NCBI 


### Step 1: Place all MAGs into one directory

Each of the MAGs used in this dataset can be found [here](dx.doi.org/10.5281/zenodo.3874996). Drag and drop or use the command line (whatever makes you happiest!) to move all MAGs into a newly created folder. Navigate to this folder in your terminal.

### Step 2: Make a file with all names of MAGs in them
Now that we are in this folder, we want to create a text file with the names of each of our MAGs. We can do this with a simple command:
```
ls > mags.txt
```
So now we have created a new file with a list of everything in this folder. However, that also includes the name of this newly created text file, in this case mags.txt. We do not want to have this text file name within this file, so we can simply "nano" into the text file 

```
nano mags.txt
```
and remove that one line (mags.txt). Make sure to save the modified file.

### Step 3: Create a database of all MAGs

Next, we want to create a BLASTable database for each of our MAGs. This following example is an array job where we are submitting many similar jobs using a single script. In my case, I have 84 genomes, so I am creating 84 jobs in a single execution. 

:exclamation: An important note, we are creating protein databases out of our NUCLEOTIDE MAG files to blast for PROTEIN queries. :exclamation: 
```
#!/bin/bash
#SBATCH -p your/partition/here      
#SBATCH -N 1                          
#SBATCH -n 1            
#SBATCH --mem=10G                     
#SBATCH --time=1:00:00                
#SBATCH -J build_db            
#SBATCH --output=dbbuild_%a.out  
#SBATCH --error=dbbuild__%a.err
#SBATCH --array=1-84 #change this last number to reflect your number of MAGs or genomes


FILENAME=$(awk "NR==$SLURM_ARRAY_TASK_ID" mags.txt) #change to reflect your job scheduling system

conda activate blast #load BLAST here

makeblastdb -in $FILENAME -parse_seqids -dbtype nucl -out $FILENAME  
```
You can visit [this site](https://www.unix.com/shell-programming-and-scripting/118605-awk-whats-nf-nr.html) to understand more about awk and NR.

### Step 4: Create a fasta file for each desired protein you are BLASTing

For example, I acquired the sulfotransferase PROTEIN sequence that I wanted to search for within each of my MAGs. 

```
>sulfo_PF13469
PIFIVGLPRSGTTLLHRLLAAHPQVRPPEETVIPILALLQSGRELRRLLDALTREDAELPHGPEECWQLLRQAFASFILEALARVSYARWLCDKSPSHLFYLDLLLRLFPDAKFIHLVRPDVISSYCSLSYSDFLDQIlARWARAYMAARARLPPDRFLDVRYEDLVADPEGTLRRIYDFLGLPWDD
```
I created a new file

```
nano sulfo.fasta
```
and pasted the above sequence into the file.

My file was in the same directory as my MAGs and mags.txt file.

### Step 5: BLAST proteins agains (nucleotide) MAGs

Once again, I created a job array while blasting (just makes things faster).
```
#!/bin/bash
#SBATCH -p your/partition/here      
#SBATCH -N 1                          
#SBATCH -n 1            
#SBATCH --mem=10G                     
#SBATCH --time=1:00:00                
#SBATCH -J blasts              
#SBATCH --output=blasts_%a.out     
#SBATCH --error=blasts__%a.err      
#SBATCH --array=1-84

conda activate blast #load BLAST here

FILENAME=$(awk "NR==$SLURM_ARRAY_TASK_ID" mags.txt) 

tblastn -query sulfo.fasta -db $FILENAME -out $FILENAME+"hits_sulfotransferase" -outfmt "7 std"
```

BLAST is now going through each of our MAG databases that we created in Step 3 and searching for our protein query, in this case sulfotransferase. And each sulfotransferase hit for each MAG will be outputted into a single file that will include the filename for each MAG and append the "+hits_sulfotransferase" to the file name. Here we chose the standard output format which includes percent identity, e-value, and many other details for each hit found. The output files will have the following structure: 


```
# TBLASTN 2.6.0+
# Query: sulfo_PF13469
# Database: SB_biofilm_MAG_10.fa
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 13 hits found
sulfo_PF13469   k141_434019     32.381  210     117     6       2       186     62714   63343   1.64e-14        68.6
sulfo_PF13469   k141_117837     28.497  193     123     5       1       183     87658   87095   8.31e-09        52.0
sulfo_PF13469   k141_117837     31.111  45      31      0       50      94      67863   67729   1.3     27.3
sulfo_PF13469   k141_221753     46.667  30      16      0       156     185     22597   22508   0.011   33.5
sulfo_PF13469   k141_51103      45.161  31      17      0       1       31      12128   12220   0.24    29.6
sulfo_PF13469   k141_51103      54.545  22      10      0       4       25      6887    6952    2.6     26.2
sulfo_PF13469   k141_51103      27.586  29      21      0       157     185     12686   12772   3.7     25.8
sulfo_PF13469   k141_51103      45.455  22      12      0       4       25      5825    5890    3.8     25.8
sulfo_PF13469   k141_81607      39.394  33      20      0       153     185     52643   52545   0.27    29.3
sulfo_PF13469   k141_386313     57.895  19      8       0       2       20      95397   95453   2.8     26.2
sulfo_PF13469   k141_546304     33.333  24      16      0       2       25      32506   32435   3.1     26.2
sulfo_PF13469   k141_291400     56.250  16      7       0       2       17      186     139     4.5     25.4
sulfo_PF13469   k141_282143     44.000  25      14      0       2       26      22593   22667   8.9     24.6
# BLAST processed 1 queries
```
## :memo: Preparing acquired data for visualization
Now that we have all of our hits for each gene for each MAG, it is time to put it all in a table that we can input into R. This first requires some file manipulation/ bash magic.

### Step 1: Move all output files into a single directory

I blasted several different proteins, so I have ended up with many outputted hit files. We first want to move them to one directory. 

```
mkdir HITS #made a directory called HITS
mv *+hits_* HITS #moved all outputted files with the `+hits_` "motif" into the HITS directory 
```
### Step 2: Move HITS directory to local computer
This is the point at which I openned a new terminal window and moved my HITS directory to my local computer. 
```
scp -r path/to/HITS/directory/on/HPC path/to/local/computer/directory/where/you/want/HITS/folder
```

### Step 3: Now for file manipulation/restructuring

I first navigated to my local HITS directory and removed the first 5 lines of each file which were consistent in all output files and looked like the following:
```
# TBLASTN 2.6.0+
# Query: sulfo_PF13469
# Database: SB_biofilm_MAG_10.fa
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 13 hits found
```
Remove first 5 lines:
```
for f in *; do
    tail -n +6 "$f" > "${f}".tmp && mv "${f}".tmp "$f"
    echo "Processing $f"
done
```
I then removed the last line of each file which was also consistent:
```
# BLAST processed 1 queries
```
Remove last line:


```
for f in *; do
    sed '$d' "$f" > "${f}".tmp && mv "${f}".tmp "$f"
    echo "Processing $f"
done
```
Now we have all of our output files that only have out hits like so:
```
sulfo_PF13469   k141_434019     32.381  210     117     6       2       186     62714   63343   1.64e-14        68.6
sulfo_PF13469   k141_117837     28.497  193     123     5       1       183     87658   87095   8.31e-09        52.0
sulfo_PF13469   k141_117837     31.111  45      31      0       50      94      67863   67729   1.3     27.3
sulfo_PF13469   k141_221753     46.667  30      16      0       156     185     22597   22508   0.011   33.5
sulfo_PF13469   k141_51103      45.161  31      17      0       1       31      12128   12220   0.24    29.6
sulfo_PF13469   k141_51103      54.545  22      10      0       4       25      6887    6952    2.6     26.2
sulfo_PF13469   k141_51103      27.586  29      21      0       157     185     12686   12772   3.7     25.8
sulfo_PF13469   k141_51103      45.455  22      12      0       4       25      5825    5890    3.8     25.8
sulfo_PF13469   k141_81607      39.394  33      20      0       153     185     52643   52545   0.27    29.3
sulfo_PF13469   k141_386313     57.895  19      8       0       2       20      95397   95453   2.8     26.2
sulfo_PF13469   k141_546304     33.333  24      16      0       2       25      32506   32435   3.1     26.2
sulfo_PF13469   k141_291400     56.250  16      7       0       2       17      186     139     4.5     25.4
sulfo_PF13469   k141_282143     44.000  25      14      0       2       26      22593   22667   8.9     24.6

```
### Step 4: Combining all files

#### Here we merge all files together and include the file name to each line so that we can track which MAG each hit was from. 

```
awk -F, 'NR==1{print "SB_biofilm_MAG_1.fa+hits_FGE",$0;next}{print FILENAME , $0}' OFS=, * > combined_hits.csv
```

Note that I have 'SB_biofilm_MAG_1.fa+hits_FGE' because that is the name of my first file. You could also do
```
awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' OFS=, * > combined_hits.csv
```
but then once you take a look at the file, the first file name is "file_name" like the below image demonstrates:

![](https://i.imgur.com/218AW2x.png)

When really what we want is:

![](https://i.imgur.com/05btC1N.png)

This is why I used "SB_biofilm_MAG_1.fa+hits_FGE" directly, because I know that that is the first file to be read, so it saves having to go through and edit it myself. :heavy_exclamation_mark: But for files that you are using that are named differently, you would need to make sure to use the proper name accordingly. 

The combined.csv file was then openned in excel and the headers "MAG" and "info" were added like so:

![](https://i.imgur.com/DbYwT3i.png)

This file was saved as a csv file and imported into R.

## :one: JGI IMG annotatoins :arrow_right: NCBI

### Step 1: Acquired JGI IMG annotations 

First, JGI IMG was used to annotate all MAGs. A bulk download of PFAM annotations for each MAG was done. A file like such was created for each MAG:
![](https://i.imgur.com/nU7NMxn.png)

### Step 2: Modify files to prepare for merging

First, navigate to your terminal and 'cd' into your directory that contains all of these PFAM hit files (see above image of example file). Before merging all files, the first line of each file (the header) is removed:
```
for f in *; do
    tail -n +2 "$f" > "${f}".tmp && mv "${f}".tmp "$f"
    echo "Processing $f"
done
```
### Step 3: Merge all files 

This can be done similar to the previous workflow:
```
awk -F, 'NR==1{print "MAG_1.txt",$0;next}{print FILENAME , $0}' OFS=, * > combined_pfams.txt
```

### Step 4: Prepare a different file for each of the glycosyltransferase and glycosyl hydrolase gene hits 

For the glycosyltransferase gene hits, we create a new file named combined_pfams_GT.txt which only includes glycosyltransferase PFAM hits:

```
sed '/Glyco.*tr/!d' combined_pfams.txt > combined_pfams_GT.txt
```
For the glycosyl hydrolase gene hits, we create a new file named combined_pfams_GH.txt which only includes glycosyl hydrolase PFAM hits:
```
sed '/Glyco.*hyd/!d' combined_pfams.txt > combined_pfams_GH.txt
```
### Step 5: Modify files to prepare for loading into R

I converted each file into csv format.

convert to csv and add 2 column headers : MAG and info in excel
```
sed 's/ \+/,/g' combined_pfams_GT.txt > combined_pfams_GT.csv

sed 's/ \+/,/g' combined_pfams_GH.txt > combined_pfams_GH.csv
```
### Step 6: Identify unique pfam hits in R

We then loaded these files into R and idenitified the unique PFAM hits (using the 'unique' funtion in R) for each of the glycosyltransferase and glycosyl hydrolase genes (since there are many subfamilies). 

```
GT_hits <- read.csv("~/Desktop/combined_pfams_GT.csv") 
GT_df <- GT_hits %>% 
  separate(info, sep="\t", into = c("gene_oid",	"gene_length", "query_start",	"query_end",	"subj_start",	"subj_end", "evalue",	"bit_score",	"pfam_id",	"pfam_name",	"pfam_length"), remove=TRUE, convert = FALSE)%>%
  data.frame
GT_unique_pfams <- unique(GT_df$pfam_id)
View(GT_unique_pfams)
```
>Note: The column headers are the headers that were removed earlier in Step 2. 

The same was executed for the glycosyl hydrolase hits. 

```
GH_hits <- read.csv("~/Desktop/combined_pfams_GH.csv") 
GH_df <- GH_hits %>% 
  separate(info, sep="\t", into = c("gene_oid",	"gene_length", "query_start",	"query_end",	"subj_start",	"subj_end", "evalue",	"bit_score",	"pfam_id",	"pfam_name",	"pfam_length"), remove=TRUE, convert = FALSE)%>%
  data.frame
GH_unique_pfams <- unique(GH_df$pfam_id)
View(GH_unique_pfams)
```
We see that for glycosyltransferase, there are 38 unique glycosyltransferase (left) and 79 unique glycosyl hydrolase (right) hits/families that were identified. (Only first few shown)

![](https://i.imgur.com/7O7DeWf.png)     ![](https://i.imgur.com/CPIylLE.png)

### Step 7: Acquiring all PFAM hits

Similar to how the four previous genes - sulfotransferase, sulfatase, FGE, and polysaccharide export genes -  were gathered and placed into individual fast files and blasted using NCBI BLAST+, each of these PFAMs were acquired, placed in their own fast file, and then all moved into one folder each for the glycosyltransferase and glycosyl hydrolase PFAM queries.

### Step 8: Moving PFAM queries over to HPC and preparing to blast

```
scp -r path/to/local/computer/PFAM/HITS/folder path/to/HPC/folder/where/Blasting/will/happen
```

### Step 9: Preparing for BLASTing

Within this folder on the HPC, you should now have a list of each fasta file that contains the PFAM hits. In the glycosyltransferase directoy, there should be 38 fasta files and within the glycosyl hydrolase directory, there should be 79 fasta files. In each deirectory, execute the following:

```
ls > GTs.txt
ls > GHs.txt
```
Within the glycosyltransferase and glycosyl hydrolase directories, respectively. 

Open each of these text files and delete the line that includes the file name (either `GTs.txt` or `GHs.txt`).

Next, move all contents of each glycosyltransferase and glycosyl hydrolase directory into the directory that contains all of the MAG BLAST databases that was created in the earlier BLAST workflow.

```
mv * ../ #for me, I moved it up one directory because I placed the glycosyltransferase and glycosyl hydrolase PFAM hit directories into the direcotry that contained all of the MAGs and MAG blast databases. 
```

### Step 10: BLASTing

Now you should be sitting in a directory that contains the following:
* MAG BLAST databases
* mags.txt file
* GTs.txt file
* GHs.txt file
* all of your glycosyltransferase fasta files
* all of your glycosyl hydrolase fasta files

Within this directory, create the following script. It is very similar to the one used for the other 4 genes but includes a for loop since there are many (38 and 79) genes to go through. You will want to execute this twice - once for all the glycosyltransferase pfams and once for the glycosyl hydrolase pfams.

Starting with glycosyltransferase..
```
#!/bin/bash
#SBATCH -p your/partition/here      
#SBATCH -N 1                          
#SBATCH -n 1            
#SBATCH --mem=10G                     
#SBATCH --time=1:00:00                
#SBATCH -J blasts              
#SBATCH --output=blasts_%a.out     
#SBATCH --error=blasts__%a.err      
#SBATCH --array=1-84

module add BLAST here 

FILENAME=$(awk "NR==$SLURM_ARRAY_TASK_ID" mags.txt) 

for i in $(cat GTs.txt)
do
tblastn -query $i -db $FILENAME -out $FILENAME+"hits_$i" -outfmt "7 std"
done
```
Now move all hits into a new folder created for glycosyltransferase hits.

```
mkdir GT_hits #makes new directory

mv *+hits_* GT_hits/
```

Do this once more with glycosyl hydrolase hits:

```
#!/bin/bash
#SBATCH -p your/partition/here      
#SBATCH -N 1                          
#SBATCH -n 1            
#SBATCH --mem=10G                     
#SBATCH --time=1:00:00                
#SBATCH -J blasts              
#SBATCH --output=blasts_%a.out     
#SBATCH --error=blasts__%a.err      
#SBATCH --array=1-84

module add BLAST here 

FILENAME=$(awk "NR==$SLURM_ARRAY_TASK_ID" mags.txt) 

for i in $(cat GHs.txt)
do
tblastn -query $i -db $FILENAME -out $FILENAME+"hits_$i" -outfmt "7 std"
done
```


```
mkdir GH_hits #makes new directory

mv *+hits_* GH_hits/
```

### Step 11: More file manipulation

Now navigate to your `GT_hits/` directory and remove the first 5 lines and the last line from each file (just like before!):

```
for f in *; do
    tail -n +6 "$f" > "${f}".tmp && mv "${f}".tmp "$f"
    echo "Processing $f"
done
```
And the last line.
```
for f in *; do
    sed '$d' "$f" > "${f}".tmp && mv "${f}".tmp "$f"
    echo "Processing $f"
done
```
And once again, combine the files!
```
awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' OFS=, * > GT_combined_hits.csv
```

Do this for both the glycosyltransferase and glycosyl hydrolase hits.

### Step 12: Move over files

Again, open your local terminal and move these two files to your computer.

* GT_combined_hits.csv
* GH_combined_hits.csv

```
scp path/to/HPC/file/GT_combined_hits.csv path/to/local/computer/
scp path/to/HPC/file/GH_combined_hits.csv path/to/local/computer/
```

### Step 13: Final file edits 

Once again, I opened excel and modified the files by adding the column names "MAG" and "info", and changing the `file_name` in line 1. 

![](https://i.imgur.com/jlF6WtY.png)

These tables were then imported into R.
 
## :memo: Data manipulation and visulaization in R

#### View ==this== README for detailed explanation of how I did this in R.






