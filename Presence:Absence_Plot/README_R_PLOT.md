# README
#### Here you can find my workflow in R that I used to create the "presence/absence" plot of genes involved in the production and degradation of sulfated polysaccharides in these Shark Bay pustular microbial mat MAGs. 

## :computer: Getting started using R
 
>Note: I am using R version 3.6.0. 

I gathered gene count data for 6 genes I was interested in that are involved in the processes of production and degradation of sulfated polysaccharides. <b>Bolded</b> are the genes whose PFAM IDs were initially gathered from JGI IMG and the used to blast. All other genes had PFAM hits directly taken from NCBI BLAST. 

| Gene          | Production/ Degradation               |
| ----------------- |:----------------------- |
| Sulfotransferase      | Production   |
| <b>Glycosyltransferase</b> | Production    |
| Polysaccharide exporter         | Production     |
| Sulfatase       | Degradation    | 
| <b>Glycosyl hydrolase</b>   | Degradation |
| Formylglycine-generating enzyme (FGE)  | Degradation |

 The workflow for creating the different files (combined_hits.csv, GT_combined_hits.csv, and GH_combined_hits.csv) we are importing into R can be found [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/Presence:Absence_Plot/README_BLAST.md). 

### Step 1: Importing libraries and files

```
library(reshape)
library(ggplot2)
library(readr)
library(tidyselect)
library(tidyverse)
library(dplyr)

NCBI_gene_hits <- read.csv("~/Desktop/combined_hits.csv") 
GT_hits <- read.csv("~/Desktop/GT_combined_hits.csv") 
GH_hits <- read.csv("~/Desktop/GH_combined_hits.csv") 
```

The combined_hits.csv file can be found [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/Presence:Absence_Plot/combined_hits.csv). The GT_combined_hits.csv can be found [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/Presence:Absence_Plot/GT_combined_hits.csv). The GH_combined_hits.csv can be found [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/Presence:Absence_Plot/GH_combined_hits.csv). 

Now that we have imported our files of both of our NCBI and JGI hits, we can take a look at them by executing the following lines:

```
View(NCBI_gene_hits)
View(GT_hits)
View(GH_hits)
```

### Step 2: Manipulation, manipulation, and more manipulation

## NCBI gene hits (Sulfotransferase, sulfatase, FGE, and polysaccharide exporter genes)

*First focusing on the "NCBI gene hits"...*

MAG 78 was removed from the entire dataset due to inconsistencies with taxonomy assignment and phylogenetic placement. Thus, it was removed from this data.

```
NCBI_gene_hits_no_MAG_78 = NCBI_gene_hits[-(8491:8586),] 
```
Now, we will split the 'info' column into columns based on tab separation and assign names to each of these columns AND split the 'MAG' column into 'MAG' name and 'hit' specificity. 
```
df <- NCBI_gene_hits_no_MAG_78 %>% 
  separate(info, sep="\t", into = c("query acc.ver", "subject acc.ver", "% identity",
                                    "alignment length", "mismatches", "gap opens", 
                                    "q. start", "q. end", "s. start", "s. end", "evalue_orig",
                                    "bit score"), remove=TRUE, convert = FALSE)%>%
  separate(MAG, sep = "[+]", into = c("MAG", "hit", remove=TRUE, convert = FALSE)) %>% 
  data.frame
```

If this doesn't make sense you can take a look at it with the ```View()``` function.

We went from our table looking like this:

![](https://i.imgur.com/hAAtebI.png)

to this:

![](https://i.imgur.com/oL3uIuT.png)

Next, I needed to convert the e-values (which were being read as characters) to numeric. I did so by creating a new dataframe with the old e-values and converting them to numeric and then binding the old and new dataframes together.

>For example, 
(numeric values): 9e^-10 < 2.01
(character values): 9e^-10 > 2.01 because "e^-10" is not perceived by R to be numeric so it might as well read 9 > 2.01

```
eval_df <- as.numeric(as.character(df$evalue_orig))
eval_df2 <- as.data.frame(eval_df)
colnames(eval_df2)[1] <- "evalue"
newest_df <- cbind(df, eval_df2$evalue)
colnames(newest_df)[17] <- "evalue"
View(newest_df)
```

For my purposes, I was mainly interested in having a table with only the MAG ID (ex. MAG 1, MAG 24, etc.), the gene hit (ex. sulfatase, sulfotransferase, etc.), percent identity, and e-values. Therefore, I removed all other columns.

```
clean_df <- newest_df[ -c(3:6,8:16) ]
```
![](https://i.imgur.com/C1uwA4j.png)

I also changed the second column name to "gene":
```
colnames(clean_df)[2] <- "gene"
```

This leaves us with a table that looks like this:
![](https://i.imgur.com/f2gAEff.png)

## JGI gene hits (glycosyltransferase and glycosyl hydrolase)

These same steps were taken with the glycosyltransferase and glycosyl hydrolase data.

MAG 78 was removed from the dataset:
```
GT_hits_no_MAG_78 = GT_hits[-(16882:17115),]
GH_hits_no_MAG_78 = GH_hits[-(28179:28691),]
```
Assigning column names by splitting current headers:
```
GT_df <- GT_hits_no_MAG_78 %>% 
  separate(info, sep="\t", into = c("query acc.ver", "subject acc.ver", "% identity",
                                    "alignment length", "mismatches", "gap opens", 
                                    "q. start", "q. end", "s. start", "s. end", "evalue_orig",
                                    "bit score"), remove=TRUE, convert = FALSE)%>%
  separate(MAG, sep = "[+]", into = c("MAG", "hit", remove=TRUE, convert = FALSE)) %>% 
  data.frame
  
################################################################################  

GH_df <- GH_hits_no_MAG_78 %>% 
  separate(info, sep="\t", into = c("query acc.ver", "subject acc.ver", "% identity",
                                    "alignment length", "mismatches", "gap opens", 
                                    "q. start", "q. end", "s. start", "s. end", "evalue_orig",
                                    "bit score"), remove=TRUE, convert = FALSE)%>%
  separate(MAG, sep = "[+]", into = c("MAG", "hit", remove=TRUE, convert = FALSE)) %>% 
  data.frame
```

"Restoring" evalue column values as numeric and simplifying the dataframe
```
GT_eval_df <- as.numeric(as.character(GT_df$evalue_orig))
GT_eval_df2 <- as.data.frame(GT_eval_df)
colnames(GT_eval_df2)[1] <- "evalue"
GT_newest_df <- cbind(GT_df, GT_eval_df2$evalue)
colnames(GT_newest_df)[17] <- "evalue"
View(GT_newest_df)

GT_clean_df <- GT_newest_df[ -c(3:6,8:16) ]
GT_clean_df <- GT_clean_df[ -c(2) ]

################################################################################  

GH_eval_df <- as.numeric(as.character(GH_df$evalue_orig))
GH_eval_df2 <- as.data.frame(GH_eval_df)
colnames(GH_eval_df2)[1] <- "evalue"
GH_newest_df <- cbind(GH_df, GH_eval_df2$evalue)
colnames(GH_newest_df)[17] <- "evalue"
View(GH_newest_df)

GH_clean_df <- GH_newest_df[ -c(3:6,8:16) ]
GH_clean_df <- GH_clean_df[ -c(2) ]
```
Each table will look like the following:
![](https://i.imgur.com/RJpSRij.png)

This nest step is unique to these two tables. A new column was added to each table with either the glycosyltransferase or glycosyl hydrolase labels. This was done because ultimately, regardless of what the pfam subfamily was, each hit was a glycosyltransferase or glycosyl hydrolase gene, and that is what we want to plot (not each individual family). This was done so:
```
GT_clean_df$gene <- "glycosyltransferase"
GH_clean_df$gene <- "glycosyl hydrolase"
```
![](https://i.imgur.com/60W6CC7.png)


![](https://i.imgur.com/Oe1Vucm.png)

Now we combine both dataframes:

```
combined_GT_GH_df <- dplyr::bind_rows(GT_clean_df, GH_clean_df)
View(combined_GT_GH_df)
```

![](https://i.imgur.com/YNvVct6.png)

### Step 3: Combining all dataframes and creating presence/absence plot
## Combining
Now we combine all dataframes regardless of  
```
combined_all_genes_df <- dplyr::bind_rows(clean_df, combined_GT_GH_df)
View(combined_all_genes_df)
```
And then clean up the MAG IDs...
```
clean_combined_all_genes_df <- combined_all_genes_df %>% 
  mutate(MAG = gsub('SB_biofilm_|\\.fa', '', MAG)) %>% 
  data.frame
```
## Presence/absence

Now, we are creating a new column called 'presence', where we are basically saying, if the column 'evalue' has a value of 1e^-10 or smaller AND the column 'X..identity' (% identity) is greater than or equal to 30% for any given gene hit, we are assigning a value of '1' (pretty much for TRUE). If a gene hit has an e-value larger than 1e^-10 and a percent identity smaller than 30%, OR if only one of these is true, then a value of '0' will be applied for that hit within this new column called presence. 

```
clean_combined_all_genes_df$presence = ifelse(clean_combined_all_genes_df$evalue<=0.0000000001 & clean_combined_all_genes_df$X..identity >= 30,1,0)
```

An example of these results can be seen below.

![](https://i.imgur.com/X149zCP.png)

Now we use the summarize function to add the values in the presence column and group by each unique MAG and gene (group_by function) and we put those summed values into a new column called "present". The number in the "present" column tells us how many times there was a hit for that gene in that MAG that met the e-value and percent identity cutt-off values. 
```
summarized_df <- clean_combined_all_genes_df %>% 
  group_by(MAG, gene) %>%
  summarize(present= sum(presence)) %>% 
  data.frame
View(summarized_df)
```
![](https://i.imgur.com/dtywGaj.png)

Now since we are just interested in presence and absence, we convert those values to 0 and 1 which reflect absence and presence, respectively.
```
summarized_df$pres_abs = ifelse(summarized_df$present >= 1,1,0)
```
![](https://i.imgur.com/POqNL89.png)

We can then remove the presence row:
```
clean_summarized_df <- summarized_df[ -c(3) ]
```
![](https://i.imgur.com/qkPVWTB.png)

Next, we add taxonomy to the table by first importing the text file which can be found [here](https://github.com/emilieskoog/SharkBay2020-analysis/blob/master/Presence:Absence_Plot/taxa.txt).
```
taxa_info <- read.delim("~/Desktop/taxa.txt") 
```
And then modifying it by breaking up the "info" header and specifying columns names according to taxonomic rank.

```
taxa_info <- taxa_info %>% 
  separate(info, sep=";", into = c("domain","phylum", "class", "order", "family", "genus", "species"), remove=TRUE, convert = FALSE)%>%
  data.frame
```
![](https://i.imgur.com/39Jb77D.png)

>Notice this table and the clean_summarized_df table both have the column "MAG" with the same labeling scheme.


Then we merge it with our dataframe with our presence and absence information. 
```
df_with_taxa <- merge(taxa_info, clean_summarized_df, by="MAG")
```
and remove the "__p" in front of each phylum name (will become handy later)
```
cleaned_df_with_taxa <- df_with_taxa %>% 
  mutate(phylum = gsub('p__', '', phylum)) %>% 
  data.frame
View(cleaned_df_with_taxa)
```
![](https://i.imgur.com/CQK1dH5.png)

I wanted my table to be summarized at the order level for each gene. I also wanted to average the amount of MAGs in each order that did contain at least one gene (present = 1), and I assigned that value to a new column called "avg".
```
making_it_happen <- cleaned_df_with_taxa %>% 
  group_by(phylum, class, order, gene) %>%
  summarize(avg= mean(pres_abs)) %>%
  data.frame
```
I created a new column called "class_rder" that combined the "class" and "order" columns for easy grouping by higher level of taxonomy during plotting.

```
making_it_happen$class_order <- paste(making_it_happen$class,making_it_happen$order)
```
I then wanted to create a new column called "metabolic_potential" where certain genes are defined as either production or degradation genes. 
```
making_it_happen <- making_it_happen %>% 
  mutate(metabolic_potential = (gene == "hits_sulfatase" | gene == "hits_FGE" | gene == "glycosyl hydrolase"), metabolic_potential = ifelse(metabolic_potential == TRUE, "degradation", "production")) %>%
  group_by(phylum, class, order, gene, metabolic_potential) %>%
  data.frame
```
This code says, if a gene is annotated as "hits_sulfatase" or "hits_FGE" or "glycosyl hydrolase", assign it that "value" of "degradation" and if it is not, then assign it as "production." We are now adding this value of production or degradation as something to group these genes by in addition to the taxonomy and gene.

I then wanted to change the order that these labels would appear in (it would have been alphabetical, but I wanted production on the left and degradation on the right), and I wanted to capitalize them. So here I use the "factor" function to order them and the "labels" function to write exactly how I wanted them to read.
```
making_it_happen$metabolic_potential_order <- factor(making_it_happen$metabolic_potential, levels = c("production", "degradation"), labels = c("Production", "Degradation"))
```
I also wanted to rearrange the order in which the phyla appearred so I did the following:
```
making_it_happen$phylum <- factor(making_it_happen$phylum, levels = c("Hydrogenedentota", "Cyanobacteriota", "BRC1", "Proteobacteria", "Planctomycetota", "Bacteroidota", "Myxococcota", "Chloroflexota", "Verrucomicrobiota"))
```
I also renamed the six genes and rearranged the order of them as well (notice I assigned this gene_order value which will be important to call this when when plotting).
```
making_it_happen$gene_order <- factor(making_it_happen$gene, levels = c("hits_sulfotransferase3", "glycosyltransferase", "hits_poly_ex", "hits_sulfatase", "hits_FGE", "glycosyl hydrolase" ), labels = c("Sulfotransferase", "Glycosyltransferase", "Polysaccharide exporter", "Sulfatase", "FGE", "Glycosyl hydrolase"))
```
And finally, I wanted to change the decimal value of the "avg" column to a percentage value. Again, notice how I assigned it a new value of "avg_percentage."
```
making_it_happen$avg_percentage <- making_it_happen$avg*100
```
AND NOW PLOT!

```
plot <- ggplot(data = making_it_happen) +
  geom_tile(mapping = aes(x = gene_order, y = class_order, fill= avg_percentage, colour=I("black"))) +
  labs(x = "Percentage of MAGs with at least one copy of gene") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_gradientn(colours = c( "white", "#95D840FF", "#2D708EFF"), values = c(0, 0.1, 1)) +
  facet_grid(phylum ~ metabolic_potential_order, scales = "free", space = "free", drop = FALSE) +
  theme(strip.text.y = element_text(angle = 0)) + theme(strip.background = element_rect(fill=NA))  + theme(panel.spacing.y=unit(0.1, "lines"))
plot
```
If all went well (and with a little illustrator magic), you should have something that looks like this!:tada: :confetti_ball: :tada:

![](https://i.imgur.com/mOBqn0U.png)
