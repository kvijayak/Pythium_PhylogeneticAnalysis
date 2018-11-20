## Assignment 2, Part C

## Introduction and Hypotiesis ------

## My focus is on a taxa called Pythium.  Pythium are a collection of fungi that are most known for their tendency to cause great damage to agriculture fields and greenhouse soils.  While the vast majority of Pythium species are pathogenic to plants, there are a few that are atually pathogenic to other fungi (including other pythium).  This small subset of pythium are considered "mycoparasites", and are accordingly used by farmers as a bioagent in soil treatment.  One well-known mycoparasite species is pythium oligandrum, which will be the focus of this phylogenetic analysis.  The goal of this analysis is to create a tree using the pythium dataset from the BOLD database.  I'm going to "track" the p.oligandrum species throughout the alignment, clustering, and dendrogram-construction.  Once I have the dendrogram, I'd like to see which species are clustered with my species of interest.  My hypothesis is that if additional pythium species  are clustered togeather with p.oligandrum, then they too might also possess mycoparasitic properties that could qualify them as bioagents for soil health.  To do this, I will download all of the pythium data available in the BOLD database. I'll clean and process the data to create a dendrogram.  Then, I'll see which other species were clustered with p.oligandrum.  Finally, I'll do a literature search to see if these species do in-fact possess mycoparasitic properties. 

## It's important to note that my dendrogram will depict the evolutionary relationships between pythium species based on a short standarized fragment of DNA (in this case the ITS marker).  If other OTUs are clsutered with the mycoparasitic species p.oligandrum, this doesn't necessarily make them mycoparasites themselves.  This would only be true if the "parasite phenotype" (which researchers believe is caused by the structure of their spiny oogonia) was a direct result of the transcripts/proteins created by the ITS mitochondrial gene (which is very likely not the case).




##SOURCES
## "Comparison of the mycoparasites Pythium periplocum, P. acanthicum and P. oligandrum" (Riberio et. al. 1995)
## "A new mycoparasite, Pythium lycopersicum, isolated in Isparta, Turkey: morphology, molecular characteristics, and its antagonism with phytopathogenic fungi" (karaca, 2008)
## "A new species of Pythium with ¢lamentous sporangia having pectinolytic activities, isolated in the Burgundy region of France" (Paul, 2001)
## "Mycoparasitism by Pythium oligandrum and P. acanthicum" (Deacon, 1978)





## Import data and Load Packages------

#install.packages("ape")
library(ape)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("plotrix")
library(plotrix)

#source("https://bioconductor.org/biocLite.R")
#biocLite("DECIPHER")
library(DECIPHER)
#biocLite("Biostrings")
library(Biostrings)


## BOLD Data tsv File
pythium <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Pythium&format=tsv")



## I'm going to use the markercode with the most data values
## 1249 ITS markers, 885 for COI
markerCount <- pythium %>%
  group_by(markercode) %>%
  summarize(n = length(markercode)) %>%
  arrange(desc(n)) %>%
  print()


## ITS Analysis------

## Because there are slightly more records for the ITS markercode, I'll filter for database entries with an ITS markercode sequence.  Also,  I'll remove all data without a bin_uri or species name

pythiumITS <- pythium %>%
  filter(markercode == "ITS") %>%
  filter(str_detect(nucleotides, "[ATCG]")) %>%
  filter(str_detect(bin_uri, ".")) %>%
  filter(str_detect(species_name, ".")) 

## Because the bin_uri is considered a more rigourous proxy for the species_name, I'd like to focus on the bin_uri when I build my dendrogram.  However, the basis of my analysis is to track a particular group of species (p.oligandrum and it's potentially related taxa) based on the species names. To relate species names to bins, I'm going to create a dataframe that includes the bin_uri in one column and it's associated species_name in a second column.

bin_to_species.name <- unique(pythiumITS[,c(8,22)]) 

## I also removed duplicate rows with the unique() function.  Now, I can relate any bin to any species name, or vice versa. Let's search for the bin associated with our species of interest, p.oligandrum: 

bin_to_species.name[which(str_detect(bin_to_species.name$species_name, "Pythium oligandrum")),]


## the bin associated with p.oligandrum is BOLD:AA06345.  We can now use this identifier to track our species of interest throughout the remainder of the analysis.  For our alignment and subsequent clustering, we'll take one record (at random) from each unique bin_uri, using the sample_n() function from the dplyr package. I'm also using set.seed in hopes that when this script is re-run, the same "random" collection of records will be subsetted.  If the script is re-run and sample_n takes a different sample of records, the remainder of the downtream analysis may not be reproducable

set.seed(86)
pythiumITS_sample <- pythiumITS %>%
  group_by(bin_uri) %>%
  sample_n(1)


## Here, we have 149 bin_uri's and 97 species names
length(unique(pythiumITS_sample$bin_uri))
length(unique(pythiumITS_sample$species_name))

## We have an outlier in our new dataframe.  The BOLD:ACE4105 bin is highly dissimilar from the other 148 sequences in our set (based on a dendrogram I created not shown in this script).  When I BLASTed this sequence, I actually found that although it is genetic information from the ITS region, it does not belong to any species of pythium.  BLAST returned high-scoring hits for other fungal species, such as Sarocladium kiliense and Nectria mauritiicola.
nucleotides <- DNAStringSet(pythiumITS_sample$nucleotides)
Pythium_musclealignment <- DNAStringSet(muscle::muscle(nucleotides, quiet = FALSE, maxiters = 4, diags = FALSE))
Pythium_musclealignment

dnaBin.pythium <- as.DNAbin(Pythium_musclealignment)

distanceMatrix <- dist.dna(dnaBin.pythium, model = "TN93", as.matrix = TRUE, 
                           pairwise.deletion = TRUE)

Dendogram.pythium <- IdClusters(distanceMatrix,
                                method = "UPGMA",
                                cutoff= 0.03,
                                showPlot = TRUE,
                                type = "dendrogram",
                                verbose = TRUE)
mean(unlist(lapply(Pythium_musclealignment, str_count, ("-"))))
min(unlist(lapply(Pythium_musclealignment, str_count, ("-"))))
which.min(lapply(Pythium_musclealignment, str_count, ("-")))

max(unlist(lapply(Pythium_musclealignment[-129], str_count, ("-"))))
which.max(lapply(Pythium_musclealignment[-129], str_count, ("-")))


pythiumITS_sample <- pythiumITS_sample[-129,]

## I'm going to use DECIPHER to align the sequences in our pythium sample. I chose DECIPHER for a few  reasons: first, it has a very nice visualization tool for aligned sequences.  Also, it's IDClusters() function can cluster data as well as create dendrograms that can be touched-up for better visualization.  This data set isn't particularily large, but DECIPHER also alows users to create a web-based database when running analysis, so that large and unwieldy files need not be loaded into the workspace.

## Before we can perfrom the msa, our nucleotides must be converted into an XStringSet format, and DECIPHER requires the gaps be removed:

pythiumITS_sample$nucleotides <- DNAStringSet(pythiumITS_sample$nucleotides)
pythiumITS_msa <- AlignSeqs(RemoveGaps(pythiumITS_sample$nucleotides))


## Let's change the names of our allignment to the bin_uri for easier identification.  Then we can view our alignment
names(pythiumITS_msa) <- pythiumITS_sample$bin_uri
BrowseSeqs(pythiumITS_msa)

## Now, let's cluster our reads and build a dendrogram (again using DECIPHER's framework for clustering).  First the data must be converted to an object of class "DNAbin"

pythiumITS_asDNAbin <- as.DNAbin(pythiumITS_msa)
class(pythiumITS_asDNAbin)

## Now, we'll create a distance matrix using the Tamura and Nei 1993 (TN93) evolutionary model.  The TN93 model stands out from earlier models by giving transitions a higher probability of occurance than transversions.  Further, the TN93 distinguishes between the two types of transitions, where A-G transitions can be assigned a different probability of occurance than T-C transitions. 

pythiumITS_asdistDNA <- dist.dna(pythiumITS_asDNAbin, 
                                 as.matrix = TRUE, 
                                 pairwise.deletion = TRUE,
                                 model = "TN93")

## we'll cluster with DECIPHER's IDClusters() function. I had difficulty choosing an algorithm for clustering my data.  I initially clustered the data using each of the methods (single, complete, UPGMA, WPGMA, and neighbor joining) to see the extent of the differences between them:

my_algorithms <- c("single", "complete", "UPGMA", "WPGMA", "NJ")
graphs <- list()

par(mfrow=c(3,2))
for (i in 1:length(my_algorithms)) {  
  algorithm_test <- as.DNAbin(pythiumITS_msa)
  algorithm_test <- dist.dna(algorithm_test, 
                             as.matrix = TRUE, 
                             pairwise.deletion = TRUE,
                             model = "TN93")
  graphs[[i]] <- IdClusters(algorithm_test,
                            method = my_algorithms[i],
                            cutoff= 0.04,
                            type = "dendrogram")
  
  par(cex=0.1, mar=c(5, 6, 10, 2))
  plot(graphs[[i]][[2]], xlab="",ylab="",main="",axes=FALSE)
  par(cex=.75)
  title(main = my_algorithms[i])
}

dev.copy2pdf(file="Pythium Dendrograms - Different Clustering Algorithms")

## UPGMA and WPGMA (and maybe "complete") produced similar-looking trees, however, there is quite a bit of variability.  This is especially true with the neighbor-joing algorithm, possibly becuase that tree isn't rooted.  I ended up using the UPGMA method for the remainder of the analysis.  One of the limitations of UPGMA is that it considers all lineages evolving at the same rate.  However, since my data include members from a somewhat-specific taxa (pyhthium), I figured this wouldn't be a major problem.  

## Choosing a threshold also required some experimenting.  I wanted to find a balance between a diversity of clusters, while avoiding too many singleton and doubleton OTUs.  For this reason, I chose a threshold level of .04

pythiumITS_cluster <- IdClusters(pythiumITS_asdistDNA,
                                 method = "UPGMA",
                                 cutoff= 0.04,
                                 type = "dendrogram")


## I had a difficult time formating my dendrogram with the plot() function.  I took this code from a fourm I found online.  Basically, the code is creating a blank plot and setting paramaters before drawing the dendrogram onto the page.  Then resets the axis parameters once the dendrogram is created.

par(mfrow=c(1,1))
par(cex=0.5, mar=c(5, 6, 4, 2))
plot(pythiumITS_cluster, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(main="Pythium Phylogenetic Tree (ITS Marker)\nUPGMA; Threshold = .04; TN93 Evolutionary Model")
draw.circle(13.75,0,2.75)
axis(2)

dev.copy2pdf(file="Pythium (ITS Marker) UPGMA")



## My species of interest (p.oligandrum, BOLD:AAO6345) was classified with five other taxa: BOLD:AAO6828, BOLD:AAO6829, BOLD:ACE4754, BOLD:ACE4756, and BOLD:AAO5255.  This clade is circled in black on my dendrogram. I'm going to reference these bins in the tsv file to derrive species names for these taxa

pythiumITS_mycoparasites <- pythium %>%
  filter(str_detect(bin_uri,
      "AAO6345|AAO6828|AAO6829|ACE4754|ACE4756|AAO5255")) %>%
  print()

## The new dataframe I created contains all of the species names associated with the bin_uris that were clustered with the p.oligandrum species when I created the dendrogram.  In total, these represent 7 potential mycoparacytes (including my species of interest) on the basis that they are closely related to the known mycoparasyte p.oligandrum.  They are:

# Pythium acanthicum
# Pythium amasculinum
# Pythium hydnosporum
# Pythium lycopersicum
# Pythium oligandrum
# Pythium ornamentatum
# Pythium periplocum

## Interestngly, if we revisit the graph with the 5 plots of the ITS dendrograms based on different clustering algorithms, we see that regardless of the algorithm used (single, complete, UPGMA, WPGMA, or NJ), these species are clustered togeather (with the exception of complete and NJ where two clusters are created) and within the same clade in every plot.  In each of the 5 plots here, our clade of interest is circled in black:


par(mfrow=c(3,2))
for (i in 1:length(my_algorithms)) {  
  algorithm_test <- as.DNAbin(pythiumITS_msa)
  algorithm_test <- dist.dna(algorithm_test, 
                             as.matrix = TRUE, 
                             pairwise.deletion = TRUE,
                             model = "TN93")
  graphs[[i]] <- IdClusters(algorithm_test,
                            method = my_algorithms[i],
                            cutoff= 0.04,
                            type = "dendrogram")
  
  par(cex=0.1, mar=c(5, 6, 10, 2))
  plot(graphs[[i]][[2]], xlab="",ylab="",main="",axes=FALSE)
  par(cex=.75)
  title(main = my_algorithms[i])
  
  if (my_algorithms[i] == "single")
    draw.circle(49.75,0,2.75)
  if (my_algorithms[i] == "complete")
    draw.circle(14.25,0,2.75)
  if ( my_algorithms[i] == "UPGMA")
    draw.circle(13.75,0,2.75)
  if (my_algorithms[i] == "WPGMA")
    draw.circle(12.25,0,2.75)
  if (my_algorithms[i] == "NJ")
    draw.circle(104.5,0.065,2.25)
  
}

dev.copy2pdf(file="Pythium Dendrograms - Different Clustering Algorithms w/ Clade Circled")


## I'd like to do a quick search into the literature to see if there's been any research into mycoparasitic properties of these species and if they've been used as biocontrol agents in soil.  But first, I'm going to re-run my analysis using the COI markercode instead of the ITS.  I'm doing this to test if there are any other OTUs closely related to p.oligandrum bin_uri that wern't captured by the ITS markercode analysis (potentially because of missing data)








## COI Analysis ------

## Filter for COI (and filter out missing bin_uri & species_names)
pythiumCOI <- pythium %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(bin_uri, ".")) %>%
  filter(str_detect(species_name, ".")) 


## Create a bin-to-species name index
bin_to_species.name2 <- unique(pythiumCOI[,c(8,22)]) 

## Sample one record per bin_uri
set.seed(43)
pythiumCOI_sample <- pythiumCOI %>%
  group_by(bin_uri) %>%
  sample_n(1)

## 169 bins, 93 species names
length(unique(pythiumCOI_sample$bin_uri))
length(unique(pythiumCOI_sample$species_name))

## Align our Sequences
pythiumCOI_sample$nucleotides <- DNAStringSet(pythiumCOI_sample$nucleotides)
pythiumCOI_msa <- AlignSeqs(RemoveGaps(pythiumCOI_sample$nucleotides))
names(pythiumCOI_msa) <- pythiumCOI_sample$bin_uri


## cluster our reads and build a dendrogram

pythiumCOI_asDNAbin <- as.DNAbin(pythiumCOI_msa)
pythiumCOI_asdistDNA <- dist.dna(pythiumCOI_asDNAbin, 
                                 as.matrix = TRUE, 
                                pairwise.deletion = TRUE,
                                 model = "TN93")
pythiumCOI_cluster <- IdClusters(pythiumCOI_asdistDNA,
                                 method = "UPGMA",
                                 cutoff= 0.04,
                                 type = "dendrogram")


## Plot our COI Dendrogram
par(mfrow=c(1,1))
par(cex=0.5, mar=c(5, 6, 4, 2))
plot(pythiumCOI_cluster, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(main="Pythium Phylogenetic Tree (COI Marker)\nUPGMA; Threshold = .04; TN93 Evolutionary Model")
draw.circle(12.75,0,2.75)
axis(2)

dev.copy2pdf(file="Pythium (COI Marker) UPGMA")


## My new dendrogram based on the COI markercode is a bit more crowded beucase we clustered 169 sequences instead of 149.  Importantly though, the dendrogram shows the same 6 bins within the same clade again (again, circled in black).  On one hand, I was hoping to discover more bins associated with mycoparasite p.oligandrum. However, it's nice to see that both analyses returned the same results, indicating that these bins are indeed closely related and not by chance.

## since no new bins were discovered by my COI-based analysis, there will be no new entries added to the dataframe "pythiumITS_mycoparacyte" and the species names contained in that dataframe suggest all of my canadite taxa for mycoparasitism.  Here are the species:

pythiumITS_mycoparasites %>%
  arrange(species_name) %>%
  print()


## Literature Search --------

## Doing a quick literature search, here's what I found:

#1 p. acanthicum: a reported mycoparacyte with varying degrees of antaganosim towatds other fungi species.  p. acanthicum is a well-studied, oftentimes with p.oligandrum

#2 p. amasculinum: a few studies indicate that p.amasculinum have mycoparasitic properties.  Researchers in Turkey found that this species infects fungal pests B. cinerea and R. solani.  However, very little research exists on p.amasculinum, and studies that do exist tend to group it with it's well-known relitave p.oligancrum

#3 p. hydnosporum: Uniprot classifies p.hydnosporum as a mycoparasitic fungus.  However, research specifically devoted to p.hydnosporum as a mycoparasite is lacking.  Most research on this species is focused on it's closely related taxa p.oligandrum

#4 p. lycopersicum: more recent research (2008) has identified p. lycopersicum as a mycoparasite, stating that it shows a pronounced antagonism towards other fungal species (the researhcers actually mentioned how the ITS region is highly similar to that of p.oligandrum's)

#5 p. oligandrum: the central species in this analysis, p.oligandrum is the best-researched and most effective known mycoparacyte of the pythium taxa.  Samples of p.oligandrum have been used (with success) as biocontrol agents to antagonize soil-pathogens that migh otherwise infect plants/crops.

#6 p. ornamentatum: very little research exists on p.ornamentatium.  I haven't found any information linking p.ornamentatum and mycoparasitic properties; I generally haven't found any useful information on this species at all, although it is commonly inlcuded in databases like MycoBank, the BOLD Database, and  The Global Catalogue of Microorganisms.

#7 p. periplocum: a known mycoparacyte, is often studied alongside p. oligandrum in mycoparacytic studies.  this taxa is not as antagonistic or common as the p. oligandrum species, but is often cited as a known mycoparacyte and has been used with some success in agriculture


## 4 of these 7 species have been well-researched and their mycoparacytic properties have been validated.  Two of them (p.hydnosporum and p.amasculinum) have alledged mycoparasitic properties, but research is scarce, and their mycoparasitic abilities may be only attributed to the fact that they are highly similar to the p.oligandrum species.  I couldn't find any research on p.ornamentatum.    




## Conclusion -------

## to conclude my script, I wanted to run one more analysis, this time includeing all records from the pythium data that contain a bin_uri.  I thought that including all records with bin_uri information would greatly expand the scope of pythium records, and possibly link other OTUs to the p.oligandrum cluster.  However, when I filtered for all unique bins, I found that only 171 existed in the entire BOLD database.  Considering that my COI analysis contained 169 unique bins, I nearly explored all 171 bins available in the database with my COI analysis, so another analysis including the additional 2 bin_uri's is liklely unnecessary.   

## In this script, I discovered 3 additional pythium species that are known to exhibit a mycoparacitic phenotype, and an additional 3 with possible mycoparacitic properties.  To answer my origional question: Yes, focusing a known mycoparasite, p.oligandrum, I found 6 additional pythium species with a similar mycoparasitic capabilities.  As I said earlier however, It's not necessarily the ITS and COI-5P sequences that "give" these species their mycoparasitic phenotypes.  An interesting next step would be examine what exactly produces this phenotype, which would likely require an analysis of the pythium genome on a genome-wide scale.


