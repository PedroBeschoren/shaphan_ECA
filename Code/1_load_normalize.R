
#load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)


# increases memory limit used by R
memory.limit(size = 350000)


##### Loading 16S  data ##########

# loads the OTU table, immediatetly saving it as a phyloseq OTU table (lighter than a df)
raw_bac_otutab<-otu_table(object = read.table(file = "./Data/16S_ECA_seqtab_final.txt", 
                                              header = TRUE,
                                              row.names = 1, # first column has row names (ASV names)
                                              check.names = FALSE), # prevents "X" to be added to column names, such as X49_16S,
                          taxa_are_rows = TRUE)


#updated tax
sklearn_bac_taxtab<-read.table(file = "./Data/16S_ECA_taxonomy.tsv", 
                               header = TRUE,
                               sep = "\t",
                               row.names = 1, # first column has row names (ASV names)
                               check.names = FALSE) # prevents "X" to be added to column names, such as X49_16S,

# adjusts number and name of columns
sklearn_bac_taxtab<- separate(data = sklearn_bac_taxtab,
                              col = Taxon,
                              into =c("Kingdom", 
                                      "Phylum", 
                                      "Class", 
                                      "Order", 
                                      "Family",
                                      "Genus", 
                                      "Species"),
                              sep = ";")

# change from d__Bacteria to k__bacteria, mathcing fungi dataset
sklearn_bac_taxtab$Kingdom<-gsub("d__", "k__", sklearn_bac_taxtab$Kingdom)

# saves taxa as phyloseq object
sklearn_bac_taxtab<-tax_table(object = as.matrix(sklearn_bac_taxtab))

#change name and remove old object
raw_bac_taxtab<-sklearn_bac_taxtab
rm(sklearn_bac_taxtab)

# loads the representative sequences table, immediatetly saving it as a phyloseq OTU table (lighter than a df)
raw_bac_refseq<-refseq(physeq = Biostrings::readDNAStringSet(filepath = "./Data/16S_ECA_repset.fasta", use.names = TRUE)) 
taxa_names(raw_bac_refseq)<-gsub(" .*", "", taxa_names(raw_bac_refseq)) # drops taxonomy from ASV names

# loads the mapping tile, immediatetly saving it as a phyloseq OTU table (lighter than a df)
raw_bac_metadata<-sample_data(object = read.table(file = "./Data/16S_ECA_Mapping_file.txt", 
                                                  header = TRUE,
                                                  sep = "\t",
                                                  row.names = 1, # first column has row names (ASV names)
                                                  check.names = FALSE)) # prevents "X" to be added to column names, such as X49_16S,

# build the main 16S phyloseq object by putting all these phyloseq-class objects (otu table, tax table, ref seq, metadta) into a single ps
raw_bac_ps<-merge_phyloseq(raw_bac_otutab,
                           raw_bac_taxtab,
                           raw_bac_refseq,
                           raw_bac_metadata)

# remove old objects to reduce memory use
rm(raw_bac_otutab,raw_bac_taxtab,raw_bac_refseq,raw_bac_metadata)

# change names from ASV to fASV, so they can be distinguished from bacterial ASVs
taxa_names(raw_bac_ps)<-paste("b", taxa_names(raw_bac_ps), sep = "")
tax_table(raw_bac_ps)<-gsub(" ", "", tax_table(raw_bac_ps)) # drops a space charather from taxa names. I don't know how that charather ended up in there

# let's check the imported objects. Often errors will arise from typos when filling up the data sheets 
otu_table(raw_bac_ps)[1:10,1:10]
sample_data(raw_bac_ps)[1:10,1:10]
tax_table(raw_bac_ps)[1:10,1:6]
refseq(raw_bac_ps)



# run garbage collection after creating large objects
gc()
