

#**********************************************#
####### remove_Chloroplast_Mitochondria ########
#**********************************************#

# This function was written by Pedro Costa in Jul/2021
# it will remove identified o__Chloroplast and f__Mitochondria taxa from the dataset
# it's input in'ts a phyloseq object with these plant sequences
# the output will not contain these platn sequences
# there are many options and plots you could make to check how much plant DNA you ahve amplified - something we won't explore ehre


##### Chloroplast and mitochondrial DNA removal ##########
remove_Chloroplast_Mitochondria<- function(physeq_object){
  
  #Removes Chloroplast
  #ATTENTION: if you just do >subset_taxa(physeq_object, Rank4!="o__Chloroplast") ; you will also remove NAs in the identification. thus you have to turn the ASV id into a factor before removing them, check https://hypocolypse.github.io/16s-data-prep.html
  
  # generate a df with Chloroplast ASVs
  CH1 <- subset_taxa(physeq_object, Order == "o__Chloroplast" | Order == "Chloroplast") # get all Chloroplasts...
  CH1 <-  as(tax_table(CH1), "matrix")
  CH1 <- row.names(CH1) # get ASV IDs...
  CH1df <- as.factor(CH1) # set IDs as factors
  goodTaxa <- setdiff(taxa_names(physeq_object), CH1df) #define taxa you should keep
  ps_no_chloro <- prune_taxa(goodTaxa, physeq_object) # your new physeq object is now chloroplast-free, but retains NA in identification
  
  
  
  #Removes Mitochondria
  #ATTENTION: if you just do >subset_taxa(physeq_chloro_mito, Rank5!="f__Mitochondria") ; you will also remove NAs in the identification. thus you have to turn the ASV id into a factor before removing them, check https://hypocolypse.github.io/16s-data-prep.html
  MT1 <- subset_taxa(physeq_object, Family == "f__Mitochondria" | Family == "Mitochondria")
  MT1 <-  as(tax_table(MT1), "matrix")
  MT1 <- row.names(MT1)
  MT1df <- as.factor(MT1)
  goodTaxa <- setdiff(taxa_names(physeq_object), MT1df)
  ps_no_chloro_mito <- prune_taxa(goodTaxa, ps_no_chloro)
  
  #excellent! let's save this a new phyloseq object
  physeq_clean<-ps_no_chloro_mito
  return(physeq_clean)
  
}

# ********************************************** Done! #













#**********************************************#
############ decontaminate_and_plot ############
#**********************************************#

# the input is a phyloseq object with "Blank"samples annotated in the "Stress" metadata column
# The function will perform decontamination, report n of contaminants, and provide a plot of contaminating ASVs
# the output is a list with the decontaminated phyloseq object, along with the plot and the number of true contaminants
##### NOTE: parameters of isContaminant() can be optimized for each dataset


decontaminate_and_plot<-function(phyloseq_object){
  
  require(decontam)
  
  # set negatives as TRUE in a  new column 
  sample_data(phyloseq_object)$is.neg <- sample_data(phyloseq_object)$sample  == "MOCK" 
  
  # add library size to metadata
  sample_data(phyloseq_object)$library_size<-sample_sums(phyloseq_object)
  
  # decontaminate! read details on decontam package. use at least 3 blank samples for minimal decontamination. grouping your blanks in batches (such as sampling blanks, DNA extraction blanks, PCR blanks) can help you remove contaminants according different sources. Try to obtain the DNA concentration from the PCR product to include this variable in the decontamination process.
  decontam_output <- isContaminant(transform_sample_counts(phyloseq_object, function(OTU) OTU/sum(OTU)),
                                   neg="is.neg",
                                   threshold= 0.15)  
  n_contaminants<- table(decontam_output$contaminant) # this shows the number of  contaminates (=TRUE)
  head(decontam_output) # checks decontam output
  
  contaminants <- rownames(subset(decontam_output, contaminant%in%c("TRUE")))
  contaminants #this is a list of otus classified as contaminants
  
  
  # Make phyloseq object of presence-absence (pa) in negative controls and true samples
  physeq_norm_pa <- transform_sample_counts(transform_sample_counts(phyloseq_object, function(OTU) OTU/sum(OTU)),
                                            function(abund) 1*(abund>0))
  
  physeq_norm_pa_neg <- prune_samples(sample_data(physeq_norm_pa)$sample == "MOCK", 
                                      physeq_norm_pa)
  physeq_norm_pa_pos <- prune_samples(sample_data(physeq_norm_pa)$sample != "MOCK", 
                                      physeq_norm_pa)
  
  # Make data.frame of prevalence in positive and negative samples
  physeq_norm_df_pa <- data.frame(pa.pos=taxa_sums(physeq_norm_pa_pos), 
                                  pa.neg=taxa_sums(physeq_norm_pa_neg),
                                  contaminant=decontam_output$contaminant)
  
  decontamination_plot<-ggplot(data=physeq_norm_df_pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
    geom_point() +
    xlab("Prevalence (Negative Controls)") + 
    ylab("Prevalence (True Samples)")+
    ggtitle("Decontamination plot")
  
  
 
  
  #clean physeq object by removing contaminant OTUS, then remove blank sample 
  phyloseq_object <- prune_taxa(!taxa_names(phyloseq_object) %in% contaminants, phyloseq_object) # keep only taxa that are not in the list of contaminants
  #phyloseq_object <- subset_samples(phyloseq_object, colSums(otu_table(phyloseq_object))>0)
  
  #clean physeq object by removing blank samples
  physeq_decontaminated <- subset_samples(phyloseq_object, sample != "MOCK")
  

  
  
  output<-list(physeq_decontaminated,
               decontamination_plot,
               decontam_output,
               n_contaminants)
  
  
  return(output)
}

# ********************************************** Done! #














#**********************************************#
################## load_16S_biom ###############
#**********************************************#

# this function tunrs a biom file into a filtered phyloseq object. it should be used for very large biom files that cannot be loaded with phyloseq::import_biom
# such pre-filtering reduces size of the object, allowing R to work outside the HPC

# the input is  the file path to a biom file. 
# the function will then create a phyloseq object (with feature table, metadata, and taxonomy) out of this biom file
# then it will remove ASVs appearing less than 8 times
# then it runs the custom function backup_and_rename
# then it runs the custom function remove_Chloroplast_Mitochondria
# the output is a phyloseq object without host DNA and wihtout rare taxa


load_16S_biom<- function(file){
  infile <- file
  biom <- read.biom(infile, tree=FALSE, prune =TRUE)
  
  physeq_16S<-phyloseq(otu_table(as.matrix(biom$counts), taxa_are_rows = TRUE),
                       tax_table(biom$taxonomy),
                       sample_data(biom$metadata))
  
  # removes taxa occuring less than 8 times in the dataser, as VSEARCH standard
  otu_table(physeq_16S) <- otu_table(physeq_16S)[which (rowSums(otu_table(physeq_16S)) > 7),] 
  
  # this will adjust the ASV names
  physeq_16S<-backup_and_rename(physeq_16S)
  
  # This will remove plastid and mithocndrial DNA
  physeq_16S<- remove_Chloroplast_Mitochondria(physeq_16S)
  
  # removes undedntified seqeunces and acheaea
  physeq_16S<-subset_taxa(physeq_16S, Domain=="d__Bacteria")
  
  return(physeq_16S)
  
}

# ********************************************** Done! #








#**********************************************#
################### rename_otu #################
#**********************************************#


# from https://rdrr.io/github/cpauvert/psadd/man/rename_otu.html
#' OTU renaming: from hash to rank abundance prefix
#'
#' OTU names when generated with vsearch for example corresponds
#' to sequence SHA1 hash value. Despite being unique, this notation
#' is not human readable.
#' This function transform these hashes into easier OTU names.
#' Ex: `e443b18` to `OTU_1` for the most abundant OTU
#'
#' This renaming also has the advantage to ease subsetting with OTU names
#' because OTU names are the \code{\link[base]{row.names}} of several
#' \code{\link[phyloseq]{phyloseq}} object.
#'
#' @param physeq \code{\link{phyloseq-class}}
#' @param hash2id_file A character giving the filename to write the
#' correspondence table between the hash and the new ID or "" for output to
#' the console.
#' @return physeq \code{\link{phyloseq-class}} with OTU renamed.
#'
#' @importFrom phyloseq taxa_sums taxa_names
#' @export
#' @examples
#' require(phyloseq)
#' data(esophagus)
#' taxa_names(esophagus)[1:6]
#' taxa_names( rename_otu( esophagus, hash2id_file = "" ) )[1:6]
rename_otu<-function(physeq, hash2id_file){
  # OTU Abundance vector
  otu_abd<-taxa_sums(physeq)
  # Get OTU rank after sorting by decreasing abundance
  rank_otus<-names(sort(otu_abd, decreasing = T))
  # Fetch hash otus names
  hash_otus<-taxa_names(physeq)
  # Associate hash to rank abundance prefix
  taxa_names(physeq)<-paste("ASV", match(hash_otus, rank_otus), sep = "_")
  # Create correspondence table for tracability
  correspondence<-data.frame(HashID = hash_otus,
                             RankID = taxa_names(physeq))
  # Write correspondence table to file
  write.csv(x = correspondence,file = hash2id_file,row.names = FALSE)
  return(physeq)
}
#' Convert taxa/ID through a correspondence table
#' @param physeq \code{\link{phyloseq-class}}
#' @param hashTable  A \code{\link{data.frame}} with at least a **UsearchID** and **HashID** columns.
#' @param ... The subsetting expression that should be applied to the
#'  correspondance table \code{hashTable}. This is passed on to \code{\link[base]{subset}}.
#'
#'
#'
#' @return physeq-renamed \code{\link{phyloseq-class}} with OTU renamed.
#' @export
#'
#' @examples
id2hash<-function(physeq, hashTable, ...){
  # Select only rows on interest, assembly method or min overlap for ex.
  # keep only ID columns
  intermediate.df<-subset(hashTable, subset = ..., select = c("UsearchID", "HashID"))
  # Drop unused levels
  intermediate.df<-droplevels(intermediate.df)
  # Replace rownames by USEARCH ID to speed up replacement.
  rownames(intermediate.df)<-intermediate.df$UsearchID
  # Discard unused column
  intermediate.df$UsearchID<-NULL
  
  # Fetch OTU ids
  otu_id<-taxa_names(physeq)
  # Replace USEARCH Id by Hash ID
  taxa_names(physeq)<-as.vector(intermediate.df[ otu_id, ])
  # Return modified phyloseq object
  return(physeq)
}
#' Taxa/OTU abundance in a given sample coupled with rank and taxonomic information
#'
#' @param physeq \code{\link{phyloseq-class}}
#' @param smp A single sample of interest. This integer or sample name will be
#' passed to \code{\link[phyloseq]{get_taxa}}.
#' @param level A character indicating the taxonomic rank of interest. Should be in
#' \code{\link[phyloseq]{rank_names}}. Default is "Species".
#'
#' @return An abundance sorted \code{\link{data.frame}} with the following components
#' \describe{
#' \item{\code{Sample}}{Sample of interest}
#' \item{\code{Abundance}}{Taxa/OTU abundance in the sample of interest}
#' \item{\code{Rank}}{Taxa/OTU rank in the dataset. The higher the rank, the less abundant the OTU.}
#' \item{\code{level}}{Taxonomic rank of interest. Note that the column is named after
#' the provided value from \code{level}.}
#' }
#' @export
#' @seealso \code{\link[phyloseq]{get_taxa}}
#' @examples
get_taxa_nomy<-function(physeq, smp, level="Species"){
  # Fetch OTU/taxa in samples
  taxa_list<-get_taxa(physeq, smp)
  # Subset list to taxa detected
  taxa_list<-taxa_list[ taxa_list > 0 ]
  df<-data.frame(
    Sample = smp,
    Abundance = taxa_list,
    Rank = match( names(taxa_list), names(sort(taxa_sums(physeq),decreasing = T))),
    TaxLevel = as(tax_table(physeq)[names(taxa_list),level],"matrix"))
  # Order by decreasing abundance
  df<-df[order(-df$Abundance),]
  # Include taxonomic level
  colnames(df)<-c('Sample','Abundance','Rank',level)
  return(df)
}

# ********************************************** Done! #










#**********************************************#
################ backup_and_rename #############
#**********************************************#

# this function uses the rename_otus function , makes a backup for dada2 IDs, and then givens new ASV names based on taxa abundance
# then it changes the names of the taxonomy from rank1, rank2 to phylum, class, etc
# it's input is a phyloseq object with taxa names such as 7f09498fb455a333197f691f08b55bc7
# it's output is a phyloseq object with taxa names such as ASV_1, ASV_2, ASV_3, according the ASV abundance
backup_and_rename<-function(phyloseq_object){ 
  
  tax_table(phyloseq_object)<- cbind(tax_table(phyloseq_object), "DADA2_ID"= row.names(tax_table(phyloseq_object))) #makes a backup of the dada2 IDs
  taxa_names(phyloseq_object)<-taxa_names(rename_otu(phyloseq_object, ""))
  colnames(tax_table(phyloseq_object)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "DADA2_ID") # this renames the taxonomy from Rank1, Rank2 to phylum, class etc
  
  return(phyloseq_object)}

# ********************************************** Done! #









#**********************************************#
########## check_n_taxa_NA_percentage ##########
#**********************************************#

# This functions gives a the percentage of NA in the taxonomy annotations for every taxa level
# it has 2 arguments as input:
# ps_object = a phyloseq objects with a taxonomy table
# ntaxa = a numeric value to indicate how many taxa you want to look at (only top 10 or top 10.000)

# check % of unindentified taxa from the top X taxa
check_n_taxa_NA_percentage<-function(ps_object, ntaxa){
  # get top taxa, as defined in ntaxa
  pruned_ps <-prune_taxa(taxa = taxa_names(ps_object)[1:ntaxa],
                         x = ps_object) 
  
  #calculate % of missing taxa
  percentage_NA<-colSums(is.na(pruned_ps@tax_table)/ntaxa(pruned_ps)*100)
  return(percentage_NA)
}

# ********************************************** Done! #
