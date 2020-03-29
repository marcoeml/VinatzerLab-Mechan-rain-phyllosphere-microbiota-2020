## Experimental evidence for rain being a reservoir of tomato phyllosphere microbiota
## Marco E. Mechan Llontop - 2020
## Vinatzer Laboratory at Virginia Tech

## Set working directory
setwd("~/Desktop/Rain_tomato paper")

##Instal libraries
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library("phyloseq") 
library("scales") 
library("ggplot2") 
library("Rcpp") 
library("dplyr")
library("magrittr") ## need to run every time you start R and want to use %>%
library("RColorBrewer")
library("microbiome")
library("knitr")
library("vegan")
library("DESeq2")

## Load Mapping file and OTU biom file from QIIME
map <- import_qiime_sample_data("Map_Merged_RSF_Rain_Tomato.txt") 
Sys.setlocale('LC_ALL','C')
otu <- import_biom(BIOMfilename = "filtered_otu_table.biom","rep_set.tre") 

## Merge map file and OTU table and let it be a phyloseq object 
run <- merge_phyloseq(otu,map) # Use the official taxonomic ranks names 
## Use the official taxonomic ranks names 
colnames(tax_table(run)) <- c("Kingdom", "Phylum", "Class",
                                   "Order", "Family", "Genus","Species")


######## Here, we analyze several set of samples: 1) Rain samples (Rain microbiome), 2) Rain as ######## bacterial inoculum (Rain as source of the phyllosphere microbiome), 3) tomato plants never ######## exposed to rain(Greenhouse)

####################################################################################################
################################## 1) RAIN MICROBIOME ############################################

## Subsample-Only rain samples:
run_rain <- subset_samples(run, Source%in%c("Rain"))
run_rain1 <- prune_taxa(taxa_sums(run_rain) > 0, run_rain)
summary(sample_data(run_rain1)$Source)
print(run_rain1)
###### 10 Samples ########
###### 13104 taxa ########
write.csv(otu_table(run_rain1), 'run_rain_OTU_table.csv')
sample_sums(run_rain1)

########### Relative abundance of the Rain Microbiome

## Phylum relative abundance
run_phylum_rain <- run_rain1 %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Phylum)

write.csv(run_phylum_rain, 'run_phylum_rain_abun01.csv')

n <- dim(run_phylum_rain)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_phylum_rain$Phylum <- factor(run_phylum_rain$Phylum, levels = rev(levels(run_phylum_rain$Phylum)))
run_phylum_rain$Time <- factor(run_phylum_rain$Time,levels = c("Apr15","Aug15","Mar16","Apr16","May16","Jul16","Oct16","Dec16","Jun19" ))
ggplot(run_phylum_rain,aes(x=Time,y=Abundance,fill=Phylum)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Phylum > 1%) \n") +
  xlab("Rain Colection") +
  theme(axis.text.x=element_text(size=16,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Phylum Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))


##Class relative abundance
run_class_rain <- run_rain1 %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Class)

write.csv(run_class_rain, 'run_class_rain_abun01.csv')

library(RColorBrewer)
n <- dim(run_class_rain)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_class_rain$Class <- factor(run_class_rain$Class, levels = rev(levels(run_class_rain$Class)))
run_class_rain$Time <- factor(run_class_rain$Time,levels = c("Apr15","Aug15","Mar16","Apr16","May16","Jul16","Oct16","Dec16","Jun19" ))
ggplot(run_class_rain,aes(x=Time,y=Abundance,fill=Class)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Class > 1%) \n") +
  xlab("Rain Colection") +
  theme(axis.text.x=element_text(size=16,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Class Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))

##Genus relative abundance
run_genus_rain <- run_rain1 %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.001) %>%
  arrange(Genus)

write.csv(run_genus_rain, 'run_genus_rain_abun001.csv')

library(RColorBrewer)
n <- dim(run_genus_rain)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_genus_rain$Genus <- factor(run_genus_rain$Genus, levels = rev(levels(run_genus_rain$Genus)))
run_genus_rain$Time <- factor(run_genus_rain$Time,levels = c("Apr15","Aug15","Mar16","Apr16","May16","Jul16","Oct16","Dec16","Jun19" ))
ggplot(run_genus_rain,aes(x=Time,y=Abundance,fill=Genus)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus) > 2%) \n") +
  xlab("Rain Colection") +
  theme(axis.text.x=element_text(size=14,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=11),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Genus Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))

################ Rarefaction Curves #########
## Calculate rarefaction Curves
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures) # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none')) # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}

## Sequencing depth - Rarefaction curves plots
depth = 100000
step = 1000
rarefaction_curve_data <- calculate_rarefaction_curves(run_rain1, c('Observed'), rep(seq(1,depth,by=step), each = 10))
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(run_rain)), by.x = 'Sample', by.y = 'row.names')
rarefaction_curve_data_summary_verbose$SampleID <- factor(rarefaction_curve_data_summary_verbose$SampleID,levels = map$SampleID)
n <- dim(rarefaction_curve_data_summary_verbose)[1]
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot( data = rarefaction_curve_data_summary_verbose,
        mapping = aes( x = Depth, y = Alpha_diversity_mean,
                       ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                       ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                       colour = SampleID,
                       group = SampleID)
) +
  scale_fill_manual(values=col_vector) +
  geom_line( ) +
  geom_point(size=0.5 ) +
  guides(fill=guide_legend(title="Sample")) +
  theme(legend.text = element_text(size=14)) +
  theme(axis.text.x=element_text(size=12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size=16, face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10))+
  xlab("Sequences per sample") +
  ylab("Observed OTU")+
  ggtitle("Rarefaction Curves Rain Samples")


################ The Core Rain Microbiome

## Rarefaction to calculate core microbiome
run_rain.rarefied = rarefy_even_depth(run_rain1, rngseed=1, sample.size=1*min(sample_sums(run_rain1)), replace=F)
sample_sums(run_rain1)
sample_sums(run_rain.rarefied)
  ## RESULTS: 
    ## Samples were rarefied to 7133 reads per sample
    ## 6325 taxa OTU were identified


## Check the data 
print(run_rain.rarefied) 
## keep only taxa with positive sums
rain_rare.1 <- prune_taxa(taxa_sums(run_rain.rarefied) > 0, run_rain.rarefied)

## Relative abundances
rain_rare.rel <- microbiome::transform(rain_rare.1, "compositional")

## Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(rain_rare.rel, detection = 1, sort = TRUE))

#Phyloseq object of the core microbiota:
rain.core <- core(rain_rare.rel, detection = 0, prevalence = .5)

#Retrieving the associated taxa names from the phyloseq object:
rain_core.taxa <- taxa(rain.core)
class(rain_core.taxa)
# get the taxonomy data
tax.mat <- tax_table(rain.core)
tax.df <- as.data.frame(tax.mat)

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core members based on the thresholds that were used.
rain_core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% rain_core.taxa)
knitr::kable(head(rain_core.taxa.class))

#Core heatmaps

# Core with absolute counts and vertical view:
# and minimum population prevalence (given as percentage)
#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-1), log10(max(abundances(rain_rare.1))/10), length = 10)

Raincore <- plot_core(rain_rare.1, plot.type = "heatmap", 
                      prevalences = prevalences,
                      detections = detections,
                      colours = rev(brewer.pal(5, "Spectral")),
                      min.prevalence = .9)

# Data used for plotting 
df <- Raincore$data 
# list of OTUs
list <- df$Taxa 
#check the OTU ids
print(list) 
# Taxonomy data
tax <- tax_table(rain_rare.1)
tax <- as.data.frame(tax)
# Add the OTus to last column
tax$OTU <- rownames(tax)
# select taxonomy of only 
# those OTUs that are used in the plot
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 

# Merge all the column into one except the Domain as all is bacteria in this case
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Family","Genus","OTU"), sep = "_;", remove = TRUE)
### RESULT: 23 OTUs in the Core Rain microbiome
tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
# add this new information into the plot data df
df$Taxa <- tax.unit$Taxa_level
# Taxonomic information
knitr::kable(head(df))
# replace the data in the plot object
Raincore$data <- df
## Detection Threshold is the Relative Abundance in %
plot(Raincore + theme(axis.text.x=element_text(size=12,angle=0,hjust =1),
                      axis.text.y = element_text(size = 12),
                      strip.text.x = element_text(size=10,colour = "black"), 
                      strip.text.y = element_text(size=12, face = 'bold'),
                      plot.title = element_text(size = rel(2)),
                      axis.title=element_text(size=18,face="bold", vjust = 10),
                      legend.text = element_text(size=10),
                      plot.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                xlab("Detection threshold (Relative abundance %) "))


###############################################################################################
###################### 2) RAIN AS A SOURCE OF THE PHYLLOSPHERE MICROBIOME #####################

######## Number of OTUs ########
OTU_CR_day0 <- subset_samples(run, System%in%c("CR.d0"))
OTU_CR_d0 <- prune_taxa(taxa_sums(OTU_CR_day0) > 0, OTU_CR_day0)
summary(sample_data(OTU_CR_d0)$System)
print(OTU_CR_d0) ######### 8 samples #### 7994 TAXA ########

OTU_CR_day7 <- subset_samples(run, System%in%c("CR.d7"))
OTU_CR_d7 <- prune_taxa(taxa_sums(OTU_CR_day7) > 0, OTU_CR_day7)
summary(sample_data(OTU_CR_d7)$System)
print(OTU_CR_d7) ######### 9 samples #### 10672 TAXA ########

OTU_FR_day0 <- subset_samples(run, System%in%c("FR.d0"))
OTU_FR_d0 <- prune_taxa(taxa_sums(OTU_FR_day0) > 0, OTU_FR_day0)
summary(sample_data(OTU_FR_d0)$System)
print(OTU_FR_d0) ######### 8 samples #### 6644 TAXA ########

OTU_FR_day7 <- subset_samples(run, System%in%c("FR.d7"))
OTU_FR_d7 <- prune_taxa(taxa_sums(OTU_FR_day7) > 0, OTU_FR_day7)
summary(sample_data(OTU_FR_d7)$System)
print(OTU_FR_d7) ######### 8 samples #### 8982 TAXA ########

OTU_W_day0 <- subset_samples(run, System%in%c("W.d0"))
OTU_W_d0 <- prune_taxa(taxa_sums(OTU_W_day0) > 0, OTU_W_day0)
summary(sample_data(OTU_W_d0)$System)
print(OTU_W_d0) ######### 6 samples #### 5827 TAXA ########

OTU_W_day7 <- subset_samples(run, System%in%c("W.d7"))
OTU_W_d7 <- prune_taxa(taxa_sums(OTU_W_day7) > 0, OTU_W_day7)
summary(sample_data(OTU_W_d7)$System)
print(OTU_W_d7) ######### 6 samples #### 7882 TAXA ########

## Subsample-Only laboratory innoculated tomato plants with rain samples:
run_rain_Tom <- subset_samples(run, Source%in%c("Concentrated.Rain","Filtered.Rain","Sterile.Water"))
run_rain_Tom1 <- prune_taxa(taxa_sums(run_rain_Tom) > 0, run_rain_Tom)
summary(sample_data(run_rain_Tom1)$Source)
print(run_rain_Tom1)
###### 13104 taxa ########
write.csv(otu_table(run_rain_Tom1), 'run_Rain_Tomato_OTU_table.csv')
sample_sums(run_rain_Tom1)

################### Relative abundance 
## Phylum relative abundance
run_phylum_rain_Tomato <- run_rain_Tom1 %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Phylum)

write.csv(run_phylum_rain_Tomato, 'run_phylum_rain_Tomato_abun01.csv')

library(RColorBrewer)
n <- dim(run_phylum_rain_Tomato)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_phylum_rain_Tomato$Phylum <- factor(run_phylum_rain_Tomato$Phylum, levels = rev(levels(run_phylum_rain_Tomato$Phylum)))
run_phylum_rain_Tomato$Description <- factor(run_phylum_rain_Tomato$Description,levels = c("Apr15.CR.d0","Aug15.CR.d0","Mar16.CR.d0","Apr16.CR.d0","May16.CR.d0","Jul16.CR.d0","Apr15.CR.d7","Aug15.CR.d7","Mar16.CR.d7","Apr16.CR.d7","May16.CR.d7","Jul16.CR.d7","Apr15.FR.d0","Aug15.FR.d0","Mar16.FR.d0","Apr16.FR.d0","May16.FR.d0","Jul16.FR.d0","Apr15.FR.d7","Aug15.FR.d7","Mar16.FR.d7","Apr16.FR.d7","May16.FR.d7","Jul16.FR.d7","Mar16.W.d0","Apr16.W.d0","May16.W.d0","Jul16.W.d0","Mar16.W.d7","Apr16.W.d7","May16.W.d7","Jul16.W.d7"))
ggplot(run_phylum_rain_Tomato,aes(x=Description,y=Abundance,fill=Phylum)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Phylum > 1%) \n") +
  xlab("Rain Treatment") +
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Phylum Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))


##Class relative abundance
run_class_rain_Tomato <- run_rain_Tom1 %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Class)

write.csv(run_class_rain_Tomato, 'run_class_rain_Tomato_abun01.csv')

library(RColorBrewer)
n <- dim(run_class_rain_Tomato)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_class_rain_Tomato$Class <- factor(run_class_rain_Tomato$Class, levels = rev(levels(run_class_rain_Tomato$Class)))
run_class_rain_Tomato$Description <- factor(run_class_rain_Tomato$Description,levels = c("Apr15.CR.d0","Aug15.CR.d0","Mar16.CR.d0","Apr16.CR.d0","May16.CR.d0","Jul16.CR.d0","Apr15.CR.d7","Aug15.CR.d7","Mar16.CR.d7","Apr16.CR.d7","May16.CR.d7","Jul16.CR.d7","Apr15.FR.d0","Aug15.FR.d0","Mar16.FR.d0","Apr16.FR.d0","May16.FR.d0","Jul16.FR.d0","Apr15.FR.d7","Aug15.FR.d7","Mar16.FR.d7","Apr16.FR.d7","May16.FR.d7","Jul16.FR.d7","Mar16.W.d0","Apr16.W.d0","May16.W.d0","Jul16.W.d0","Mar16.W.d7","Apr16.W.d7","May16.W.d7","Jul16.W.d7"))
ggplot(run_class_rain_Tomato,aes(x=Description,y=Abundance,fill=Class)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Class > 1%) \n") +
  xlab("Rain Treatment") +
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Class Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))


## Family relative abundance
run_family_rain_Tomato <- run_rain_Tom1 %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.03) %>%
  arrange(Family)

write.csv(run_family_rain_Tomato, 'run_family_rain_Tomato_abun01.csv')

library(RColorBrewer)
n <- dim(run_family_rain_Tomato)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_family_rain_Tomato$Family <- factor(run_family_rain_Tomato$Family, levels = rev(levels(run_family_rain_Tomato$Family)))
run_family_rain_Tomato$Description <- factor(run_family_rain_Tomato$Description,levels = c("Apr15.CR.d0","Aug15.CR.d0","Mar16.CR.d0","Apr16.CR.d0","May16.CR.d0","Jul16.CR.d0","Apr15.CR.d7","Aug15.CR.d7","Mar16.CR.d7","Apr16.CR.d7","May16.CR.d7","Jul16.CR.d7","Apr15.FR.d0","Aug15.FR.d0","Mar16.FR.d0","Apr16.FR.d0","May16.FR.d0","Jul16.FR.d0","Apr15.FR.d7","Aug15.FR.d7","Mar16.FR.d7","Apr16.FR.d7","May16.FR.d7","Jul16.FR.d7","Mar16.W.d0","Apr16.W.d0","May16.W.d0","Jul16.W.d0","Mar16.W.d7","Apr16.W.d7","May16.W.d7","Jul16.W.d7"))
ggplot(run_family_rain_Tomato,aes(x=Description,y=Abundance,fill=Family)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family > 3%) \n") +
  xlab("Rain Treatment") +
  theme(axis.text.x=element_text(size=12,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Family Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))


## Genus relative abundance
run_genus_rain_Tomato <- run_rain_Tom1 %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.001) %>%
  arrange(Genus)

write.csv(run_genus_rain_Tomato, 'run_genus_rain_Tomato_abun001.csv')

library(RColorBrewer)
n <- dim(run_genus_rain_Tomato)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_genus_rain_Tomato$Genus <- factor(run_genus_rain_Tomato$Genus, levels = rev(levels(run_genus_rain_Tomato$Genus)))
run_genus_rain_Tomato$Description <- factor(run_genus_rain_Tomato$Description,levels = c("Apr15.CR.d0","Aug15.CR.d0","Mar16.CR.d0","Apr16.CR.d0","May16.CR.d0","Jul16.CR.d0","Apr15.CR.d7","Aug15.CR.d7","Mar16.CR.d7","Apr16.CR.d7","May16.CR.d7","Jul16.CR.d7","Apr15.FR.d0","Aug15.FR.d0","Mar16.FR.d0","Apr16.FR.d0","May16.FR.d0","Jul16.FR.d0","Apr15.FR.d7","Aug15.FR.d7","Mar16.FR.d7","Apr16.FR.d7","May16.FR.d7","Jul16.FR.d7","Mar16.W.d0","Apr16.W.d0","May16.W.d0","Jul16.W.d0","Mar16.W.d7","Apr16.W.d7","May16.W.d7","Jul16.W.d7"))
ggplot(run_genus_rain_Tomato,aes(x=Description,y=Abundance,fill=Genus)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus > 3%) \n") +
  xlab("Rain Treatment") +
  theme(axis.text.x=element_text(size=12,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Genus Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))

########################## MICROBIAL DIVERSITY ANALYSIS

## Subsample: laboratory innoculated tomato plants and rain samples:
run_Rain_Tomato <- subset_samples(run, Source%in%c("Rain","Concentrated.Rain","Filtered.Rain","Sterile.Water"))
summary(sample_data(run_Rain_Tomato)$Source)

######### ALPHA DIVERSITY ANALYSIS
## Calculate rarefaction Curves
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures) # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none')) # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}

## Sequencing depth - Rarefaction curves plots
depth = 100000
step = 1000
rarefaction_curve_data <- calculate_rarefaction_curves(run_Rain_Tomato, c('Observed'), rep(seq(1,depth,by=step), each = 10))
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(run_Rain_Tomato)), by.x = 'Sample', by.y = 'row.names')
rarefaction_curve_data_summary_verbose$SampleID <- factor(rarefaction_curve_data_summary_verbose$SampleID,levels = map$SampleID)
n <- dim(rarefaction_curve_data_summary_verbose)[1]
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot( data = rarefaction_curve_data_summary_verbose,
        mapping = aes( x = Depth, y = Alpha_diversity_mean,
                       ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                       ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                       colour = System,
                       group = SampleID)
) +
  scale_fill_manual(values=col_vector) +
  geom_line( ) +
  geom_point(size=0.5 ) +
  guides(fill=guide_legend(title="Sample")) +
  theme(legend.text = element_text(size=14)) +
  theme(axis.text.x=element_text(size=12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size=16, face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10))+
  xlab("Sequences per sample") +
  ylab("Observed OTU")+
  ggtitle("Rarefaction Curves Rain-Tomato")

############# ALPHA DIVERSITY INDICES 
##### Rarefaction  
run.rare_Rain_Tomato = rarefy_even_depth(run_Rain_Tomato, rngseed=1, sample.size=1*min(sample_sums(run_Rain_Tomato)), replace=F)

sample_sums(run_Rain_Tomato)
sample_sums(run.rare_Rain_Tomato)
### RESULTS: Samples rarefied to 7133 reads

############ PLOT INDICES
plot_richness(run.rare_Rain_Tomato,x="System",measures=c("Observed","Shannon","Simpson")) + 
  geom_boxplot() +
  ylab("Alpha Diversity") + 
  theme(axis.text.x=element_text(angle=90,hjust=0.5,size=10),
        axis.title = element_text(size=18),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=18,face="bold"))

rich = estimate_richness(run.rare_Rain_Tomato)
write.csv(rich, 'richness_Rain_tomato_experiment.csv')

pairwise.wilcox.test(rich$Observed, sample_data(run.rare_Rain_Tomato)$System)
### RESULTS: CR.d0 vs CR.d7 = 0.020,   FR.d0 vs FR.d7 = 0.135,   W.d0 vs W.d7 = 0.312
pairwise.wilcox.test(rich$Shannon, sample_data(run.rare_Rain_Tomato)$System)
### RESULTS: CR.d0 vs CR.d7 = 0.049,   FR.d0 vs FR.d7 = 0.023,   W.d0 vs W.d7 = 0.281
pairwise.wilcox.test(rich$Simpson, sample_data(run.rare_Rain_Tomato)$System)
### RESULTS: CR.d0 vs CR.d7 = 0.58,   FR.d0 vs FR.d7 = 1.0,   W.d0 vs W.d7 = 0.91

############### BETA DIVERSITY ANALYSIS

#### Pairwise distances between samples are calculated: Weighted-Unifrac and unweighted-Unifrac.
dm_weighted_unifrac <- phyloseq::distance(run.rare_Rain_Tomato, method = "wUniFrac")
#### Plot
ordWU_rain_Tom <- ordinate(run.rare_Rain_Tomato,method='PCoA',distance=dm_weighted_unifrac) 
plot_ordination(run.rare_Rain_Tomato, ordWU_rain_Tom,,color = 'Source', shape='DayPoint') + 
  geom_point(size=3) + 
  theme(axis.text.x=element_text(size=14,angle=0,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ggtitle("PCoA: Weighted Unifrac") +
  scale_y_continuous(labels=percent_format(),expand=c(0,0))

adonis(dm_weighted_unifrac ~ sample_data(run.rare_Rain_Tomato)$System)
### RESULTS: Permutations: 999, Df:6, R2:1.9671 Pr(>F): 0.017

dm_unweighted_unifrac <- phyloseq::distance(run.rare_Rain_Tomato, method='Unifrac')
#### Plot
ordUU_rain_Tom <- ordinate(run.rare_Rain_Tomato,method='PCoA',distance=dm_unweighted_unifrac) 
plot_ordination(run.rare_Rain_Tomato, ordUU_rain_Tom,,color = 'Source',shape='DayPoint') + 
  geom_point(size=3) + 
  theme(axis.text.x=element_text(size=14,angle=0,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ggtitle("PCoA: Unweighted Unifrac") +
  scale_y_continuous(labels=percent_format(),expand=c(0,0))
  
adonis(dm_unweighted_unifrac ~ sample_data(run.rare_Rain_Tomato)$System)
### RESULTS: Permutations: 999, Df:6, F.Model: 1.6631; R2:0.17214; Pr(>F): 0.001

bray_diss_rain_Tom = phyloseq::distance(run.rare_Rain_Tomato, method="bray")
#### Plot
ordination_rain_Tom = ordinate(run.rare_Rain_Tomato, method="PCoA", distance=bray_diss_rain_Tom)
plot_ordination(run.rare_Rain_Tomato, ordination_rain_Tom, color = 'Source',shape='DayPoint') + theme(aspect.ratio=1)+
  geom_point(size=3) + 
  theme(axis.text.x=element_text(size=14,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ggtitle("PCoA: Bray-Curtis")

adonis(bray_diss_rain_Tom ~ sample_data(run.rare_Rain_Tomato)$System)
### RESULTS: Permutations: 999, Df:6, F.Model: 2.2065; R2:0.21618; Pr(>F): 0.001


############## DIFERENCTIAL ABUNDANCES RAIN AS SOURCE OF PHYLLOSPHERE MICROBIOME
### Tomato plants treated with Rain: day0 vs day 7
run_RAIN <- subset_samples(run, System%in%c("CR.d0","CR.d7"))
run_RAIN1 <- prune_taxa(taxa_sums(run_RAIN) > 0, run_RAIN)
summary(sample_data(run_RAIN1)$System)
print(run_RAIN1)
###### 17 samples
###### 13904 taxa
sample_sums(run_RAIN1)

run_RAIN1 <- prune_samples(sample_sums(run_RAIN1) > 500, run_RAIN1)
head(sample_data(run_RAIN1)$System, 25)

deseq_rain = phyloseq_to_deseq2(run_RAIN1, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_rain), 1, gm_mean)
deseq_rain = estimateSizeFactors(deseq_rain, geoMeans = geoMeans)
deseq_rain = DESeq(deseq_rain, fitType="local")

res = results(deseq_rain)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(run_RAIN1)[rownames(sigtab), ], "matrix"))
head(sigtab)
##### To write all OTUs that were significant different: positives and negatives
sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab, 'DEseq_all_values_RAINd0_vs_RAINd7.csv')
###### RESULT: 116 TAXA differentially abundant ##############
##### To subset positives values
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab, 'Differential_abundance_RAINd0_vs_RAINd7.csv')
###### RESULT: 104 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Tomato treated with CR d0 vs d7")


### Tomato plants treated with Filtered Rain: day0 vs day 7
run_FRain <- subset_samples(run, System%in%c("FR.d0","FR.d7"))
run_FRAIN1 <- prune_taxa(taxa_sums(run_FRain) > 0, run_FRain)
summary(sample_data(run_FRAIN1)$System)
print(run_FRAIN1)
###### 16 samples
###### 11298 taxa
sample_sums(run_FRAIN1)

run_FRAIN1 <- prune_samples(sample_sums(run_FRAIN1) > 500, run_FRAIN1)
head(sample_data(run_FRAIN1)$System, 25)
deseq_FR = phyloseq_to_deseq2(run_FRAIN1, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_FR), 1, gm_mean)
deseq_FR = estimateSizeFactors(deseq_FR, geoMeans = geoMeans)
deseq_FR = DESeq(deseq_FR, fitType="local")

res_FR = results(deseq_FR)
res_FR = res_FR[order(res_FR$padj, na.last=NA), ]
alpha = 0.01
sigtab_FR = res_FR[(res_FR$padj < alpha), ]
sigtab_FR = cbind(as(sigtab_FR, "data.frame"), as(tax_table(run_FRAIN1)[rownames(sigtab_FR), ], "matrix"))
head(sigtab_FR)
##### To write all OTUs that were significant different: positives and negatives
sigtab_FR = sigtab_FR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_FR, 'DEseq_all_values_FilRAINd0_vs_FilRAINd7.csv')
###### RESULT: 0 TAXA differentially abundant ##############
##### To subset positives values
posigtab_FR = sigtab_FR[sigtab_FR[, "log2FoldChange"] > 0, ]
posigtab_FR = posigtab_FR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab_FR, 'Differential_abundance_FR_d0_vs_FR_d7.csv')
############ RESULT: 0 TAXA that are differentially abundant ###########


### Tomato plants treated with Sterile Water: day0 vs day 7
run_SWater <- subset_samples(run, System%in%c("W.d0","W.d7"))
run_SW1 <- prune_taxa(taxa_sums(run_SWater) > 0, run_SWater)
summary(sample_data(run_SW1)$System)
print(run_SW1)
###### 12 samples
###### 9929 taxa
sample_sums(run_SW1)
run_SW1 <- prune_samples(sample_sums(run_SW1) > 500, run_SW1)
head(sample_data(run_SW1)$System, 25)
deseq_FR = phyloseq_to_deseq2(run_SW1, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_W), 1, gm_mean)
deseq_W = estimateSizeFactors(deseq_W, geoMeans = geoMeans)
deseq_W = DESeq(deseq_W, fitType="local")

res_W = results(deseq_W)
res_W = res_W[order(res_W$padj, na.last=NA), ]
alpha = 0.01
sigtab_W = res_W[(res_W$padj < alpha), ]
sigtab_W = cbind(as(sigtab_W, "data.frame"), as(tax_table(run_SW1)[rownames(sigtab_W), ], "matrix"))
head(sigtab_W)
##### To write all OTUs that were significant different: positives and negatives
sigtab_W = sigtab_W[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_W, 'DEseq_all_values_SWaterd0_vs_SWaterd7.csv')
###### RESULT: 0 TAXA differentially abundant ##############
##### To subset positives values
posigtab_W = sigtab_W[sigtab_W[, "log2FoldChange"] > 0, ]
posigtab_W = posigtab_W[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab_W, 'Differential_abundance_water_0_7.csv')
########## 12 TAXA are differentially abundant ###########

library("ggplot2")
theme_set(theme_bw())
sigtabgen_W = subset(sigtab_W, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_W$log2FoldChange, sigtabgen_W$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_W$Phylum = factor(as.character(sigtabgen_W$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_W$log2FoldChange, sigtabgen_W$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_W$Genus = factor(as.character(sigtabgen_W$Genus), levels=names(x))
ggplot(sigtabgen_W, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Tomato treated with Water d0 vs d7")


### Rain vs Tomato treated with C.Rain.day7
run_R_R7 <- subset_samples(run, System%in%c("CR.d7","Atm.Rain"))
run_R_CR7 <- prune_taxa(taxa_sums(run_R_R7) > 0, run_R_R7)
summary(sample_data(run_R_CR7)$System)
print(run_R_CR7)
###### 19 samples
###### 18870 taxa
sample_sums(run_R_CR7)

run_R_CR7 <- prune_samples(sample_sums(run_R_CR7) > 500, run_R_CR7)
head(sample_data(run_R_CR7)$System, 25)
deseq_Rain_7 = phyloseq_to_deseq2(run_R_CR7, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Rain_7), 1, gm_mean)
deseq_Rain_7 = estimateSizeFactors(deseq_Rain_7, geoMeans = geoMeans)
deseq_Rain_7 = DESeq(deseq_Rain_7, fitType="local")

res_Rain_7 = results(deseq_Rain_7)
res_Rain_7 = res_Rain_7[order(res_Rain_7$padj, na.last=NA), ]
alpha = 0.01
sigtab_Rain_7 = res_Rain_7[(res_Rain_7$padj < alpha), ]
sigtab_Rain_7 = cbind(as(sigtab_Rain_7, "data.frame"), as(tax_table(run_R_CR7)[rownames(sigtab_Rain_7), ], "matrix"))
head(sigtab_Rain_7)
##### To write all OTUs that were significant different: positives and negatives
sigtab_Rain_7 = sigtab_Rain_7[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_Rain_7, 'DEseq_all_values_Rain_vs_CRd7.csv')
###### RESULT: 187 TAXA differentially abundant ##############
##### To subset positives values
posigtab_Rain_7 = sigtab_Rain_7[sigtab_Rain_7[, "log2FoldChange"] > 0, ]
posigtab_Rain_7 = posigtab_Rain_7[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab_Rain_7, 'Differential_abundance_Rain_vs_rain7.csv')
###### RESULT: 66 TAXA differentially abundant ##############
theme_set(theme_bw())
sigtabgen_Rain_7 = subset(sigtab_Rain_7, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_Rain_7$log2FoldChange, sigtabgen_Rain_7$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_7$Phylum = factor(as.character(sigtabgen_Rain_7$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_Rain_7$log2FoldChange, sigtabgen_Rain_7$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_7$Genus = factor(as.character(sigtabgen_Rain_7$Genus), levels=names(x))
ggplot(sigtabgen_Rain_7, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Rain vs Tomato treated with CR.d7")


### Rain vs FilteredRain day 7
run_Rain_FR7 <- subset_samples(run, System%in%c("FR.d7","Atm.Rain"))
run_R_FR7 <- prune_taxa(taxa_sums(run_Rain_FR7) > 0, run_Rain_FR7)
summary(sample_data(run_R_FR7)$System)
print(run_R_FR7)
###### 18 samples
###### 18145 taxa
sample_sums(run_R_FR7)

run_R_FR7 <- prune_samples(sample_sums(run_R_FR7) > 500, run_R_FR7)
head(sample_data(run_R_FR7)$System, 25)
deseq_Rain_FR7 = phyloseq_to_deseq2(run_R_FR7, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Rain_FR7), 1, gm_mean)
deseq_Rain_FR7 = estimateSizeFactors(deseq_Rain_FR7, geoMeans = geoMeans)
deseq_Rain_FR7 = DESeq(deseq_Rain_FR7, fitType="local")

res_Rain_FR7 = results(deseq_Rain_FR7)
res_Rain_FR7 = res_Rain_FR7[order(res_Rain_FR7$padj, na.last=NA), ]
alpha = 0.01
sigtab_Rain_FR7 = res_Rain_FR7[(res_Rain_FR7$padj < alpha), ]
sigtab_Rain_FR7 = cbind(as(sigtab_Rain_FR7, "data.frame"), as(tax_table(run_R_FR7)[rownames(sigtab_Rain_FR7), ], "matrix"))
head(sigtab_Rain_FR7)
##### To write all OTUs that were significant different: positives and negatives
sigtab_Rain_FR7 = sigtab_Rain_FR7[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_Rain_FR7, 'DEseq_all_values_Rain_vs_FRd7.csv')
###### RESULT: 294 TAXA differentially abundant ##############
##### To subset positives values
posigtab_Rain_FR7 = sigtab_Rain_FR7[sigtab_Rain_FR7[, "log2FoldChange"] > 0, ]
posigtab_Rain_FR7 = posigtab_Rain_FR7[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_Rain_FR7, 'Differential_abundance_Rain_vs_FR7.csv')
###### RESULT: 91 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_Rain_FR7 = subset(sigtab_Rain_FR7, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_Rain_FR7$log2FoldChange, sigtabgen_Rain_FR7$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_FR7$Phylum = factor(as.character(sigtabgen_Rain_FR7$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_Rain_FR7$log2FoldChange, sigtabgen_Rain_FR7$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_FR7$Genus = factor(as.character(sigtabgen_Rain_FR7$Genus), levels=names(x))
ggplot(sigtabgen_Rain_FR7, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 9),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Rain vs tomato treated with FR.d7")


### Rain vs Tomato treated with Sterile Water day 7
run_R_SWater7 <- subset_samples(run, System%in%c("W.d7","Atm.Rain"))
run_R_SW7 <- prune_taxa(taxa_sums(run_R_SWater7) > 0, run_R_SWater7)
summary(sample_data(run_R_SW7)$System)
print(run_R_SW7)
###### 16 samples
###### 17829 taxa
sample_sums(run_R_SW7)

run_R_SW7 <- prune_samples(sample_sums(run_R_SW7) > 500, run_R_SW7)
head(sample_data(run_R_SW7)$System, 25)
deseq_Rain_W7 = phyloseq_to_deseq2(run_R_SW7, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Rain_W7), 1, gm_mean)
deseq_Rain_W7 = estimateSizeFactors(deseq_Rain_W7, geoMeans = geoMeans)
deseq_Rain_W7 = DESeq(deseq_Rain_W7, fitType="local")

res_Rain_W7 = results(deseq_Rain_W7)
res_Rain_W7 = res_Rain_W7[order(res_Rain_W7$padj, na.last=NA), ]
alpha = 0.01
sigtab_Rain_W7 = res_Rain_W7[(res_Rain_W7$padj < alpha), ]
sigtab_Rain_W7 = cbind(as(sigtab_Rain_W7, "data.frame"), as(tax_table(run_R_SW7)[rownames(sigtab_Rain_W7), ], "matrix"))
head(sigtab_Rain_W7)
##### To write all OTUs that were significant different: positives and negatives
sigtab_Rain_W7 = sigtab_Rain_W7[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_Rain_W7, 'DEseq_all_values_Rain_vs_SWd7.csv')
###### RESULT: 364 TAXA differentially abundant ##############
##### To subset positives values
posigtab_Rain_W7 = sigtab_Rain_W7[sigtab_Rain_W7[, "log2FoldChange"] > 0, ]
posigtab_Rain_W7 = posigtab_Rain_W7[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_Rain_W7, 'Differential_abundance_Rain_Water7.csv')
###### RESULT: 154 TAXA differentially abundant ##############
theme_set(theme_bw())
sigtabgen_Rain_W7 = subset(sigtab_Rain_W7, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_Rain_W7$log2FoldChange, sigtabgen_Rain_W7$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_W7$Phylum = factor(as.character(sigtabgen_Rain_W7$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_Rain_W7$log2FoldChange, sigtabgen_Rain_W7$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_W7$Genus = factor(as.character(sigtabgen_Rain_W7$Genus), levels=names(x))
ggplot(sigtabgen_Rain_W7, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 9),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Rain vs tomato treated with W.d7")


###### Rain vs Tomato Treated with rain day 0
run_R_R0 <- subset_samples(run, System%in%c("CR.d0","Atm.Rain"))
run_R_CR0 <- prune_taxa(taxa_sums(run_R_R0) > 0, run_R_R0)
summary(sample_data(run_R_CR0)$System)
print(run_R_CR0)
###### 18 samples
###### 16917 taxa
sample_sums(run_R_CR0)

run_R_CR0 <- prune_samples(sample_sums(run_R_CR0) > 500, run_R_CR0)
head(sample_data(run_R_CR0)$System, 25)
deseq_Rain_0 = phyloseq_to_deseq2(run_R_CR0, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Rain_0), 1, gm_mean)
deseq_Rain_0 = estimateSizeFactors(deseq_Rain_0, geoMeans = geoMeans)
deseq_Rain_0 = DESeq(deseq_Rain_0, fitType="local")

res_Rain_0 = results(deseq_Rain_0)
res_Rain_0 = res_Rain_0[order(res_Rain_0$padj, na.last=NA), ]
alpha = 0.01
sigtab_Rain_0 = res_Rain_0[(res_Rain_0$padj < alpha), ]
sigtab_Rain_0 = cbind(as(sigtab_Rain_0, "data.frame"), as(tax_table(run_R_CR0)[rownames(sigtab_Rain_0), ], "matrix"))
head(sigtab_Rain_0)
##### To write all OTUs that were significant different: positives and negatives
sigtab_Rain_0 = sigtab_Rain_0[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_Rain_0, 'DEseq_all_values_Rain_vs_CRd0.csv')
###### RESULT: 292 TAXA differentially abundant ##############
##### To subset positives values
posigtab_Rain_0 = sigtab_Rain_0[sigtab_Rain_0[, "log2FoldChange"] > 0, ]
posigtab_Rain_0 = posigtab_Rain_0[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_Rain_0, 'Differential_abundance_Rain_vs_rain_tomato0.csv')
###### RESULT: 148 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_Rain_0 = subset(sigtab_Rain_0, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_Rain_0$log2FoldChange, sigtabgen_Rain_0$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_0$Phylum = factor(as.character(sigtabgen_Rain_0$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_Rain_0$log2FoldChange, sigtabgen_Rain_0$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_0$Genus = factor(as.character(sigtabgen_Rain_0$Genus), levels=names(x))
ggplot(sigtabgen_Rain_0, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 9),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Rain vs Tomato treated with CR.d0")


### Rain vs tomato treated with Filtered Rain day 0
run_Rain_FR0 <- subset_samples(run, System%in%c("FR.d0","Atm.Rain"))
run_R_FR0 <- prune_taxa(taxa_sums(run_Rain_FR0) > 0, run_Rain_FR0)
summary(sample_data(run_R_FR0)$System)
print(run_R_FR0)
###### 18 samples
###### 16784 taxa
sample_sums(run_R_FR0)

run_R_FR0 <- prune_samples(sample_sums(run_R_FR0) > 500, run_R_FR0)
head(sample_data(run_R_FR0)$System, 25)
deseq_Rain_FR0 = phyloseq_to_deseq2(run_R_FR0, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Rain_FR0), 1, gm_mean)
deseq_Rain_FR0 = estimateSizeFactors(deseq_Rain_FR0, geoMeans = geoMeans)
deseq_Rain_FR0 = DESeq(deseq_Rain_FR0, fitType="local")

res_Rain_FR0 = results(deseq_Rain_FR0)
res_Rain_FR0 = res_Rain_FR0[order(res_Rain_FR0$padj, na.last=NA), ]
alpha = 0.01
sigtab_Rain_FR0 = res_Rain_FR0[(res_Rain_FR0$padj < alpha), ]
sigtab_Rain_FR0 = cbind(as(sigtab_Rain_FR0, "data.frame"), as(tax_table(run_R_FR0)[rownames(sigtab_Rain_FR0), ], "matrix"))
head(sigtab_Rain_FR0)
##### To write all OTUs that were significant different: positives and negatives
sigtab_Rain_FR0 = sigtab_Rain_FR0[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_Rain_FR0, 'DEseq_all_values_Rain_vs_FRd0.csv')
###### RESULT: 472 TAXA differentially abundant ##############
##### To subset positives values
posigtab_Rain_FR0 = sigtab_Rain_FR0[sigtab_Rain_FR0[, "log2FoldChange"] > 0, ]
posigtab_Rain_FR0 = posigtab_Rain_FR0[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_Rain_FR0, 'Differential_abundance_Rain_vs_tomato_FR0.csv')
###### RESULT: 178 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_Rain_FR0 = subset(sigtab_Rain_FR0, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_Rain_FR0$log2FoldChange, sigtabgen_Rain_FR0$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_FR0$Phylum = factor(as.character(sigtabgen_Rain_FR0$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_Rain_FR0$log2FoldChange, sigtabgen_Rain_FR0$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_FR0$Genus = factor(as.character(sigtabgen_Rain_FR0$Genus), levels=names(x))
ggplot(sigtabgen_Rain_FR0, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Rain vs Tomato treated with FR.d0")


### Rain vs Tomato treated with sterile Water day 0
run_R_SWater0 <- subset_samples(run, System%in%c("W.d0","Atm.Rain"))
run_R_SW0 <- prune_taxa(taxa_sums(run_R_SWater0) > 0, run_R_SWater0)
summary(sample_data(run_R_SW0)$System)
print(run_R_SW0)
###### 16 samples
###### 16333 taxa
sample_sums(run_R_SW0)

run_R_SW0 <- prune_samples(sample_sums(run_R_SW0) > 500, run_R_SW0)
head(sample_data(run_R_SW0)$System, 25)
deseq_Rain_W0 = phyloseq_to_deseq2(run_R_SW0, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Rain_W0), 1, gm_mean)
deseq_Rain_W0 = estimateSizeFactors(deseq_Rain_W0, geoMeans = geoMeans)
deseq_Rain_W0 = DESeq(deseq_Rain_W0, fitType="local")

res_Rain_W0 = results(deseq_Rain_W0)
res_Rain_W0 = res_Rain_W0[order(res_Rain_W0$padj, na.last=NA), ]
alpha = 0.01
sigtab_Rain_W0 = res_Rain_W0[(res_Rain_W0$padj < alpha), ]
sigtab_Rain_W0 = cbind(as(sigtab_Rain_W0, "data.frame"), as(tax_table(run_R_SW0)[rownames(sigtab_Rain_W0), ], "matrix"))
head(sigtab_Rain_W0)
##### To write all OTUs that were significant different: positives and negatives
sigtab_Rain_W0 = sigtab_Rain_W0[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(sigtab_Rain_W0, 'DEseq_all_values_Rain_vs_SWd0.csv')
###### RESULT: 212 TAXA differentially abundant ##############
##### To subset positives values
posigtab_Rain_W0 = sigtab_Rain_W0[sigtab_Rain_W0[, "log2FoldChange"] > 0, ]
posigtab_Rain_W0 = posigtab_Rain_W0[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
###### RESULT: 77 TAXA differentially abundant ##############
write.csv(posigtab_Rain_W0, 'Differential_abundance_Rain_vs_tomato_Water0.csv')

theme_set(theme_bw())
sigtabgen_Rain_W0 = subset(sigtab_Rain_W0, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_Rain_W0$log2FoldChange, sigtabgen_Rain_W0$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_W0$Phylum = factor(as.character(sigtabgen_Rain_W0$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_Rain_W0$log2FoldChange, sigtabgen_Rain_W0$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_Rain_W0$Genus = factor(as.character(sigtabgen_Rain_W0$Genus), levels=names(x))
ggplot(sigtabgen_Rain_W0, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=16,angle=0,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("Rain vs Tomato treated with W.d0")


#################################################################################################
###################### 3) TOMATO LEAF MICROBIOME GROWN ON GREENHOUSE SYSTEM #####################

##### Number of OTUs #######
RSF_Hydroponics <- subset_samples(run, System%in%c("Hydroponic"))
RSF_H <- prune_taxa(taxa_sums(RSF_Hydroponics) > 0, RSF_Hydroponics)
summary(sample_data(RSF_H)$System)
print(RSF_H)  ###### 29 samples #### 13910 TAXA ######

RSF_Organics <- subset_samples(run, System%in%c("Organic"))
RSF_O <- prune_taxa(taxa_sums(RSF_Organics) > 0, RSF_Organics)
summary(sample_data(RSF_O)$System)
print(RSF_O)  ###### 18 Samples ####10157 TAXA ######

RSF_OutRain <- subset_samples(run, System%in%c("Out.Rain"))
RSF_OUTR <- prune_taxa(taxa_sums(RSF_OutRain) > 0, RSF_OutRain)
summary(sample_data(RSF_OUTR)$System)
print(RSF_OUTR)  ###### 7 samples #### 10190 TAXA ######

RSF_OutNoRain <- subset_samples(run, System%in%c("Out.No.Rain"))
RSF_OUTNR <- prune_taxa(taxa_sums(RSF_OutNoRain) > 0, RSF_OutNoRain)
summary(sample_data(RSF_OUTNR)$System)
print(RSF_OUTNR)  ###### 3 samples #### 7691 TAXA ######

## Subsample-Only Greenhouse and plants grown outside samples:
run_RSF_TOM <- subset_samples(run, System%in%c("Hydroponic","Organic","Out.No.Rain","Out.Rain"))
run_RSF_Tomato <- prune_taxa(taxa_sums(run_RSF_TOM) > 0, run_RSF_TOM)
summary(sample_data(run_RSF_Tomato)$System)
print(run_RSF_Tomato)
###### 57 samples
###### 23146 taxa
sample_sums(run_RSF_Tomato)

################ Relative abundance 
#### Phylum Relative abundance
run_RSF_phylum <- run_RSF_Tomato %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Phylum)

n <- dim(run_RSF_phylum)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_RSF_phylum$Phylum <- factor(run_RSF_phylum$Phylum, levels = rev(levels(run_RSF_phylum$Phylum)))
run_RSF_phylum$SampleID <- factor(run_RSF_phylum$SampleID, levels = c("RSF1.H.T1","RSF1.H.T2","RSF1.H.T3","RSF1.H.B1","RSF1.H.B2","RSF1.H.B3","RSF2.H.T1","RSF2.H.T2","RSF2.H.T3","RSF3.H.1","RSF3.H.2","RSF3.H.3","RSF4.H.1","RSF4.H.2","RSF4.H.3","RSF4.H.4","RSF5.HK.1","RSF5.HK.2","RSF5.HK.3","RSF5.HQ.1","RSF5.HQ.2","RSF5.HQ.3","RSF6.HQ.1","RSF6.HQ.2","RSF6.HQ.3","RSF6.HK.1","RSF6.HK.2","RSF6.HK.3","RSF6.HK.4","RSF1.S.T1","RSF2.S.T1","RSF2.SK.T1","RSF2.S.T2","RSF2.SK.T2","RSF2.S.T3","RSF2.SK.T3","RSF3.S.1","RSF3.S.2","RSF3.S.3","RSF3.S.4","RSF4.S.1","RSF4.S.2","RSF4.S.3","RSF4.S.4","RSF5.SK.1","RSF5.SK.2","RSF5.SK.3","RioG.Roof.1","RioG.Roof.2","Kom.Roof","TR.6.17.Latham","TR.6.21.Latham","TR.7.15.AREC","TR.7.25.AREC","TNR.6.17.Latham","TNR.6.21.Latham","TNR.7.15.AREC"))

ggplot(run_RSF_phylum,aes(x=SampleID,y=Abundance,fill=Phylum)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Phylum > 1%) \n") +
  xlab("Sample") +
  theme(axis.text.x=element_text(size=11,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=10),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Phylum Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))


### Class Relative abundance
run_RSF_class <- run_RSF_Tomato %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Class)

n <- dim(run_RSF_class)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_RSF_class$Class <- factor(run_RSF_class$Class, levels = rev(levels(run_RSF_class$Class)))
run_RSF_class$SampleID <- factor(run_RSF_class$SampleID, levels = c("RSF1.H.T1","RSF1.H.T2","RSF1.H.T3","RSF1.H.B1","RSF1.H.B2","RSF1.H.B3","RSF2.H.T1","RSF2.H.T2","RSF2.H.T3","RSF3.H.1","RSF3.H.2","RSF3.H.3","RSF4.H.1","RSF4.H.2","RSF4.H.3","RSF4.H.4","RSF5.HK.1","RSF5.HK.2","RSF5.HK.3","RSF5.HQ.1","RSF5.HQ.2","RSF5.HQ.3","RSF6.HQ.1","RSF6.HQ.2","RSF6.HQ.3","RSF6.HK.1","RSF6.HK.2","RSF6.HK.3","RSF6.HK.4","RSF1.S.T1","RSF2.S.T1","RSF2.SK.T1","RSF2.S.T2","RSF2.SK.T2","RSF2.S.T3","RSF2.SK.T3","RSF3.S.1","RSF3.S.2","RSF3.S.3","RSF3.S.4","RSF4.S.1","RSF4.S.2","RSF4.S.3","RSF4.S.4","RSF5.SK.1","RSF5.SK.2","RSF5.SK.3","RioG.Roof.1","RioG.Roof.2","Kom.Roof","TR.6.17.Latham","TR.6.21.Latham","TR.7.15.AREC","TR.7.25.AREC","TNR.6.17.Latham","TNR.6.21.Latham","TNR.7.15.AREC"))

ggplot(run_RSF_class,aes(x=SampleID,y=Abundance,fill=Class)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Class > 1%) \n") +
  xlab("Sample") +
  theme(axis.text.x=element_text(size=11,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=10),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Class Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))


#### Family Relative Abundance
run_RSF_family <- run_RSF_Tomato %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.03) %>%
  arrange(Family)

n <- dim(run_RSF_family)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_RSF_family$Family <- factor(run_RSF_family$Family, levels = rev(levels(run_RSF_family$Family)))
run_RSF_family$SampleID <- factor(run_RSF_family$SampleID, levels = c("RSF1.H.T1","RSF1.H.T2","RSF1.H.T3","RSF1.H.B1","RSF1.H.B2","RSF1.H.B3","RSF2.H.T1","RSF2.H.T2","RSF2.H.T3","RSF3.H.1","RSF3.H.2","RSF3.H.3","RSF4.H.1","RSF4.H.2","RSF4.H.3","RSF4.H.4","RSF5.HK.1","RSF5.HK.2","RSF5.HK.3","RSF5.HQ.1","RSF5.HQ.2","RSF5.HQ.3","RSF6.HQ.1","RSF6.HQ.2","RSF6.HQ.3","RSF6.HK.1","RSF6.HK.2","RSF6.HK.3","RSF6.HK.4","RSF1.S.T1","RSF2.S.T1","RSF2.SK.T1","RSF2.S.T2","RSF2.SK.T2","RSF2.S.T3","RSF2.SK.T3","RSF3.S.1","RSF3.S.2","RSF3.S.3","RSF3.S.4","RSF4.S.1","RSF4.S.2","RSF4.S.3","RSF4.S.4","RSF5.SK.1","RSF5.SK.2","RSF5.SK.3","RioG.Roof.1","RioG.Roof.2","Kom.Roof","TR.6.17.Latham","TR.6.21.Latham","TR.7.15.AREC","TR.7.25.AREC","TNR.6.17.Latham","TNR.6.21.Latham","TNR.7.15.AREC"))

ggplot(run_RSF_family,aes(x=SampleID,y=Abundance,fill=Family)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family > 3%) \n") +
  xlab("Sample") +
  theme(axis.text.x=element_text(size=11,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Family Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))


##Genus relative abundance
run_RSF_genus <- run_RSF_Tomato %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.05) %>%
  arrange(Genus)

n <- dim(run_RSF_genus)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

run_RSF_genus$Genus <- factor(run_RSF_genus$Genus, levels = rev(levels(run_RSF_genus$Genus)))
run_RSF_genus$SampleID <- factor(run_RSF_genus$SampleID, levels = c("RSF1.H.T1","RSF1.H.T2","RSF1.H.T3","RSF1.H.B1","RSF1.H.B2","RSF1.H.B3","RSF2.H.T1","RSF2.H.T2","RSF2.H.T3","RSF3.H.1","RSF3.H.2","RSF3.H.3","RSF4.H.1","RSF4.H.2","RSF4.H.3","RSF4.H.4","RSF5.HK.1","RSF5.HK.2","RSF5.HK.3","RSF5.HQ.1","RSF5.HQ.2","RSF5.HQ.3","RSF6.HQ.1","RSF6.HQ.2","RSF6.HQ.3","RSF6.HK.1","RSF6.HK.2","RSF6.HK.3","RSF6.HK.4","RSF1.S.T1","RSF2.S.T1","RSF2.SK.T1","RSF2.S.T2","RSF2.SK.T2","RSF2.S.T3","RSF2.SK.T3","RSF3.S.1","RSF3.S.2","RSF3.S.3","RSF3.S.4","RSF4.S.1","RSF4.S.2","RSF4.S.3","RSF4.S.4","RSF5.SK.1","RSF5.SK.2","RSF5.SK.3","RioG.Roof.1","RioG.Roof.2","Kom.Roof","TR.6.17.Latham","TR.6.21.Latham","TR.7.15.AREC","TR.7.25.AREC","TNR.6.17.Latham","TNR.6.21.Latham","TNR.7.15.AREC"))

ggplot(run_RSF_genus,aes(x=SampleID,y=Abundance,fill=Genus)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus > 5%) \n") +
  xlab("Sample") +
  theme(axis.text.x=element_text(size=11,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=8),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Genus Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))

### Alpha diversity
#### Calculate rarefaction
#First we create a function that can iterate the BIOM file, i.e. rarefaction
### Alpha diversity
#### Calculate rarefaction
#First we create a function that can iterate the BIOM file, i.e. rarefaction

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures) # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none')) # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}

#One thing worth to note is the choice of depth. For the most of the time, rarefaction curve is used to compare the richness of each sample, so the depth is the minimal number of OTU of the samples. However, if you want to see the raw sample size, use the maximum.

depth = 150000
step = 1000
rarefaction_curve_data <- calculate_rarefaction_curves(run_RSF_Tomato, c('Observed'), rep(seq(1,depth,by=step), each = 10))
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(run)), by.x = 'Sample', by.y = 'row.names')
rarefaction_curve_data_summary_verbose$SampleID <- factor(rarefaction_curve_data_summary_verbose$SampleID,levels = map$SampleID)
n <- dim(rarefaction_curve_data_summary_verbose)[1]

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot( data = rarefaction_curve_data_summary_verbose,
        mapping = aes( x = Depth, y = Alpha_diversity_mean,
                       ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                       ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                       colour = System,
                       group = SampleID)
) +
  scale_fill_manual(values=col_vector) +
  geom_line( ) +
  geom_point(size=0.5 ) +
  guides(fill=guide_legend(title="Sample")) +
  theme(legend.text = element_text(size=14)) +
  theme(axis.text.x=element_text(size=12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size=16, face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=16,face="bold", vjust = 10))+
  xlab("Sequences per sample") +
  ylab("Observed OTU")+
  ggtitle("Rarefaction Curves RSF-Outside")

########## Rarefaction for Diversity analysis
run_RSF_rarefied = rarefy_even_depth(run_RSF_Tomato, rngseed=1, sample.size=1*min(sample_sums(run_RSF_Tomato)), replace=F)

sample_sums(run_RSF_Tomato)
print(run_RSF_Tomato)
####### 57 samples ###### 23146 taxa ########
sample_sums(run_RSF_rarefied)
print(run_RSF_rarefied)
####### 57 samples ###### 10525 taxa ########
####### Samples were rarefied to 2546 reads/sample #####

####### Alpha diversity indices######
plot_richness(run_RSF_rarefied,x="System",measures=c("Observed","Shannon","Simpson")) + 
  geom_boxplot() +
  ylab("Alpha Diversity") + 
  theme(axis.text.x=element_text(angle=90,hjust=0.5,size=10),
        axis.title = element_text(size=18),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=20,face="bold"))


rich = estimate_richness(run_RSF_rarefied)
write.csv(rich, 'richness_RSF_Outside_samples.csv')

pairwise.wilcox.test(rich$Observed, sample_data(run_RSF_rarefied)$System)
pairwise.wilcox.test(rich$Shannon, sample_data(run_RSF_rarefied)$System)
pairwise.wilcox.test(rich$Simpson, sample_data(run_RSF_rarefied)$System)

################# Beta diversity
dm_weighted_unifrac <- phyloseq::distance(run_RSF_rarefied, method = "wUniFrac")
ordWU <- ordinate(run_RSF_rarefied,method='PCoA',distance=dm_weighted_unifrac) 
#### Plot
plot_ordination(run_RSF_rarefied, ordWU,,color = 'System') + 
  geom_point(size=5) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ggtitle("PCoA: Weighted Unifrac")

dm_unweighted_unifrac <- phyloseq::distance(run_RSF_rarefied, method='Unifrac')
ordUU <- ordinate(run_RSF_rarefied,method='PCoA',distance=dm_unweighted_unifrac) 
#### Plot
plot_ordination(run_RSF_rarefied, ordUU,,color = 'System') + 
  geom_point(size=5) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ggtitle("PCoA: Unweighted Unifrac") 

adonis(dm_weighted_unifrac ~ sample_data(run_RSF_rarefied)$System)
adonis(dm_unweighted_unifrac ~ sample_data(run_RSF_rarefied)$System)

bray_diss = phyloseq::distance(run_RSF_rarefied, method="bray")
ordination = ordinate(run_RSF_rarefied, method="PCoA", distance=bray_diss)
#### Plot
plot_ordination(run_RSF_rarefied, ordination, color = 'System') + theme(aspect.ratio=1)+
  geom_point(size=5) + 
  theme(axis.text.x=element_text(size=14,angle=0,hjust =1),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ggtitle("PCoA: Bray-Curtis")

adonis(bray_diss ~ sample_data(run_RSF_rarefied)$System)

#################### Core Microbiome of phyllosphere on tomato grown in Greenhouse
######### Hydroponics tomato plants 
RSF_rare_HYD <- subset_samples(run_RSF_rarefied, System == "Hydroponic")
RSF_rare_Hydro <- prune_taxa(taxa_sums(RSF_rare_HYD) > 0, RSF_rare_HYD)
summary(sample_data(RSF_rare_Hydro)$System)
print(RSF_rare_Hydro)
sample_sums(RSF_rare_Hydro)
## RESULTS: 
###### 29 samples
###### 7216 taxa
###### Samples were rarefied to 2546 reads per sample

## keep only taxa with positive sums
Hydro_rare.1 <- prune_taxa(taxa_sums(RSF_rare_Hydro) > 0, RSF_rare_Hydro)

## Relative abundances
Hydro_rare.rel <- microbiome::transform(Hydro_rare.1, "compositional")

## Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(Hydro_rare.rel, detection = 1, sort = TRUE))

## Phyloseq object of the core microbiota:
Hydro.core <- core(Hydro_rare.rel, detection = 0, prevalence = .5)

## Retrieving the associated taxa names from the phyloseq object:
Hydro_core.taxa <- taxa(Hydro.core)
class(Hydro_core.taxa)
## get the taxonomy data
Hydrotax.mat <- tax_table(Hydro.core)
Hydrotax.df <- as.data.frame(Hydrotax.mat)

## add the OTus to last column
Hydrotax.df$OTU <- rownames(Hydrotax.df)

## select taxonomy of only 
## those OTUs that are core members based on the thresholds that were used.
Hydro_core.taxa.class <- dplyr::filter(Hydrotax.df, rownames(Hydrotax.df) %in% Hydro_core.taxa)
knitr::kable(head(Hydro_core.taxa.class))

###############cCore heatmaps
## Core with absolute counts and vertical view:
## and minimum population prevalence (given as percentage)
## Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-2), log10(max(abundances(Hydro_rare.1))/10), length = 10)

Hydrocore <- plot_core(Hydro_rare.1, plot.type = "heatmap", 
                       prevalences = prevalences,
                       detections = detections,
                       colours = rev(brewer.pal(5, "Spectral")),
                       min.prevalence = .7)

### Data used for plotting 
Hydrodf <- Hydrocore$data 
### list of OTUs
Hydrolist <- Hydrodf$Taxa 
## check the OTU ids
print(Hydrolist) 

## Taxonomy data
Hydrotax <- tax_table(Hydro_rare.1)
Hydrotax <- as.data.frame(Hydrotax)
## Add the OTus to last column
Hydrotax$OTU <- rownames(Hydrotax)
## select taxonomy of only 
## those OTUs that are used in the plot
Hydrotax2 <- dplyr::filter(Hydrotax, rownames(Hydrotax) %in% Hydrolist) 

## Merge all the column into one except the Domain as all is bacteria in this case
hydrotax.unit <- tidyr::unite(Hydrotax2, Taxa_level,c("Family", "Genus","OTU"), sep = "_;", remove = TRUE)
hydrotax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", hydrotax.unit$Taxa_level)

## add this new information into the plot data df
Hydrodf$Taxa <- hydrotax.unit$Taxa_level
## Taxonomic information
knitr::kable(head(Hydrodf))
## replace the data in the plot object
Hydrocore$data <- Hydrodf

## Detection Threshold is the Relative Abundance in %
plot(Hydrocore + theme(axis.text.x=element_text(size=12,angle=0,hjust =1),
                       axis.text.y = element_text(size = 12,face="italic"),
                       strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
                       strip.text.y = element_text(size=18, face = 'bold'),
                       plot.title = element_text(size = rel(2)),
                       axis.title=element_text(size=18,face="bold", vjust = 10),
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +
       ggtitle("Core Hydroponic system")+
       xlab("Detection threshold (Relative abundance %) "))


######### Organics tomato plants 
RSF_rare_ORGA <- subset_samples(run_RSF_rarefied, System == "Organic") 
RSF_rare_Soil <- prune_taxa(taxa_sums(RSF_rare_ORGA) > 0, RSF_rare_ORGA)
summary(sample_data(RSF_rare_Soil)$System)
print(RSF_rare_Soil)
sample_sums(RSF_rare_Soil)
## RESULTS: 
###### 18 samples
###### 3628 taxa
###### Samples were rarefied to 2546 reads per sample

## keep only taxa with positive sums
Soil_rare.1 <- prune_taxa(taxa_sums(RSF_rare_Soil) > 0, RSF_rare_Soil)
## Relative abundances
Soil_rare.rel <- microbiome::transform(Soil_rare.1, "compositional")
## Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(Soil_rare.rel, detection = 1, sort = TRUE))
## Phyloseq object of the core microbiota:
Soil.core <- core(Soil_rare.rel, detection = 0, prevalence = .5)

##Retrieving the associated taxa names from the phyloseq object:
Soil_core.taxa <- taxa(Soil.core)
class(Soil_core.taxa)
## get the taxonomy data
Soiltax.mat <- tax_table(Soil.core)
Soiltax.df <- as.data.frame(Soiltax.mat)
## add the OTus to last column
Soiltax.df$OTU <- rownames(Soiltax.df)

## select taxonomy of only 
## those OTUs that are core members based on the thresholds that were used.
Soil_core.taxa.class <- dplyr::filter(Soiltax.df, rownames(Soiltax.df) %in% Soil_core.taxa)
knitr::kable(head(Soil_core.taxa.class))

####Core heatmaps
#Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-2), log10(max(abundances(Soil_rare.1))/10), length = 10)

Soilcore <- plot_core(Soil_rare.1, plot.type = "heatmap", 
                      prevalences = prevalences,
                      detections = detections,
                      colours = rev(brewer.pal(5, "Spectral")),
                      min.prevalence = .7)

## Data used for plotting 
Soildf <- Soilcore$data 
## list of OTUs
Soillist <- Soildf$Taxa 
##check the OTU ids
print(Soillist) 
## Taxonomy data
Soiltax <- tax_table(Soil_rare.1)
Soiltax <- as.data.frame(Soiltax)
## Add the OTus to last column
Soiltax$OTU <- rownames(Soiltax)
## select taxonomy of only 
## those OTUs that are used in the plot
Soiltax2 <- dplyr::filter(Soiltax, rownames(Soiltax) %in% Soillist) 
## Merge all the column into one except the Domain as all is bacteria in this case
Soiltax.unit <- tidyr::unite(Soiltax2, Taxa_level,c("Family", "Genus","OTU"), sep = "_;", remove = TRUE)

Soiltax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", Soiltax.unit$Taxa_level)
## add this new information into the plot data df
Soildf$Taxa <- Soiltax.unit$Taxa_level
## Taxonomic information
knitr::kable(head(Soildf))
## replace the data in the plot object
Soilcore$data <- Soildf
## Detection Threshold is the Relative Abundance in %
plot(Soilcore + theme(axis.text.x=element_text(size=12,angle=0,hjust =1),
                      axis.text.y = element_text(size = 12,face="italic"),
                      strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
                      strip.text.y = element_text(size=18, face = 'bold'),
                      plot.title = element_text(size = rel(2)),
                      axis.title=element_text(size=18,face="bold", vjust = 10),
                      plot.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
       ggtitle("Core Soil system")+
       xlab("Detection threshold (Relative abundance %) "))


############### Tmato Plants exposed to Rain 
RSF_rare_OUT_R <- subset_samples(run_RSF_rarefied, System == "Out.Rain") 
RSF_rare_Out.Rain <- prune_taxa(taxa_sums(RSF_rare_OUT_R) > 0, RSF_rare_OUT_R)
summary(sample_data(RSF_rare_Out.Rain)$System)
print(RSF_rare_Out.Rain)
sample_sums(RSF_rare_Out.Rain)
## RESULTS: 
###### 7 samples
###### 2422 taxa
###### Samples were rarefied to 2546 reads per sample

## keep only taxa with positive sums
Out.Rain_rare.1 <- prune_taxa(taxa_sums(RSF_rare_Out.Rain) > 0, RSF_rare_Out.Rain)
## Relative abundances
Out.Rain_rare.rel <- microbiome::transform(Out.Rain_rare.1, "compositional")
## Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(Out.Rain_rare.rel, detection = 1, sort = TRUE))
## Phyloseq object of the core microbiota:
Out.Rain.core <- core(Out.Rain_rare.rel, detection = 0, prevalence = .5)

## Retrieving the associated taxa names from the phyloseq object:
Out.Rain_core.taxa <- taxa(Out.Rain.core)
class(Out.Rain_core.taxa)
## get the taxonomy data
Out.Raintax.mat <- tax_table(Out.Rain.core)
Out.Raintax.df <- as.data.frame(Out.Raintax.mat)

## add the OTus to last column
Out.Raintax.df$OTU <- rownames(Out.Raintax.df)
## select taxonomy of only 
## those OTUs that are core members based on the thresholds that were used.
Out.Rain_core.taxa.class <- dplyr::filter(Out.Raintax.df, rownames(Out.Raintax.df) %in% Out.Rain_core.taxa)
knitr::kable(head(Out.Rain_core.taxa.class))

##### Core heatmaps
## Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-2), log10(max(abundances(Out.Rain_rare.1))/10), length = 10)

Out.Raincore <- plot_core(Out.Rain_rare.1, plot.type = "heatmap", 
                          prevalences = prevalences,
                          detections = detections,
                          colours = rev(brewer.pal(5, "Spectral")),
                          min.prevalence = .7)

## Data used for plotting 
Out.Raindf <- Out.Raincore$data 
## list of OTUs
Out.Rainlist <- Out.Raindf$Taxa 
## check the OTU ids
print(Out.Rainlist) 
## Taxonomy data
Out.Raintax <- tax_table(Out.Rain_rare.1)
Out.Raintax <- as.data.frame(Out.Raintax)
## Add the OTus to last column
Out.Raintax$OTU <- rownames(Out.Raintax)
## select taxonomy of only 
## those OTUs that are used in the plot
Out.Raintax2 <- dplyr::filter(Out.Raintax, rownames(Out.Raintax) %in% Out.Rainlist) 

## Merge all the column into one except the Domain as all is bacteria in this case
Out.Raintax.unit <- tidyr::unite(Out.Raintax2, Taxa_level,c("Family", "Genus","OTU"), sep = "_;", remove = TRUE)
Out.Raintax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", Out.Raintax.unit$Taxa_level)
## add this new information into the plot data df
Out.Raindf$Taxa <- Out.Raintax.unit$Taxa_level
## Taxonomic information
knitr::kable(head(Out.Raindf))
## replace the data in the plot object
Out.Raincore$data <- Out.Raindf
## Detection Threshold is the Relative Abundance in %
plot(Out.Raincore + theme(axis.text.x=element_text(size=12,angle=0,hjust =1),
                          axis.text.y = element_text(size = 12,face="italic"),
                          strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
                          strip.text.y = element_text(size=18, face = 'bold'),
                          plot.title = element_text(size = rel(2)),
                          axis.title=element_text(size=18,face="bold", vjust = 10),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
                    ggtitle("Core Tomatoes exposed to Rain")+
                    xlab("Detection threshold (Relative abundance %) "))

######### Tomato Plants NOT exposed to Rain 
RSF_rare_Out.NRain <- subset_samples(run_RSF_rarefied, System == "Out.No.Rain") 
RSF_rare_Out.No.Rain <- prune_taxa(taxa_sums(RSF_rare_Out.NRain) > 0, RSF_rare_Out.NRain)
summary(sample_data(RSF_rare_Out.No.Rain)$System)
print(RSF_rare_Out.No.Rain)
sample_sums(RSF_rare_Out.No.Rain)
## RESULTS: 
###### 3 samples
###### 1312 taxa
###### Samples were rarefied to 2546 reads per sample

## keep only taxa with positive sums
Out.No.Rain_rare.1 <- prune_taxa(taxa_sums(RSF_rare_Out.No.Rain) > 0, RSF_rare_Out.No.Rain)

## Relative abundances
Out.No.Rain_rare.rel <- microbiome::transform(Out.No.Rain_rare.1, "compositional")
## Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(Out.No.Rain_rare.rel, detection = 1, sort = TRUE))
## Phyloseq object of the core microbiota:
Out.No.Rain.core <- core(Out.No.Rain_rare.rel, detection = 0, prevalence = .5)
## Retrieving the associated taxa names from the phyloseq object:
Out.No.Rain_core.taxa <- taxa(Out.No.Rain.core)
class(Out.No.Rain_core.taxa)
## get the taxonomy data
Out.No.Raintax.mat <- tax_table(Out.No.Rain.core)
Out.No.Raintax.df <- as.data.frame(Out.No.Raintax.mat)

## add the OTus to last column
Out.No.Raintax.df$OTU <- rownames(Out.No.Raintax.df)

## select taxonomy of only 
## those OTUs that are core members based on the thresholds that were used.
Out.No.Rain_core.taxa.class <- dplyr::filter(Out.No.Raintax.df, rownames(Out.No.Raintax.df) %in% Out.No.Rain_core.taxa)
knitr::kable(head(Out.No.Rain_core.taxa.class))

###### Core heatmaps
## Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-1), log10(max(abundances(Out.No.Rain_rare.1))/10), length = 10)

Out.No.Raincore <- plot_core(Out.No.Rain_rare.1, plot.type = "heatmap", 
                             prevalences = prevalences,
                             detections = detections,
                             colours = rev(brewer.pal(5, "Spectral")),
                             min.prevalence = .9)

## Data used for plotting 
Out.No.Raindf <- Out.No.Raincore$data 
## list of OTUs
Out.No.Rainlist <- Out.No.Raindf$Taxa 
## check the OTU ids
print(Out.No.Rainlist) 
## Taxonomy data
Out.No.Raintax <- tax_table(Out.No.Rain_rare.1)
Out.No.Raintax <- as.data.frame(Out.No.Raintax)
## Add the OTus to last column
Out.No.Raintax$OTU <- rownames(Out.No.Raintax)
## select taxonomy of only 
## those OTUs that are used in the plot
Out.No.Raintax2 <- dplyr::filter(Out.No.Raintax, rownames(Out.No.Raintax) %in% Out.No.Rainlist) 

## Merge all the column into one except the Domain as all is bacteria in this case
Out.No.Raintax.unit <- tidyr::unite(Out.No.Raintax2, Taxa_level,c("Family", "Genus","OTU"), sep = "_;", remove = TRUE)

Out.No.Raintax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", Out.No.Raintax.unit$Taxa_level)
## add this new information into the plot data df
Out.No.Raindf$Taxa <- Out.No.Raintax.unit$Taxa_level
## Taxonomic information
knitr::kable(head(Out.No.Raindf))
## replace the data in the plot object
Out.No.Raincore$data <- Out.No.Raindf
## Detection Threshold is the Relative Abundance in %
plot(Out.No.Raincore + theme(axis.text.x=element_text(size=12,angle=0,hjust =1),
                             axis.text.y = element_text(size = 12,face="italic"),
                             strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
                             strip.text.y = element_text(size=18, face = 'bold'),
                             plot.title = element_text(size = rel(2)),
                             axis.title=element_text(size=18,face="bold", vjust = 10),
                             plot.background = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank()) +
                      ggtitle("Core Tomatoes outside but not exposed to Rain")+
                      xlab("Detection threshold (Relative abundance %) "))


######### Tomato Plants grown outside: 7 exposed to rain and 3 covered from rain
RSF_rare_OUT <- subset_samples(run_RSF_rarefied, Time == "Outside") 
RSF_rare_Outside <- prune_taxa(taxa_sums(RSF_rare_OUT) > 0, RSF_rare_OUT)
summary(sample_data(RSF_rare_Outside)$System)
print(RSF_rare_Outside)
sample_sums(RSF_rare_Outside)
## RESULTS: 
###### 10 samples
###### 3249 taxa
###### Samples were rarefied to 2546 reads per sample

## keep only taxa with positive sums
RSF_rare_Outside.1 <- prune_taxa(taxa_sums(RSF_rare_Outside) > 0, RSF_rare_Outside)
## Relative abundances
RSF_rare_Outside.rel <- microbiome::transform(RSF_rare_Outside.1, "compositional")
## Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(RSF_rare_Outside.rel, detection = 1, sort = TRUE))
## Phyloseq object of the core microbiota:
RSF_rare_Outside.core <- core(RSF_rare_Outside.rel, detection = 0, prevalence = .5)

## Retrieving the associated taxa names from the phyloseq object:
RSF_rare_Outside_core.taxa <- taxa(RSF_rare_Outside.core)
class(RSF_rare_Outside_core.taxa)
## get the taxonomy data
RSF_rare_Outsidetax.mat <- tax_table(RSF_rare_Outside.core)
RSF_rare_Outsidetax.df <- as.data.frame(RSF_rare_Outsidetax.mat)
## add the OTus to last column
RSF_rare_Outsidetax.df$OTU <- rownames(RSF_rare_Outsidetax.df)

## select taxonomy of only 
## those OTUs that are core members based on the thresholds that were used.
RSF_rare_Outside_core.taxa.class <- dplyr::filter(RSF_rare_Outsidetax.df, rownames(RSF_rare_Outsidetax.df) %in% RSF_rare_Outside_core.taxa)
knitr::kable(head(RSF_rare_Outside_core.taxa.class))

######### Core heatmaps
## Set different detection levels and prevalence
prevalences <- seq(.5, 1, .5) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-3), log10(max(abundances(RSF_rare_Outside.1))/10), length = 10)

RSF_rare_Outsidecore <- plot_core(RSF_rare_Outside.1, plot.type = "heatmap", 
                                  prevalences = prevalences,
                                  detections = detections,
                                  colours = rev(brewer.pal(5, "Spectral")),
                                  min.prevalence = .7)

## Data used for plotting 
RSF_rare_Outsidedf <- RSF_rare_Outsidecore$data 
## list of OTUs
RSF_rare_Outsidelist <- RSF_rare_Outsidedf$Taxa 
## check the OTU ids
print(RSF_rare_Outsidelist) 
## Taxonomy data
RSF_rare_Outsidetax <- tax_table(RSF_rare_Outside.1)
RSF_rare_Outsidetax <- as.data.frame(RSF_rare_Outsidetax)
## Add the OTus to last column
RSF_rare_Outsidetax$OTU <- rownames(RSF_rare_Outsidetax)
## select taxonomy of only 
## those OTUs that are used in the plot
RSF_rare_Outsidetax2 <- dplyr::filter(RSF_rare_Outsidetax, rownames(RSF_rare_Outsidetax) %in% RSF_rare_Outsidelist) 

## Merge all the column into one except the Domain as all is bacteria in this case
RSF_rare_Outsidetax.unit <- tidyr::unite(RSF_rare_Outsidetax2, Taxa_level,c("Family", "Genus","OTU"), sep = "_;", remove = TRUE)
######## RESULT: 2 OTUs in the Core hydroponic microbiome

RSF_rare_Outsidetax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", RSF_rare_Outsidetax.unit$Taxa_level)

## add this new information into the plot data df
RSF_rare_Outsidedf$Taxa <- RSF_rare_Outsidetax.unit$Taxa_level
## Taxonomic information
knitr::kable(head(RSF_rare_Outsidedf))
## replace the data in the plot object
RSF_rare_Outsidecore$data <- RSF_rare_Outsidedf
## Detection Threshold is the Relative Abundance in %
plot(RSF_rare_Outsidecore + theme(axis.text.x=element_text(size=12,angle=0,hjust =1),
                                  axis.text.y = element_text(size = 12,face="italic"),
                                  strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
                                  strip.text.y = element_text(size=18, face = 'bold'),
                                  plot.title = element_text(size = rel(2)),
                                  axis.title=element_text(size=18,face="bold", vjust = 10),
                                  plot.background = element_blank(),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            ggtitle("Core Tomatoes grown outside")+
                            xlab("Detection threshold (Relative abundance %) "))


######################### DIFERENTIAL ABUNDANCE ##########################################
######### Hydroponic vs Organic
run_RSF_HYDR_ORG <- subset_samples(run, System%in%c("Hydroponic","Organic"))
run_RSF_H_O <- prune_taxa(taxa_sums(run_RSF_HYDR_ORG) > 0, run_RSF_HYDR_ORG)
summary(sample_data(run_RSF_H_O)$System)
print(run_RSF_H_O)
sample_sums(run_RSF_H_O)
## RESULTS: 
###### 47 samples
###### 16386 taxa

RSF_Hyd_Org <- prune_samples(sample_sums(run_RSF_H_O) > 500, run_RSF_H_O)
head(sample_data(RSF_Hyd_Org)$System, 25)
deseq_H_O = phyloseq_to_deseq2(RSF_Hyd_Org, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_H_O), 1, gm_mean)
deseq_H_O = estimateSizeFactors(deseq_H_O, geoMeans = geoMeans)
deseq_H_O = DESeq(deseq_H_O, fitType="local")

res_H_O = results(deseq_H_O)
res_H_O = res_H_O[order(res_H_O$padj, na.last=NA), ]
alpha = 0.01
sigtab_H_O = res_H_O[(res_H_O$padj < alpha), ]
sigtab_H_O = cbind(as(sigtab_H_O, "data.frame"), as(tax_table(RSF_Hyd_Org)[rownames(sigtab_H_O), ], "matrix"))
head(sigtab_H_O)
##### To write all OTUs that were significant different: positives and negatives
sigtab_H_O = sigtab_H_O[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab_H_O, 'DEseq_all_values_Hydro_vs_Soil.csv')
###### RESULT: 76 TAXA differentially abundant ##############
##### To subset positives values
posigtab_H_O = sigtab_H_O[sigtab_H_O[, "log2FoldChange"] > 0, ]
posigtab_H_O = posigtab_H_O[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_H_O, 'Differential_abundance_Hydroponics vs Organics  .csv')
###### RESULT: 62 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_H_O = subset(sigtab_H_O, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_H_O$log2FoldChange, sigtabgen_H_O$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_O$Phylum = factor(as.character(sigtabgen_H_O$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_H_O$log2FoldChange, sigtabgen_H_O$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_O$Genus = factor(as.character(sigtabgen_H_O$Genus), levels=names(x))
ggplot(sigtabgen_H_O, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=12)) +
  ylab("Genus \n") +
  xlab("Log2FoldChange") +
  ggtitle("DEseq Hydroponic vs Organic")


######### Hydroponics vs tomatoes exposed to rain
run_RSF_HYDRO_RAIN <- subset_samples(run, System%in%c("Hydroponic","Out.Rain"))
run_RSF_H_R <- prune_taxa(taxa_sums(run_RSF_HYDRO_RAIN) > 0, run_RSF_HYDRO_RAIN)
summary(sample_data(run_RSF_H_R)$System)
print(run_RSF_H_R)
sample_sums(run_RSF_H_R)
## RESULTS: 
###### 36 samples
###### 19924 taxa

run_RSF_H_R <- prune_samples(sample_sums(run_RSF_H_R) > 500, run_RSF_H_R)
head(sample_data(run_RSF_H_R)$System, 25)
deseq_H_OutR = phyloseq_to_deseq2(run_RSF_H_R, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_H_OutR), 1, gm_mean)
deseq_H_OutR = estimateSizeFactors(deseq_H_OutR, geoMeans = geoMeans)
deseq_H_OutR = DESeq(deseq_H_OutR, fitType="local")

res_H_OutR = results(deseq_H_OutR)
res_H_OutR = res_H_OutR[order(res_H_OutR$padj, na.last=NA), ]
alpha = 0.01
sigtab_H_OutR = res_H_OutR[(res_H_OutR$padj < alpha), ]
sigtab_H_OutR = cbind(as(sigtab_H_OutR, "data.frame"), as(tax_table(run_RSF_H_R)[rownames(sigtab_H_OutR), ], "matrix"))
head(sigtab_H_OutR)
##### To write all OTUs that were significant different: positives and negatives
sigtab_H_OutR = sigtab_H_OutR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(sigtab_H_OutR, 'DEseq_all_values_Hydro_vs_OutRain.csv')
###### RESULT: 80 TAXA differentially abundant ##############
##### To subset positives values
posigtab_H_OutR = sigtab_H_OutR[sigtab_H_OutR[, "log2FoldChange"] > 0, ]
posigtab_H_OutR = posigtab_H_OutR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_H_OutR, 'Differential_abundance_Hydroponics vs_tomato_outsideRain.csv')
###### RESULT: 40 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_H_OutR = subset(sigtab_H_OutR, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_H_OutR$log2FoldChange, sigtabgen_H_OutR$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_OutR$Phylum = factor(as.character(sigtabgen_H_OutR$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_H_OutR$log2FoldChange, sigtabgen_H_OutR$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_OutR$Genus = factor(as.character(sigtabgen_H_OutR$Genus), levels=names(x))
ggplot(sigtabgen_H_OutR, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=12)) +
  ylab("Genus \n") +
  xlab("Log2FoldChange") +
  ggtitle("DEseq Hydroponic vs Outside_Rain")


######### Organics vs tomatoes exposed to rain
run_RSF_ORG_RAIN <- subset_samples(run, System%in%c("Organic","Out.Rain"))
run_RSF_O_R <- prune_taxa(taxa_sums(run_RSF_ORG_RAIN) > 0, run_RSF_ORG_RAIN)
summary(sample_data(run_RSF_O_R)$System)
print(run_RSF_O_R)
sample_sums(run_RSF_O_R)
## RESULTS: 
###### 25 samples
###### 17305 taxa

run_RSF_O_R <- prune_samples(sample_sums(run_RSF_O_R) > 500, run_RSF_O_R)
head(sample_data(run_RSF_O_R)$System, 25)
deseq_O_OutR = phyloseq_to_deseq2(run_RSF_O_R, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_O_OutR), 1, gm_mean)
deseq_O_OutR = estimateSizeFactors(deseq_O_OutR, geoMeans = geoMeans)
deseq_O_OutR = DESeq(deseq_O_OutR, fitType="local")

res_O_OutR = results(deseq_O_OutR)
res_O_OutR = res_O_OutR[order(res_O_OutR$padj, na.last=NA), ]
alpha = 0.01
sigtab_O_OutR = res_O_OutR[(res_O_OutR$padj < alpha), ]
sigtab_O_OutR = cbind(as(sigtab_O_OutR, "data.frame"), as(tax_table(run_RSF_O_R)[rownames(sigtab_O_OutR), ], "matrix"))
head(sigtab_O_OutR)
##### To write all OTUs that were significant different: positives and negatives
sigtab_O_OutR = sigtab_O_OutR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(sigtab_O_OutR, 'DEseq_all_values_Org_vs_OutRain.csv')
###### RESULT: 159 TAXA differentially abundant ##############
##### To subset positives values
posigtab_O_OutR = sigtab_O_OutR[sigtab_O_OutR[, "log2FoldChange"] > 0, ]
posigtab_O_OutR = posigtab_O_OutR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_O_OutR, 'Differential_abundance_Organics vs_tomato_outsideRain.csv')
###### RESULT: 58 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_O_OutR = subset(sigtab_O_OutR, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_O_OutR$log2FoldChange, sigtabgen_O_OutR$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_O_OutR$Phylum = factor(as.character(sigtabgen_O_OutR$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_O_OutR$log2FoldChange, sigtabgen_O_OutR$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_O_OutR$Genus = factor(as.character(sigtabgen_O_OutR$Genus), levels=names(x))
ggplot(sigtabgen_O_OutR, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=12)) +
  ylab("Genus \n") +
  xlab("Log2FoldChange") +
  ggtitle("DEseq Organic vs Outside_Rain")



######### Hydroponics vs tomatoes outside  not exposed to rain
run_RSF_Hydro_NR <- subset_samples(run, System%in%c("Hydroponic","Out.No.Rain"))
run_RSF_H_NR <- prune_taxa(taxa_sums(run_RSF_Hydro_NR) > 0, run_RSF_Hydro_NR)
summary(sample_data(run_RSF_H_NR)$System)
print(run_RSF_H_NR)
sample_sums(run_RSF_H_NR)
## RESULTS: 
###### 32 samples
###### 17625 taxa

run_RSF_H_NR <- prune_samples(sample_sums(run_RSF_H_NR) > 500, run_RSF_H_NR)
head(sample_data(run_RSF_H_NR)$System, 25)
deseq_H_OutNR = phyloseq_to_deseq2(run_RSF_H_NR, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_H_OutNR), 1, gm_mean)
deseq_H_OutNR = estimateSizeFactors(deseq_H_OutNR, geoMeans = geoMeans)
deseq_H_OutNR = DESeq(deseq_H_OutNR, fitType="local")

res_H_OutNR = results(deseq_H_OutNR)
res_H_OutNR = res_H_OutNR[order(res_H_OutNR$padj, na.last=NA), ]
alpha = 0.01
sigtab_H_OutNR = res_H_OutNR[(res_H_OutNR$padj < alpha), ]
sigtab_H_OutNR = cbind(as(sigtab_H_OutNR, "data.frame"), as(tax_table(run_RSF_H_NR)[rownames(sigtab_H_OutNR), ], "matrix"))
head(sigtab_H_OutNR)
##### To write all OTUs that were significant different: positives and negatives
sigtab_H_OutNR = sigtab_H_OutNR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(sigtab_H_OutNR, 'DEseq_all_values_Hydro_vs_OutNORain.csv')
###### RESULT: 21 TAXA differentially abundant ##############
##### To subset positives values
posigtab_H_OutNR = sigtab_H_OutNR[sigtab_H_OutNR[, "log2FoldChange"] > 0, ]
posigtab_H_OutNR = posigtab_H_OutNR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_H_OutNR, 'Differential_abundance_Hydroponics vs_tomato_outsideNoRain.csv')
###### RESULT: 14 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_H_OutNR = subset(sigtab_H_OutNR, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_H_OutNR$log2FoldChange, sigtabgen_H_OutNR$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_OutNR$Phylum = factor(as.character(sigtabgen_H_OutNR$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_H_OutNR$log2FoldChange, sigtabgen_H_OutNR$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_OutNR$Genus = factor(as.character(sigtabgen_H_OutNR$Genus), levels=names(x))
ggplot(sigtabgen_H_OutNR, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=12)) +
  ylab("Genus \n") +
  xlab("Log2FoldChange") +
  ggtitle("DEseq Hydroponic vs Outside_No_Rain")


######### Organics vs tomatoes not exposed to rain
run_RSF_Org_NR <- subset_samples(run, System%in%c("Organic","Out.No.Rain"))
run_RSF_O_NR <- prune_taxa(taxa_sums(run_RSF_Org_NR) > 0, run_RSF_Org_NR)
summary(sample_data(run_RSF_O_NR)$System)
print(run_RSF_O_NR)
sample_sums(run_RSF_O_NR)
## RESULTS: 
###### 21 samples
###### 15111 taxa

run_RSF_O_NR <- prune_samples(sample_sums(run_RSF_O_NR) > 500, run_RSF_O_NR)
head(sample_data(run_RSF_O_NR)$System, 25)
deseq_O_OutNR = phyloseq_to_deseq2(run_RSF_O_NR, ~ System)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_O_OutNR), 1, gm_mean)
deseq_O_OutNR = estimateSizeFactors(deseq_O_OutNR, geoMeans = geoMeans)
deseq_O_OutNR = DESeq(deseq_O_OutNR, fitType="local")

res_O_OutNR = results(deseq_O_OutNR)
res_O_OutNR = res_O_OutNR[order(res_O_OutNR$padj, na.last=NA), ]
alpha = 0.01
sigtab_O_OutNR = res_O_OutNR[(res_O_OutNR$padj < alpha), ]
sigtab_O_OutNR = cbind(as(sigtab_O_OutNR, "data.frame"), as(tax_table(run_RSF_O_NR)[rownames(sigtab_O_OutNR), ], "matrix"))
head(sigtab_O_OutNR)
##### To write all OTUs that were significant different: positives and negatives
sigtab_O_OutNR = sigtab_O_OutNR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(sigtab_O_OutNR, 'DEseq_all_values_Org_vs_OutNORain.csv')
###### RESULT: 91 TAXA differentially abundant ##############
##### To subset positives values
posigtab_O_OutNR = sigtab_O_OutNR[sigtab_O_OutNR[, "log2FoldChange"] > 0, ]
posigtab_O_OutNR = posigtab_O_OutNR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_O_OutNR, 'Differential_abundance_Organics vs_tomato_outsideNoRain.csv')
###### RESULT: 49 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_O_OutNR = subset(sigtab_O_OutNR, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_O_OutNR$log2FoldChange, sigtabgen_O_OutNR$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_O_OutNR$Phylum = factor(as.character(sigtabgen_O_OutNR$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_O_OutNR$log2FoldChange, sigtabgen_O_OutNR$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_O_OutNR$Genus = factor(as.character(sigtabgen_O_OutNR$Genus), levels=names(x))
ggplot(sigtabgen_O_OutNR, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=12)) +
  ylab("Genus \n") +
  xlab("Log2FoldChange") +
  ggtitle("DEseq Organic vs Outside_No_Rain")


######### Hydroponics vs tomatoes outside
run_RSF_Hydro_Outside <- subset_samples(run, Source%in%c("Hydroponic","Outside"))
run_RSF_H_Out <- prune_taxa(taxa_sums(run_RSF_Hydro_Outside) > 0, run_RSF_Hydro_Outside)
summary(sample_data(run_RSF_H_Out)$System)
print(run_RSF_H_Out)
sample_sums(run_RSF_H_Out)
## RESULTS: 
###### 39 samples
###### 21245 taxa

run_RSF_H_Out <- prune_samples(sample_sums(run_RSF_H_Out) > 500, run_RSF_H_Out)
head(sample_data(run_RSF_H_Out)$Source, 25)

deseq_H_Out = phyloseq_to_deseq2(run_RSF_H_Out, ~ Source)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_H_Out), 1, gm_mean)
deseq_H_Out = estimateSizeFactors(deseq_H_Out, geoMeans = geoMeans)
deseq_H_Out = DESeq(deseq_H_Out, fitType="local")

res_H_Out = results(deseq_H_Out)
res_H_Out = res_H_Out[order(res_H_Out$padj, na.last=NA), ]
alpha = 0.01
sigtab_H_Out = res_H_Out[(res_H_Out$padj < alpha), ]
sigtab_H_Out = cbind(as(sigtab_H_Out, "data.frame"), as(tax_table(run_RSF_H_Out)[rownames(sigtab_H_Out), ], "matrix"))
head(sigtab_H_Out)
##### To write all OTUs that were significant different: positives and negatives
sigtab_H_Out = sigtab_H_Out[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(sigtab_H_Out, 'DEseq_all_values_Hydro_vs_Outside.csv')
###### RESULT: 179 TAXA differentially abundant ##############
##### To subset positives values
posigtab_H_Out = sigtab_H_Out[sigtab_H_Out[, "log2FoldChange"] > 0, ]
posigtab_H_Out = posigtab_H_Out[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_H_Out, 'Differential_abundance_Hydroponics vs outside.csv')
###### RESULT: 116 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_H_Out = subset(sigtab_H_Out, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_H_Out$log2FoldChange, sigtabgen_H_Out$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_Out$Phylum = factor(as.character(sigtabgen_H_Out$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_H_Out$log2FoldChange, sigtabgen_H_Out$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_H_Out$Genus = factor(as.character(sigtabgen_H_Out$Genus), levels=names(x))
ggplot(sigtabgen_H_Out, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=12)) +
  ylab("Genus \n") +
  xlab("Log2FoldChange") +
  ggtitle("DEseq Hydroponic vs Outside")


######### Organics vs tomatoes outside
run_RSF_Org_Outside <- subset_samples(run, Source%in%c("Organic","Outside"))
run_RSF_O_Out <- prune_taxa(taxa_sums(run_RSF_Org_Outside) > 0, run_RSF_Org_Outside)
summary(sample_data(run_RSF_O_Out)$System)
print(run_RSF_O_Out)
sample_sums(run_RSF_O_Out)
## RESULTS: 
###### 28 samples
###### 19195 taxa

run_RSF_O_Out <- prune_samples(sample_sums(run_RSF_O_Out) > 500, run_RSF_O_Out)
head(sample_data(run_RSF_O_Out)$System, 25)

deseq_O_Out = phyloseq_to_deseq2(run_RSF_O_Out, ~ Source)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_O_Out), 1, gm_mean)
deseq_O_Out = estimateSizeFactors(deseq_O_Out, geoMeans = geoMeans)
deseq_O_Out = DESeq(deseq_O_Out, fitType="local")

res_O_Out = results(deseq_O_Out)
res_O_Out = res_O_Out[order(res_O_Out$padj, na.last=NA), ]
alpha = 0.01
sigtab_O_Out = res_O_Out[(res_O_Out$padj < alpha), ]
sigtab_O_Out = cbind(as(sigtab_O_Out, "data.frame"), as(tax_table(run_RSF_O_Out)[rownames(sigtab_O_Out), ], "matrix"))
head(sigtab_O_Out)
##### To write all OTUs that were significant different: positives and negatives
sigtab_O_Out = sigtab_O_Out[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(sigtab_O_Out, 'DEseq_all_values_Org_vs_Outside.csv')
###### RESULT: 255 TAXA differentially abundant ##############
##### To subset positives values
posigtab_O_Out = sigtab_O_Out[sigtab_O_Out[, "log2FoldChange"] > 0, ]
posigtab_O_Out = posigtab_O_Out[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
write.csv(posigtab_O_Out, 'Differential_abundance_Organics vs outside.csv')
###### RESULT: 129 TAXA differentially abundant ##############

theme_set(theme_bw())
sigtabgen_O_Out = subset(sigtab_O_Out, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen_O_Out$log2FoldChange, sigtabgen_O_Out$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_O_Out$Phylum = factor(as.character(sigtabgen_O_Out$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen_O_Out$log2FoldChange, sigtabgen_O_Out$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen_O_Out$Genus = factor(as.character(sigtabgen_O_Out$Genus), levels=names(x))
ggplot(sigtabgen_O_Out, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =1,  vjust=0.5),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=12)) +
  ylab("Genus \n") +
  xlab("Log2FoldChange") +
  ggtitle("DEseq Organic vs Outside")
