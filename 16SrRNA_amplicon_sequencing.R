

### Installing every package needed for the analysis.
### Need to also add phylosmith package. 

github.pkg = c("jfukuyama/phyloseqGraphTest", "jbisanz/qiime2R") 
bioc.pkg = c("phyloseq", "DESeq2", "decontam", "MicrobiotaProcess")
cran.pkg = c("tidyverse", "pacman", "glue", "vegan", "devtools", "ggrepel", "reshape2", 
             "ggnetwork", "DT", "intergraph","VennDiagram", "lsmeans", 
             "pheatmap", "phyloseqGraphTest")

inst.pkg = cran.pkg %in% installed.packages()
if (any(!inst.pkg)){ install.packages(cran.pkg[!inst.pkg],repos = "http://cran.rstudio.com/") } 

inst.pkg = github.pkg %in% installed.packages() 
if (any(!inst.pkg)){ devtools::install_github(github.pkg[!inst.pkg], force = TRUE) } 

inst.pkg = bioc.pkg %in% installed.packages() 
if(any(!inst.pkg)){ BiocManager::install(bioc.pkg[!inst.pkg])}

### Importing Qiime2 artifacts into a phyloseq object with qiime2R and phyloseq packages.

library("qiime2R")
library("phyloseq")
library("tidyverse")
# Setting up the working directory to import every necessary qiime artifact and a random seed for reproducibility
setwd("D:/TFM/obese_qzv")
set.seed(123)

pst = qza_to_phyloseq(features = "tableNoFilt.qza",
                      tree = "treeNoFilt.qza", 
                      taxonomy = "taxonomyNoFilt.qza",
                      metadata = "../batch_metadata.tsv")
# Merging the repseqs to the phyloseq object 
repseqs = read_qza("repseqsNoFilt.qza")$data
pst = merge_phyloseq(pst, repseqs)

# Extracting the refseqs into a FASTA file using writexstringSet function from Biostrings package
library(Biostrings)
writeXStringSet(refseq(pst), "./refseqs.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
# Converting chategorical data from the metadata file (accesed via sample_data() function from phyloseq) 
# into a factor (diet)
for(i in seq_len(ncol(sample_data(pst)))) {
  if(!is.numeric(sample_data(pst)[[i]]) && !is.logical(sample_data(pst)[[i]])) {
    sample_data(pst)[[i]] = as.factor(sample_data(pst)[[i]])  } else {
      sample_data(pst)[[i]]
    } 
}

# Removing unassigned/NA values from the Phylum taxonomic level.
pst = subset_taxa(pst, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "unassigned"))

cat("Total number of non-assigned values on the Phylum taxonomic level is : ", 
    sum(is.na(tax_table(pst)[,2]) + (tax_table(pst)[,2] == "") + 
          (tax_table(pst)[,2] == "uncharacterized") +
          (tax_table(pst)[,2] == "unassigned") + 
          (tax_table(pst)[,2] == "undefined")))


# Defining out.ASV function for removing singletones based on abundance
out.ASV = function(phyloseq, threshold =1, binwidth = 0.01) {
  
  #Loading necessary pkgs      
  pacman::p_load(glue, tidyverse, reshape2, ggrepel, S4Vectors) # nolint
  #This function requires phyloseq, tidyverse and glue packages to be loaded. 
  if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 100 ) {#making the relative abundance table
    rel_abund = as(t(otu_table(phyloseq)), "matrix")
  } else if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 1) {
    rel_abund = as(t(otu_table(phyloseq)), "matrix")
  } else {
    rel_abund = as(t(apply(otu_table(phyloseq), 
                           ifelse(taxa_are_rows(phyloseq), 1,2), 
                           function(x) x/sum(x))), "matrix")  
  } 
  
  
  names.single = apply(rel_abund, 1, function(x){ifelse(x == threshold, TRUE, ifelse(x == sum(x),
                                                                                     TRUE, FALSE))}) %>% reshape2::melt() %>% filter(value == TRUE) %>% dplyr::select(2) %>%
    pull   %>% as.vector()
  
  
  if (length(names.single) == 0 ) {
    print(glue("WOW! {length(names.single)} singletones detected in this dataset"))
    qplot.noSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
                         show.legend = F, main = "Frequency count of relative abundance, no singletones detected") +
      xlab ("Relative abundance in samples") + ylab("Frequency") + theme_bw()
    
    
    return(structure(list(qplot.noSing)))
    
  } else { 
    
    single.ASV = rel_abund[rownames(rel_abund) %in% names.single,]
    single.ASV[single.ASV == 0] <- NA # A separate dataset for annotation of singletones on the barplot
    
    qplot.withSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
                           main = "Frequency count of relative abundance with singletones") +
      geom_bar(aes(single.ASV), fill = "red",  color = NA, width = binwidth)+
      xlab ("Relative abundance in samples") + ylab("Frequency") + 
      geom_label_repel(aes(x = 1, y =length(rel_abund)/5), 
                       label.padding =  unit(0.55, "lines"), 
                       label = glue("{length(names.single)}\n Singletones"), color = "black") + theme_bw()
    
    qplot.rmSing = qplot(rel_abund[!rownames(rel_abund) %in% names.single, ], geom = "histogram",
                         binwidth = binwidth, main = "Frequency count of relative abundance without singletones") +
      xlab ("Relative abundance in samples") + ylab("Frequency")+ theme_bw()
    
    print(glue('Oh no..! {length(names.single)} singletones detected in the dataset'))
    return(structure(list(qplot.withSing, qplot.rmSing, unlist(names.single))) )
    
  }                        
  
  
}
# Plotting ASV frequencies with and without singletones.
library(gridExtra)
single.test = out.ASV(phyloseq = pst, threshold = 1, binwidth = 0.1)
singletones = single.test[[3]] # names of the singletones
singletones.plot = single.test[[1]] # plot with singletones
no.singletones.plot = single.test[[2]] # plot without singletones
grid.arrange(singletones.plot, no.singletones.plot, ncol=2)


# Removing singletones from the dataset.
# We also keep the original phyloseq object WITH singletones for some analysis. 
ps = subset_taxa(pst, !taxa_names(pst)%in% singletones) # phyloseq without singletones 

# Removing non-bacterial ASVs from the Kingdom domain (Eukarya and Archaea ASVs)
pst = subset_taxa(pst, Kingdom="d_Bacteria")
ps = subset_taxa(ps, Kingdom="d_Bacteria")
pst
ps


# After removing the singletones around 400 taxa are lost. Those are stored in the singletones_phyloseq phyloseq object. 
singletones_phyloseq = subset_taxa(pst, taxa_names(pst)%in% singletones)


# Plotting rarefaction curves 
library(MicrobiotaProcess)
ps_rar_curve <- ggrarecurve(obj = ps, 
                            linesize = 0.9,
                            indexNames = c("Observe", "Shannon"),
                            chunks = 400,
                            theme = theme(legend.spacing.y = unit(0.02, "cm"),
                                          legend.text = element_text(size = 6)),
                            show.legend = FALSE)
pst_rar_curve <- ggrarecurve(obj = pst, 
                             linesize = 0.9,
                             indexNames = c("Observe", "Shannon"),
                             chunks = 400,
                             theme = theme(legend.spacing.y = unit(0.02, "cm"),
                                           legend.text = element_text(size = 6)),
                             show.legend = FALSE)
ps_rar_curve # without singletons
pst_rar_curve # with singletons 

# Rarefying at 300000 sampling depth for both datasets
ps.rar = phyloseq::rarefy_even_depth(ps, sample.size = 50000, replace = F)
pst.rar = phyloseq::rarefy_even_depth(pst, sample.size = 50000, replace = F)
ps.rar
pst.rar

library(MBECS)
mbec.obj <- mbecProcessInput(pst, required.col = c("batch", "diet"))
mbec.obj <- mbecTransform(mbec.obj, method = "clr", offset = 0.0001)
mbecReportPrelim(input.obj=mbec.obj, model.vars=c("batch","diet"), 
                 type="clr")
mbec.obj <- mbecCorrection(mbec.obj, model.vars=c("batch","diet"), 
                           method = "bat", type = "clr")


library(edgeR)

OTU = as(otu_table(pst), "matrix")
OTU = as.data.frame(OTU)
pheno
tax = as(tax_table(pst), "matrix")
tax = as.data.frame(tax)

edgeR <- DGEList(counts = OTU, samples = pheno, genes = tax) #Change metadata.micro for your own mt object
edgeR <- calcNormFactors(edgeR)
ASV_norm <- cpm(edgeR, normalized.lib.sizes=T, log=F)

phy_OTU_norm<-otu_table(as.data.frame(ASV_norm,row.names=F), taxa_are_rows = T)
phy_taxonomy_norm<-tax_table(as.matrix(tax))
phy_metadata_norm<-sample_data(pheno) #Change metadata.micro to your own mt file

##Add taxa names
taxa_names(phy_OTU_norm)<- taxa_names(phy_taxonomy_norm)
#Check
identical(rownames(ASV_norm), rownames(tax))
#> [1] TRUE

##Merge
norm_phyloseq<-phyloseq(phy_OTU_norm,phy_taxonomy_norm,phy_metadata_norm) #Remove tree_root when working with ITS


## COMBAT BATCH EFFECT. 
library(sva)



data.batch = as(otu_table(pst), "matrix")
data.batch = as.data.frame(data.batch)
pheno = data.frame(
  batch = sample_data(pst)$batch,
  diet = sample_data(pst)$diet
)
row.names(pheno) = row.names(sample_data(pst))
mod = model.matrix(~as.factor(diet), data=pheno)

adjusted_count = ComBat_seq(data.batch, pheno$batch, group=NULL, full_mod=FALSE)

combat = ComBat(dat=data.batch, batch=pheno$batch, mod=NULL, par.prior=FALSE,mean.only = TRUE)
OTU_combat<-otu_table(as.data.frame(adjusted_count), taxa_are_rows = T)
taxa_names(OTU_combat)<- taxa_names(phy_taxonomy_norm)
identical(rownames(OTU_combat), rownames(tax))
no.batch.pst = phyloseq(OTU_combat,phy_taxonomy_norm,phy_metadata_norm) #Remove tree_root when working with ITS

#LIMMA
# add a 1 offset to the count matrix for the clr transformation
library(mixOmics)
data.batch = data.batch + 1
data.batch.clr = logratio.transfo(data.batch, logratio = "CLR")

library(limma)

removeBatchEffect(x = data.batch,
                  batch = pheno$batch,
                  covariates = pheno$diet,
                  design = pheno)


### ALPHA DIVERSITY INDEX
# Using the estimate_richness() function from phyloseq package.
# Chao1 indecex
Chao1.rarefied.no.singletons =estimate_richness(ps.rar, split = TRUE, measures = "Chao1")
Chao1.rarefied.singletons =estimate_richness(pst.rar, split = TRUE, measures = "Chao1")
# Shannon index
Shannon.rarefied.no.singletons = estimate_richness(ps.rar, split = TRUE, measures = "Shannon")
Shannon.rarefied.singletons = estimate_richness(pst.rar, split = TRUE, measures = "Shannon")
# Simpson index
Simpson.rarefied.no.singletons = estimate_richness(ps.rar, split = TRUE, measures = "Simpson")
Simpson.rarefied.singletons = estimate_richness(pst.rar, split = TRUE, measures = "Simpson")
# Phylogenetic distance                   
library(picante)
FaithPD.rarefied.no.singletons = pd(t(otu_table(ps.rar)), tree = phy_tree(ps.rar), include.root = F)$PD
FaithPD.rarefied.singletons = pd(t(otu_table(pst.rar)), tree = phy_tree(pst.rar), include.root = F)$PD

# Creating dataframes with every alpha-index
sample_data(ps.rar) <- data.frame(sample_data(ps.rar), 
                                  Chao1=Chao1.rarefied.no.singletons[[1]],
                                  Shannon = Shannon.rarefied.no.singletons$Shannon,  
                                  FaithPD = FaithPD.rarefied.no.singletons, 
                                  Simpson = Simpson.rarefied.no.singletons$Simpson)  
alpha.df.rar.nosing <- sample_data(ps.rar)
alpha.df.rar.nosing

sample_data(pst.rar) <- data.frame(sample_data(pst.rar), 
                                   Chao1=Chao1.rarefied.singletons[[1]], 
                                   Shannon = Shannon.rarefied.singletons$Shannon,  
                                   FaithPD = FaithPD.rarefied.singletons, 
                                   Simpson = Simpson.rarefied.singletons$Simpson)  
alpha.df.rar.sing <- sample_data(pst.rar)
alpha.df.rar.sing

# Dividing the tables between the HFD experiment and the Alzheimers-induced experiment.

hfd.alpha.pst = subset(alpha.df.rar.sing, diet == "SD" | diet == "HFD" | diet=="HFD_PP")
alzheimers.alpha.pst = subset(alpha.df.rar.sing, diet == "SD" | diet == "S" | diet=="S_PP")
hfd.alpha.ps = subset(alpha.df.rar.nosing, diet == "SD" | diet == "HFD" | diet=="HFD_PP")
alzheimers.alpha.ps = subset(alpha.df.rar.nosing, diet == "SD" | diet == "S" | diet=="S_PP")


# We need to format this data into dataframes to check for significant differences with ANOVA. 
format_data = function(df) {

  data_for_anova = data.frame(
    Diet = df$diet,
    Shannon = df$Shannon,
    Chao1 = df$Chao1,
    FaithPD = df$FaithPD,
    Simpson = df$Simpson)
  return(data_for_anova)
}

hfd.alpha.pst.anova.data = format_data(hfd.alpha.pst)
hfd.alpha.ps.anova.data = format_data(hfd.alpha.ps)
alzheimers.alpha.ps.anova.data = format_data(alzheimers.alpha.ps)
alzheimers.alpha.pst.anova.data = format_data(alzheimers.alpha.pst)

# Checking for variances homogeneity with levene test to later apply ANOVAS
library(car)
levene.shannon.hfd = leveneTest(Shannon ~ Diet, data = hfd.alpha.pst.anova.data)
levene.chao1.hfd = leveneTest(Shannon ~ Diet, data = hfd.alpha.pst.anova.data)
levene.chao1.hfd = leveneTest(Shannon ~ Diet, data = hfd.alpha.pst.anova.data)


# ANOVAS AND POST-HOC TUKEY FOR HFD EXPERIMENT WITH SINGLETONES
hfd.alpha.pst.anova.shannon = aov(formula = Shannon ~ Diet, data = hfd.alpha.pst.anova.data)
hfd.alpha.pst.anova.chao1 = aov(formula = Chao1 ~ Diet, data = hfd.alpha.pst.anova.data)
hfd.alpha.pst.anova.faith = aov(formula = FaithPD ~ Diet, data = hfd.alpha.pst.anova.data)

# Every alpha-index metric shows significant results between the HFD-SD and HFD+PP-SD groups.
hfd.alpha.pst.tukey.shannon = TukeyHSD(hfd.alpha.pst.anova.shannon)
hfd.alpha.pst.tukey.chao1 = TukeyHSD(hfd.alpha.pst.anova.chao1)
hfd.alpha.pst.tukey.faith = TukeyHSD(hfd.alpha.pst.anova.faith)

# ANOVAS FOR HFD EXPERIMENT WITHOUT SINGLETONES
# Results are the same as with singletones. Data with singletones will be used not to lose biological variability.
hfd.alpha.ps.anova.shannon = aov(formula = Shannon ~ Diet, data = hfd.alpha.ps.anova.data)
hfd.alpha.ps.anova.chao1 = aov(formula = Chao1 ~ Diet, data = hfd.alpha.ps.anova.data)
hfd.alpha.ps.anova.faith = aov(formula = FaithPD ~ Diet, data = hfd.alpha.ps.anova.data)

# Repeating the process for the alzheimers experiment
levene.shannon.alzheimers = leveneTest(Shannon ~ Diet, data = alzheimers.alpha.pst.anova.data)
levene.chao1.alzheimers = leveneTest(Shannon ~ Diet, data = alzheimers.alpha.pst.anova.data)
levene.chao1.alzheimers = leveneTest(Shannon ~ Diet, data = alzheimers.alpha.pst.anova.data)


# ANOVAS AND POST-HOC TUKEY FOR HFD EXPERIMENT WITH SINGLETONES
alzheimers.alpha.pst.anova.shannon = aov(formula = Shannon ~ Diet, data = alzheimers.alpha.pst.anova.data)
alzheimers.alpha.pst.anova.chao1 = aov(formula = Chao1 ~ Diet, data = alzheimers.alpha.pst.anova.data)
alzheimers.alpha.pst.anova.faith = aov(formula = FaithPD ~ Diet, data = alzheimers.alpha.pst.anova.data)

# Every alpha-index metric shows significant results between the HFD-SD and HFD+PP-SD groups.
alzheimers.alpha.pst.tukey.shannon = TukeyHSD(alzheimers.alpha.pst.anova.shannon)
alzheimers.alpha.pst.tukey.chao1 = TukeyHSD(alzheimers.alpha.pst.anova.chao1)
alzheimers.alpha.pst.tukey.faith = TukeyHSD(alzheimers.alpha.pst.anova.faith)

TukeyHSD(aov(formula = Chao1 ~ Diet, data = alzheimers.alpha.ps.anova.data))


### GRAPHS FOR ALPHA DIVERSITY
library(viridis)

labels_hfd <- c(
  "SD" = "Standard\ndiet",
  "HFD" = "High fat\ndiet",
  "HFD+PP" = "High fat\ndiet+PP"
)


labels_AD <- c(
  "SD" = "SD",
  "S" = "Sco",
  "S_PP" = "Sco+PP"
)


boxplot.shannon.hfd = ggplot(hfd.alpha.pst.anova.data, aes(x = Diet, y = Shannon, fill = Diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = expression(bold("(A) Shannon index")),
       x = "",
       y = "Shannon") +
  scale_x_discrete(labels = labels_hfd) +
  scale_fill_manual(values=c("#238A8DFF","#DCE319FF","#73D055FF")) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold")) +
  ylim(3.25, 4.75)


boxplot.chao1.hfd = ggplot(hfd.alpha.pst.anova.data, aes(x = Diet, y = Chao1, fill = Diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = expression(bold("(B) Chao1 index")),
       x = "",
       y = "Chao1") +
  scale_x_discrete(labels = labels_hfd) +
  scale_fill_manual(values=c("#238A8DFF","#DCE319FF","#73D055FF")) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold"))

library(gridExtra)
grid.arrange(boxplot.shannon.hfd, boxplot.chao1.hfd, nrow=2)
devtools::install_github("kassambara/ggpubr")
library(ggpubr)

## TFG DE JOSELUIS:

alzheimers.alpha.pst.anova.data$Diet <- factor(alzheimers.alpha.pst.anova.data$Diet, levels = c("SD", "S", "S_PP"))


# check for normality of shannon and chao1 indexes of alzheimers.alpha.pst.anova.data before t.test
shapiro.test(alzheimers.alpha.pst.anova.data$Shannon)
shapiro.test(alzheimers.alpha.pst.anova.data$Chao1)
# both > 0.05, so we can use t.test

# plot AD shannon index boxplot with manual tukey p values using ggpubr
boxplot.shannon.AD = ggplot(alzheimers.alpha.pst.anova.data, aes(x = Diet, y = Shannon, fill = Diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title="",
       x = "",
       y = "Indice de Shannon") +
  scale_x_discrete(labels = labels_AD) +
  scale_fill_manual(values=c("white","gray8","chartreuse4")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 17),
        plot.title = element_text(size = 14, face = "bold")) +
  ylim(4, 4.75)
#add p values from t.test to the graph
boxplot.shannon.AD + stat_compare_means(comparisons = list(c("SD", "S"), c("SD", "S_PP"), c("S", "S_PP")), method = "t.test", label = "p.signif", size = 5)
pairwise.t.test(alzheimers.alpha.pst.anova.data$Shannon, alzheimers.alpha.pst.anova.data$Diet, p.adjust.method = "none")
# change order of boxplots. Now SD, S, S_PP



# plot AD chao1 index boxplot with manual tukey p values using ggpubr
boxplot.chao1.AD = ggplot(alzheimers.alpha.pst.anova.data, aes(x = Diet, y = Chao1, fill = Diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title="",
       x = "",
       y = "Indice Chao1") +
  scale_x_discrete(labels = labels_AD) +
  scale_fill_manual(values=c("black","darkgreen","white")) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold"))
# add p values from t.test to the graph
boxplot.chao1.AD + stat_compare_means(comparisons = list(c("SD", "S"), c("SD", "S_PP"), c("S", "S_PP")), method = "t.test", label = "p.signif", size = 5)
pairwise.t.test(alzheimers.alpha.pst.anova.data$Chao1, alzheimers.alpha.pst.anova.data$Diet, p.adjust.method = "none")

### ESTOS SON MIS VERSIONES::


boxplot.shannon.AD = ggplot(alzheimers.alpha.pst.anova.data, aes(x = Diet, y = Shannon, fill = Diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title="",
       x = "",
       y = "Shannon index") +
  scale_x_discrete(labels = labels_AD) +
  scale_fill_manual(values=c("black","darkgreen","white")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold")) +
  ylim(3.75, 4.75)



boxplot.chao1.AD = ggplot(alzheimers.alpha.pst.anova.data, aes(x = Diet, y = Chao1, fill = Diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "",
       x = "",
       y = "Chao1 index") +
  scale_x_discrete(labels = labels_AD) +
  scale_fill_manual(values=c("black","darkgreen","white")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold"))grid.arrange(boxplot.shannon.AD,boxplot.chao1.AD)

grid.arrange(boxplot.shannon.AD, boxplot.chao1.AD, nrow=2)


### DESEQ2 DIFFERENTIAL ABUNDANCE ANALYSIS
### 

# Loading DESeq2 package
library("DESeq2")
# The function phyloseq_to_deseq2() converts the phyloseq object into a DESeqDataSet with dispersions estimated, 
# using the experimental design formula. The DESeq() fucntion does the rest of the testing (see negative binomial zero inflated)
treatment.pst.rar = subset_samples(no.batch.pst, diet=="S" | diet=="S_PP")
treatment.pst.rar.deseq = phyloseq_to_deseq2(treatment.pst.rar, ~ diet)
# calculating geometric means prior to estimating size factors.
# this normalizes all samples to the same level of reads.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatment.pst.rar.deseq), 1, gm_mean)
treatment.pst.rar.deseq = estimateSizeFactors(treatment.pst.rar.deseq, geoMeans = geoMeans)
treatment.pst.rar.deseq = DESeq(treatment.pst.rar.deseq, fitType="local")

AD.deseq.results = results(treatment.pst.rar.deseq)
AD.deseq.results = AD.deseq.results[order(AD.deseq.results$padj, na.last=NA),]
alpha = 0.01
AD.signif = AD.deseq.results[(AD.deseq.results$padj < alpha), ]
AD.signif = cbind(as(AD.signif, "data.frame"), as(tax_table(treatment.pst.rar)[rownames(AD.signif), ], "matrix"))
head(AD.signif)
# Enriched ASVs in the pectine treatment group are those with positive logfold change:
positive.AD.signif = AD.signif[AD.signif[, "log2FoldChange"] > 15, ]
positive.AD.signif = positive.AD.signif[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
positive.AD.signif # Lachnospiraceae most enriched Family, log2fold 23.3, related to SCFA production.  

negative.AD.signif = AD.signif[AD.signif[, "log2FoldChange"] < -15, ]
negative.AD.signif = negative.AD.signif[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
negative.AD.signif # Lachnospiraceae most enriched Family, log2fold 23.3, related to SCFA production.

AD.signif.lfc = rbind(positive.AD.signif, negative.AD.signif)
AD.signif.lfc = AD.signif.lfc %>% select(log2FoldChange, padj, Phylum, Class, Genus) #keeping only those columns 
write.table(x = AD.signif.lfc, file = "DESEQ2_AD_signif_ASV.txt", sep = "\t")

# Set up the plot
library(ggplot2)
library(viridis)

log2FoldChange_threshold <- 15
padj_threshold <- 0.001

# Crear el volcano plot
library(ggrepel)

# Crear el volcano plot
volcano.plot.AD <- ggplot(AD.deseq.results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(AD.deseq.results, abs(log2FoldChange) < log2FoldChange_threshold & padj > padj_threshold), 
             aes(color = "Non-significant"), size = 3, alpha = 0.6) +  # Puntos no significativos
  geom_point(data = subset(AD.deseq.results, log2FoldChange >= log2FoldChange_threshold & padj < padj_threshold), 
             aes(color = "padj < 0.001 and Log2FC > 15"), size = 3, alpha = 1) +  # Alta log2FoldChange, baja p-value
  geom_point(data = subset(AD.deseq.results, log2FoldChange <= -log2FoldChange_threshold & padj < padj_threshold), 
             aes(color = "padj < 0.001 and Log2FC < -15"), size = 3, alpha = 1) +  # Baja log2FoldChange, baja p-value
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +  # Línea horizontal para el umbral de p-value
  geom_vline(xintercept = c(-log2FoldChange_threshold, log2FoldChange_threshold), linetype = "dashed", color = "black") +  # Líneas verticales para los umbrales de log2FoldChange
  scale_color_manual(values = c("Non-significant" = "grey", 
                                "padj < 0.001 and Log2FC > 15" = "#DCE319FF", 
                                "padj < 0.001 and Log2FC < -15" = "#238A8DFF")) +  # Definir colores para los puntos
  theme_grey() +  # Estilo del tema
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +  # Etiquetas de los ejes y la leyenda de color
  guides(color = guide_legend(title = NULL)) + # Eliminar el título de la leyenda de color
  theme(axis.text = element_text(size = 12, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los textos de los ejes
        axis.title = element_text(size = 14, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los títulos de los ejes
        legend.text = element_text(size = 10)) +  # Ajustar el tamaño del texto de la leyenda
  ggtitle("S_PP vs S") +  # Añadir título al gráfico
  geom_text_repel(data = subset(AD.signif, abs(log2FoldChange) >= 15 & padj < 0.01), aes(label = Genus), size = 3)  # Añadir etiquetas taxonómicas a cada punto significativo sin solapamiento

# Mostrar el volcano plot
volcano.plot.AD

# Repeating the process for the HFD obesity group. 

obesity.pst.rar = subset_samples(ps.rar, diet=="HFD" | diet=="HFD+PP")
obesity.pst.rar.deseq = phyloseq_to_deseq2(obesity.pst.rar, ~ diet)
# calculating geometric means prior to estimating size factors.
# this normalizes all samples to the same level of reads.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(obesity.pst.rar.deseq), 1, gm_mean)
obesity.pst.rar.deseq = estimateSizeFactors(obesity.pst.rar.deseq, geoMeans = geoMeans)
obesity.pst.rar.deseq = DESeq(obesity.pst.rar.deseq, fitType="local")

HFD.deseq.results = results(obesity.pst.rar.deseq)
HFD.deseq.results = HFD.deseq.results[order(HFD.deseq.results$padj, na.last=NA),]
alpha = 0.01
HFD.signif = HFD.deseq.results[(HFD.deseq.results$padj < alpha), ]
HFD.signif = cbind(as(HFD.signif, "data.frame"), as(tax_table(obesity.pst.rar)[rownames(HFD.signif), ], "matrix"))
head(HFD.signif)
# Enriched ASVs in the pectine treatment group are those with positive logfold change:
positive.HFD.signif = HFD.signif[HFD.signif[, "log2FoldChange"] > 0, ]
positive.HFD.signif = positive.HFD.signif[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
positive.HFD.signif # Lachnospiraceae most enriched Family, log2fold 23.3, related to SCFA production.  

# Set up the plot
library(ggplot2)
library(viridis)

log2FoldChange_threshold <- 5
padj_threshold <- 0.01

# Crear el volcano plot
library(ggrepel)

# Crear el volcano plot
volcano.plot.HFD <- ggplot(HFD.deseq.results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(HFD.deseq.results, abs(log2FoldChange) < log2FoldChange_threshold & padj > padj_threshold), 
             aes(color = "Non-significant"), size = 3, alpha = 0.6) +  # Puntos no significativos
  geom_point(data = subset(HFD.deseq.results, log2FoldChange >= log2FoldChange_threshold & padj < padj_threshold), 
             aes(color = "padj < 0.01 and Log2FC > 5"), size = 3, alpha = 1) +  # Alta log2FoldChange, baja p-value
  geom_point(data = subset(HFD.deseq.results, log2FoldChange <= -log2FoldChange_threshold & padj < padj_threshold), 
             aes(color = "padj < 0.01 and Log2FC < -5"), size = 3, alpha = 1) +  # Baja log2FoldChange, baja p-value
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +  # Línea horizontal para el umbral de p-value
  geom_vline(xintercept = c(-log2FoldChange_threshold, log2FoldChange_threshold), linetype = "dashed", color = "black") +  # Líneas verticales para los umbrales de log2FoldChange
  scale_color_manual(values = c("Non-significant" = "grey", 
                                "padj < 0.01 and Log2FC > 5" = "#DCE319FF", 
                                "padj < 0.01 and Log2FC < -5" = "#238A8DFF")) +  # Definir colores para los puntos
  theme_grey() +  # Estilo del tema
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "Significance") +  # Etiquetas de los ejes y la leyenda de color
  guides(color = guide_legend(title = NULL)) + # Eliminar el título de la leyenda de color
  theme(axis.text = element_text(size = 12, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los textos de los ejes
        axis.title = element_text(size = 14, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los títulos de los ejes
        legend.text = element_text(size = 10)) +  # Ajustar el tamaño del texto de la leyenda
  ggtitle("HFD_PP vs HFD") +  # Añadir título al gráfico
  geom_text_repel(data = subset(HFD.signif, abs(log2FoldChange) >= 5 & padj < 0.01), aes(label = Genus), size = 3)  # Añadir etiquetas taxonómicas a cada punto significativo sin solapamiento

# Mostrar el volcano plot
volcano.plot.HFD


### TAXONOMY AND BETA DIVERSITY

pcoa_jaccard_alzheimers = pcoa_phyloseq(no.batch.pst, "diet", circle = TRUE, method = "bray")
pcoa_jaccard_obesity = pcoa_phyloseq(obesity.pst.rar, "diet", circle = TRUE, method = "jaccard")

pcoa_jaccard_obesity = pcoa_jaccard_obesity +
                        theme_grey() +
                        scale_fill_manual(values=c("#238A8DFF","#DCE319FF","#73D055FF")) +
  theme(axis.text = element_text(size = 12, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los textos de los ejes
        axis.title = element_text(size = 13, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los títulos de los ejes
        legend.text = element_text(size = 10, face = "bold")) +  # Ajustar el tamaño del texto de la leyenda
  ggtitle("PCoA of Jaccard distance matrix")
  
pcoa_jaccard_obesity

pcoa_jaccard_alzheimers = pcoa_jaccard_alzheimers +
  theme_grey() +
  scale_fill_manual(values=c("#238A8DFF","#DCE319FF","#73D055FF")) +
  theme(axis.text = element_text(size = 12, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los textos de los ejes
        axis.title = element_text(size = 13, face = "bold"),  # Ajustar el tamaño y el tipo de fuente de los títulos de los ejes
        legend.text = element_text(size = 10, face = "bold")) +  # Ajustar el tamaño del texto de la leyenda
  ggtitle("PCoA of Jaccard distance matrix")

pcoa_jaccard_alzheimers

# ABUNDANCE HEATMAPS
## what txonomic groups??
alzheimers.pst.rar = subset_samples(pst.rar, diet=="SD" | diet=="S" | diet =="S+PP")
alzheimers.ps.rar = subset_samples(ps.rar, diet=="SD" | diet=="S" | diet =="S+PP")
alzheimers.abundance.heatmap = abundance_heatmap(alzheimers.pst.rar, classification = 'Class',
                                                  treatment = "diet", transformation = 'log2')
obesity.pst.rar = subset_samples(pst.rar, diet =="SD" | diet =="HFD" | diet =="HFD+PP")
obesity.ps.rar = subset_samples(ps.rar, diet =="SD" | diet =="HFD" | diet =="HFD+PP")
obesity.abundance.heatmap = abundance_heatmap(obesity.pst.rar, classification = 'Class',
                                              treatment = "diet", transformation = 'log2')

# phylogeny profile

phylogeny.alzheimers = phylogeny_profile(no.batch.pst, classification = 'Phylum', 
                                         treatment = "diet", merge = TRUE, 
                                         relative_abundance = TRUE)
phylogeny.obesity = phylogeny_profile(obesity.ps.rar, classification = 'Phylum', 
                                      treatment = "diet", merge = TRUE, 
                                      relative_abundance = TRUE)

### LEFSE FOR B IOMARKER DISCOVERY

# installing package to format data for LefSe analysis.
remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)
# Fist we use the treatment.pst.rar which is a subset for the pst.rar with only AD samples. 
# transpose.otus = FALSE as ASVs are already the rownames of the object. 
phyloseq2lefse( ps = treatment.pst.rar,
                covars = "diet",
                file.name = "AD_lefse_data.txt",
                taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
                transpose.otus = FALSE)

## ANCOM-BC 
BiocManager::install("ANCOMBC")
BiocManager::install("mia")
library(ANCOMBC)
library(mia)

test = makeTreeSummarizedExperimentFromPhyloseq(treatment.pst.rar)

## Venn diagram for shared taxa between conditions.

install.packages("eulerr")
library(eulerr)
library(gplots)
library(VennDiagram)
BiocManager::install("microbiome")
library(microbiome)
devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)

# AD venn diagram of shared taxa
AD.venn.data = subset_samples(pst.rar, diet=="SD" | diet=="S" | diet=="S+PP")
table(meta(AD.venn.data)$diet)
AD.pst.relative = transform(AD.venn.data, "compositional")
AD.groups = unique(as.character(meta(AD.pst.relative)$diet))
list_core <- c() # an empty object to store information

for (n in AD.groups){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(AD.pst.relative, diet == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

mycols <- c(nonCRC="#22A884FF", CRC="#440154FF", H="#FDE725FF") 
plot(venn(list_core),
     fills = mycols)

filtered.pst = conglomerate_taxa(pst.rar, "Family")
HFD.PP.pst = subset_samples(pst.rar, diet=="HFD+PP")
filtered.HFD.pst = conglomerate_taxa(HFD.PP.pst,"Family")
S.PP.pst = subset_samples(pst.rar, diet=="S+PP")
filtered.S.pst = conglomerate_taxa(S.PP.pst,"Family")

set.seed(123)
# fast greedy algorithm for clustering and detecting communities based on modularity
co_occurrence_network(filtered.pst, treatment = "diet", subset=c("SD"),classification = "Family", cluster = TRUE, cluster_colors = brewer.pal(4, "Accent"))
co_occurrence_network(filtered.pst, treatment = "diet", subset=c("S"),classification = "Family", cluster = TRUE, cluster_colors = brewer.pal(5, "Accent"))
co_occurrence_network(filtered.S.pst, treatment = "diet",classification = "Family", cluster = TRUE)
co_occurrence_network(filtered.pst, treatment = "diet", subset=c("HFD"),classification = "Family", cluster = TRUE, cluster_colors = brewer.pal(6, "Accent"))
co_occurrence_network(filtered.HFD.pst, treatment = "diet",classification = "Family", cluster = TRUE, cluster_colors = brewer.pal(4, "Accent"))

# the prebiotic reduces in 2 the number of clusters in the co ocurrence networks. 


### BF RATIO 


phyla_obesity <- phyloseq::tax_glom(physeq = obesity.pst.rar, taxrank = "Phylum")
phyla_rel_obesity <- phyloseq::transform_sample_counts(phyla_obesity, function(x) { x/sum(x) } )
tax_table(phyla_rel_obesity)
phyla_rel_bact <- otu_table(subset_taxa(phyla_rel_obesity, Phylum == "Bacteroidota"))
phyla_rel_firm <- suppressWarnings(otu_table(subset_taxa(phyla_rel_obesity, Phylum == "Firmicutes")))
bf_ratio_obesity <- log2(phyla_rel_bact /  phyla_rel_firm)

bf_ratio_obesity # bacteroidetes/firmicutes ratio per sample
bf_ratio_obesity = bf_ratio_obesity[1,]

# Crear un vector con los nombres de las muestras


# Crear un data frame con los datos
data_bf_obesity <- data.frame(
  sampleid = names(bf_ratio_obesity),
  diet = c(rep("HFD", 8), rep("HFD+PP", 5), rep("SD", 8)),
  bf_ratio_obesity = bf_ratio_obesity
)

boxplot_bf_ratio_obesity <- ggplot(data_bf_obesity, aes(x = diet, y = bf_ratio_obesity, fill = diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = expression(bold("(A) B/F Ratio in obese mice")),
       x = "",
       y = "Ratio Bacteroidetes/Firmicutes") +
  scale_fill_manual(values = c("#238A8DFF", "#DCE319FF", "#73D055FF")) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

# Mostrar el boxplot
print(boxplot_bf_ratio_obesity) # Más bajo, mejor. 
TukeyHSD(aov(bf_ratio_obesity ~ diet, data = data_bf_obesity))

# ratio for alzheimers 

phyla_AD <- phyloseq::tax_glom(physeq = alzheimers.pst.rar, taxrank = "Phylum")
phyla_rel_AD <- phyloseq::transform_sample_counts(phyla_AD, function(x) { x/sum(x) } )
tax_table(phyla_rel_AD)
phyla_rel_bact_AD <- otu_table(subset_taxa(phyla_rel_AD, Phylum == "Bacteroidota"))
phyla_rel_firm_AD <- suppressWarnings(otu_table(subset_taxa(phyla_rel_AD, Phylum == "Firmicutes")))
bf_ratio_AD <- log2(phyla_rel_bact_AD /  phyla_rel_firm_AD)

bf_ratio_AD = bf_ratio_AD[1,] # bacteroidetes/firmicutes ratio per sample


# Crear un vector con los nombres de las muestras


# Crear un data frame con los datos
data_bf_AD <- data.frame(
  sampleid = names(bf_ratio_AD),
  diet = c(rep("SD", 8), rep("S", 8), rep("S+PP", 7)),
  bf_ratio_AD = bf_ratio_AD
)

boxplot_bf_ratio_AD <- ggplot(data_bf_AD, aes(x = diet, y = bf_ratio_AD, fill = diet)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = expression(bold("(A) B/F Ratio in AD mice")),
       x = "",
       y = "Ratio Bacteroidetes/Firmicutes") +
  scale_fill_manual(values = c("#238A8DFF", "#DCE319FF", "#73D055FF")) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

# Mostrar el boxplot
print(boxplot_bf_ratio_AD) # Más bajo, mejor. 
TukeyHSD(aov(bf_ratio_AD ~ diet, data = data_bf_AD))



## ANCOM-BC2 DA
library(ANCOMBC)
library(dplyr)
tse = mia::makeTreeSEFromPhyloseq(ps)
# Subsetting the object to the AD experiment
tse = tse[,tse$diet %in% c("SD", "S", "S_PP")]
tse$diet = factor(tse$diet, levels = c("S", "S_PP", "SD")) # S as reference.
# To control FDR from multiple testing Holm-Bonferroni method is applied.
# HB method is robust against inaccuracies in p-values, an issue often seen in DA analysis.
set.seed(123)
out = ancombc(data = tse, assay_name = "counts", 
              tax_level = "Genus", phyloseq = NULL, 
              formula = "diet", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "diet", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)
res = out$res
res_global = out$res_global
tab_lfc = res$lfc
tab_q = res$q
tab_q


df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

library(ggplot2)

# Filtrar los datos para obtener solo los taxones en el grupo dietS+PP
df_lfc_spp <- df_lfc[df_lfc$dietS_PP != 0, ]

# Ordenar los taxones por su LFC en el grupo dietS+PP
df_lfc_spp <- df_lfc_spp[order(df_lfc_spp$dietS_PP, decreasing = TRUE), ]

# Crear el gráfico de barras horizontal
ggplot(df_lfc_spp, aes(x = dietS_PP, y = reorder(taxon_id, dietS_PP))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Taxones diferenciales en el grupo dietS+PP",
       x = "Log Fold Change",
       y = "Taxón") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank())



devtools::install_github("nuriamw/micro4all")


## COMBAT BATCH EFFECT. ESTO ES RELEVANTE. 
library(sva)
data.batch = as(otu_table(pst), "matrix")
data.batch = as.data.frame(data.batch)
pheno = data.frame(
  batch = sample_data(pst)$batch,
  diet = sample_data(pst)$diet
)
row.names(pheno) = row.names(sample_data(pst))
mod = model.matrix(~as.factor(diet), data=pheno)

data.batch = as.matrix(data.batch)
adjusted_count = ComBat_seq(data.batch, pheno$batch, pheno$diet)

OTU_combat<-otu_table(as.data.frame(adjusted_count), taxa_are_rows = T)
taxa_names(OTU_combat)<- taxa_names(phy_taxonomy_norm)
identical(rownames(OTU_combat), rownames(tax))
no.batch.pst = phyloseq(OTU_combat,phy_taxonomy_norm,phy_metadata_norm) #Remove tree_root when working with ITS


library(bapred)
meancentered.otu = meancenter(as.matrix(t(OTU)), pheno$batch)
meancentered.otu = t(meancentered.otu$xadj)

adjusted_count_pca = PCA(adjusted_count)


