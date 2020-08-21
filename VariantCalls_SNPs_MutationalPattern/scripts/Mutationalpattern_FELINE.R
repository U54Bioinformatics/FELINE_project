
# 0 Install packages
if(FALSE){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("MutationalPatterns")
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
  install.packages("gridExtra")
  install.packages("ggplot2")
  install.packages("data.table")
}

prefix='FELINE'
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
# 1 R packages
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg19)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(BSgenome)
library(ref_genome, character.only = TRUE)
library(ggplot2)
library(data.table)
library(gridExtra)

# 2 Mutation files
tab_files <- list.files(path= "./FELINE_WES_variant_tab",
                        pattern = ".tab", full.names = TRUE)

patient_data <- read.table("Patient_1_46.clinical_WES.txt", sep="\t", header=T)

#Sample  Patient.Study.ID        Arm     ARM     Treatment       Ribo    TreatLab        Class   Response        Ki67_Response
#FEL013  001-102 C       C       letrozole + ribo        1       Letrozole + Ribociclib  SD      Non-responder   Acquired resistance
#FEL015  001-105 A2      A       letrozole       0       Letrozole       SD      Non-responder   Innate resistance
sample_names <- patient_data$Sample
sample_arm   <- patient_data$ARM
sample_ribo  <- gsub(" ","", patient_data$Treatment)
sample_response <- patient_data$Response

## create GRange file in a list for each sample (Hoda). vcfs is a list with all the GRanges obj
vcfs <- list()
for (i in 1:length(tab_files)){
  print(i)
  print(tab_files[i])
  Mutation_file <- fread(tab_files[i], sep="\t")
  alfa_single   <- with(Mutation_file, GRanges(CHROM, IRanges(start=Mutation_file$POS, end=Mutation_file$POS), REF= REF, ALT=ALT))
  genome(alfa_single) <- "hg19"
  vcfs[[as.vector(sample_names[i])]]<- (alfa_single)
}
vcfs

# 3 Base substitution types
muts = mutations_from_vcf(vcfs[[1]])
head(muts, 12)
types = mut_type(vcfs[[1]])
head(types, 12)
context = mut_context(vcfs[[1]], ref_genome)
head(context, 12)
type_context = type_context(vcfs[[1]], ref_genome)
lapply(type_context, head, 12)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
type_occurrences

# 4 Mutation spectrum
fig_mutation_spectrum <- paste0(prefix, '.04_mutation_spectrum.pdf')
pdf(fig_mutation_spectrum, width=10, height=7)
## Overall mutation spectrum in all patients
p1 <- plot_spectrum(type_occurrences, CT = TRUE)
## mutation spectrum in each category
p2 <- plot_spectrum(type_occurrences, CT = TRUE, by=sample_arm)
p3 <- plot_spectrum(type_occurrences, CT = TRUE, by=sample_ribo)
p4 <- plot_spectrum(type_occurrences, CT = TRUE, by=sample_response)
## mutation spectrum in each patient
p5 <- plot_spectrum(type_occurrences, CT = TRUE, by=names(vcfs))
#p3 <- plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE, by=names(vcfs))
#grid.arrange(p1, p2, p3, ncol=3, widths=c(3,3,1.75))
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()


# 5 96 mutational profile
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
head(mut_mat)

fig_96_profile_per_patient <- paste0(prefix, '.05_96_profile_per_patient.pdf')
pdf(fig_96_profile_per_patient, width=10, height=7)
for (i in c(1,3,5,7,9,11,12,13)){
  #p1 <- plot_96_profile(mut_mat[,c(i, i+1)])
  p2 <- plot_96_profile(mut_mat[,c(i, i+1)], condensed = TRUE)
  print(p2)
  #print(p2)
}
dev.off()

# 6 De novo mutational signature extraction using NMF
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 10)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B")

fig_96_profile_NMC_denovo <- paste0(prefix, '.06a_96_profile_NMC_denovo.pdf')
pdf(fig_96_profile_NMC_denovo, width=10, height=7)
p1 <- plot(estimate)
p2 <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
p3 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute", coord_flip = TRUE)
p4 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative", coord_flip = TRUE)
p5 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=TRUE)
p6 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=FALSE)
p7 <- plot_compare_profiles(mut_mat[,1],
                            nmf_res$reconstructed[,1],
                            profile_names = c("Original", "Reconstructed"),
                            condensed = TRUE)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 10)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C")

fig_96_profile_NMC_denovo <- paste0(prefix, '.06b_96_profile_NMC_denovo.pdf')
pdf(fig_96_profile_NMC_denovo, width=10, height=7)
p1 <- plot(estimate)
p2 <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
p3 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute", coord_flip = TRUE)
p4 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative", coord_flip = TRUE)
p5 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=TRUE)
p6 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=FALSE)
p7 <- plot_compare_profiles(mut_mat[,1],
                            nmf_res$reconstructed[,1],
                            profile_names = c("Original", "Reconstructed"),
                            condensed = TRUE)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

# 7 Find optimal contribution of known signatures
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
## Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
## Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
## Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
## Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

fig_96_profile_compare_COSMIC <- paste0(prefix, '.07_96_profile_compare_COSMIC.pdf')
pdf(fig_96_profile_compare_COSMIC, width=10, height=7)
p1 <- plot_96_profile(cancer_signatures[,1:2], condensed = TRUE, ymax = 0.3)
## Hierarchically cluster the COSMIC signatures based on their similarity with average linkage
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
p2 <- plot(hclust_cosmic)
## Similarity between mutational profiles and COSMIC signatures
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
p3 <- plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order, cluster_rows = TRUE)
## Find optimal contribution of COSMIC signatures to reconstruct 96 mutational profiles
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
### Select signatures with some contribution
select <- which(rowSums(fit_res$contribution) > 10)
p4 <- plot_contribution(fit_res$contribution[select,],
                        cancer_signatures[,select],
                        coord_flip = TRUE,
                        mode = "relative")
p5 <- plot_contribution_heatmap(fit_res$contribution, cluster_samples = TRUE, method = "complete")
## Compare the reconstructed mutational profile of sample 1 with its original mutational profile
p6 <- plot_compare_profiles(mut_mat[,1], fit_res$reconstructed[,1],
                            profile_names = c("Original", "Reconstructed"),
                            condensed = TRUE)
## Calculate the cosine similarity between all original and reconstructed mutational profiles
cos_sim_ori_rec <- cos_sim_matrix(mut_mat, fit_res$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
p7 <- ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
  geom_bar(stat="identity", fill = "skyblue4") +
  coord_cartesian(ylim=c(0.8, 1)) +
  ylab("Cosine similarity\n original VS reconstructed") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
  geom_hline(aes(yintercept=.95))

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()


# 8 Transcriptional strand bias analysis
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes_hg19
tissue <- as.vector(sample_names)
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
mut_mat_s[1:5,1:5]
strand_counts <- strand_occurrences(mut_mat_s, by=tissue)
head(strand_counts)
strand_bias <- strand_bias_test(strand_counts)
strand_bias

fig_strand_bias_transcriptional <- paste0(prefix, '.08_strand_bias_transcriptional.pdf')
pdf(fig_strand_bias_transcriptional, width=10, height=7)
p1 <- plot_strand(strand_counts, mode = "relative")
p2 <- plot_strand_bias(strand_bias)
grid.arrange(p1, p2)
dev.off()

# 9 Replicative strand bias analysis
repli_file = system.file("extdata/ReplicationDirectionRegions.bed",
                         package = "MutationalPatterns")
repli_strand = read.table(repli_file, header = TRUE)
repli_strand_granges = GRanges(seqnames = repli_strand$Chr,
                               ranges = IRanges(start = repli_strand$Start + 1,
                                                end = repli_strand$Stop),
                               strand_info = repli_strand$Class)
seqlevelsStyle(repli_strand_granges) = "UCSC"
repli_strand_granges
strand_rep <- mut_strand(vcfs[[1]], repli_strand_granges, mode = "replication")
mut_mat_s_rep <- mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges,
                                     mode = "replication")
mut_mat_s_rep[1:5, 1:5]
repli_strand_granges$strand_info <- factor(repli_strand_granges$strand_info,
                                           levels = c("right", "left"))
mut_mat_s_rep2 <- mut_matrix_stranded(vcfs, ref_genome, repli_strand_granges,
                                      mode = "replication")
mut_mat_s_rep2[1:5, 1:5]
strand_counts_rep <- strand_occurrences(mut_mat_s_rep, by=tissue)
head(strand_counts)
strand_bias_rep <- strand_bias_test(strand_counts_rep)

fig_strand_bias_replicative <- paste0(prefix, '.09_strand_bias_replicative.pdf')
pdf(fig_strand_bias_replicative, width=10, height=7)
p1 <- plot_strand(strand_counts_rep, mode = "relative")
p2 <- plot_strand_bias(strand_bias_rep)
grid.arrange(p1, p2)
dev.off()

# 10 Extract signatures with strand bias
nmf_res_strand <- extract_signatures(mut_mat_s, rank = 2)
colnames(nmf_res_strand$signatures) <- c("Signature A", "Signature B", "Signature C")

fig_signature_with_strand_bias <- paste0(prefix, '.10_signature_with_strand_bias.pdf')
pdf(fig_signature_with_strand_bias, width=10, height=7)
p1 <- plot_192_profile(nmf_res_strand$signatures, condensed = TRUE)
p2 <- plot_signature_strand_bias(nmf_res_strand$signatures)
grid.arrange(p1, p2, ncol = 2, widths = c(5, 1.8))
dev.off()

# 11 Genomic distribution: A rainfall plot visualizes mutation types and intermutation distance.
chromosomes <- seqnames(get(ref_genome))[1:22]
fig_genomic_distr_rainfall <- paste0(prefix, '.11_genomic_distr_rainfall.pdf')
pdf(fig_genomic_distr_rainfall, width=10, height=7)
p1 <- plot_rainfall(vcfs[[1]], title = names(vcfs[1]),
                    chromosomes = chromosomes, cex = 1.5, ylim = 1e+09)
print(p1)
dev.off()
