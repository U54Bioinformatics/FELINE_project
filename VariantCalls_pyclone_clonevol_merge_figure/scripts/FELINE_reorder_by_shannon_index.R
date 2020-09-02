library(tidyverse)

x = read.table("FELINE_patient_info.sorted.revised.txt", header=T, sep="\t")
y = read.table("FELINE_shannon_index.txt", header=T, sep="\t")
z = read.table("FELINE_dominnace_index.txt", header=T, sep="\t")
z$deltaDominance = z$Dominance - z$Dominance0

y %>% filter(!Timepoint %in% c(0)) %>% select(Patient, H0, H, delta_H) -> y
names(y) <- c("orig.ident", "H0", "H", "delta_H")

z %>% filter(!Timepoint %in% c(0)) %>% select(Patient, Dominance, Dominance0, deltaDominance) -> z
names(z) <- c("orig.ident", "Dominance", "Dominance0", "deltaDominance")

merge(x,y, on="orig.ident", all=TRUE) -> xy
merge(xy, z, on="orig.ident", all=TRUE) -> xyz
xyz %>% filter(!orig.ident %in% c("FEL016", "FEL021", "FEL029", "FEL032", "FEL040")) -> data

data$delta_H[is.na(data$delta_H)] = -10
data$deltaDominance[is.na(data$deltaDominance)] = -10
data$H0[is.na(data$H0)] = 0
data$Dominance0[is.na(data$Dominance0)] = 0

# rank by arm then shannon_index_delta
data %>% arrange(ARM, desc(delta_H), orig.ident) -> data_1
write.table(data_1, "FELINE_patient_info.sorted_by_Arm_shannon_delta.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)

# rank by arm then shannon_index_H0
data %>% arrange(ARM, desc(H0), orig.ident) -> data_2
write.table(data_2, "FELINE_patient_info.sorted_by_Arm_shannon_H0.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)

# rank by shannon_index_delta
data %>% arrange(desc(delta_H), orig.ident) -> data_3
write.table(data_3, "FELINE_patient_info.sorted_by_shannon_delta.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)

# rank by shannon_index_H0
data %>% arrange(desc(H0), orig.ident) -> data_4
write.table(data_4, "FELINE_patient_info.sorted_by_shannon_H0.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)


# rank by arm then deltaDominance
data %>% arrange(ARM, desc(deltaDominance), orig.ident) -> data_5
write.table(data_5, "FELINE_patient_info.sorted_by_Arm_deltaDominance.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)

# rank by arm then Dominance0
data %>% arrange(ARM, desc(Dominance0), orig.ident) -> data_6
write.table(data_6, "FELINE_patient_info.sorted_by_Arm_Dominance0.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)

# rank by deltaDominance
data %>% arrange(desc(deltaDominance), orig.ident) -> data_7
write.table(data_7, "FELINE_patient_info.sorted_by_deltaDominance.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)

# rank by Dominance0
data %>% arrange(desc(Dominance0), orig.ident) -> data_8
write.table(data_8, "FELINE_patient_info.sorted_by_Dominance0.revised.txt", sep="\t", quote=F, row.names=F, col.names=T)



