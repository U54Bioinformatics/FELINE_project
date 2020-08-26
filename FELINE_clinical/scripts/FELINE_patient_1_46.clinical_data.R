x = read.table("Patient_1_46.clinical_data.metadata.txt.Ki67_Response.txt", sep="\t", header=1)
#table(x$ARM, x$Response_v3b)
#x_table <- table(x$ARM, x$Response_v3b)
x$Response[with(x, Response_v3 %in% c("Progressive disease", "Stable disease", "Rebound disease"))] <- "Non-responder"
x$Response[with(x, Response_v3 %in% c("Sustained response", "Partial response"))] <- "Responder"
table(x$Response, x$Response_v3)
table(x$Response, x$Response_v3b)
table(x$Response_v3, x$Response_v3b)
write.table(x, "FELINE_patient_1_46.clinical_data.txt", quote=F, row.names=F, sep="\t")

