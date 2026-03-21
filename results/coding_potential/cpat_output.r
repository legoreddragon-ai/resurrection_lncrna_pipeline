load("references/Plant_logitModel.RData")
test <- read.table(file="results/coding_potential/cpat_output.ORF_info.tsv",sep="\t",header=T)
test$Coding_prob <- predict(mylogit,newdata=test,type="response")
write.table(test, file="results/coding_potential/cpat_output.ORF_prob.tsv", quote=F, sep="\t",row.names=FALSE, col.names=TRUE)
