source("http://www.bioconductor.org/biocLite.R")
# the below command is used to install the bioconductor package
# biocLite("SRAdb")

library(SRAdb)

ids = read.table("HMASM_test_subset.csv", header=T)

sqlfile <- getSRAdbFile()
sra_connection <- dbConnect(SQLite(), "SRAmetadb.sqlite")
### NOTE: a DB file of ~22 GB will be downloaded, you'll probably want to delete this afterwards

conversion_table = sraConvert( ids$SRS , out_type=c("run"), sra_con = sra_connection)

write.table(conversion_table , file = "SRS_SRR_ids_linked.txt" , quote=FALSE , sep="\t" , row.names=F)
