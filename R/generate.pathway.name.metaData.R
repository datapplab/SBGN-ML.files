all.pathway.info = read.table("./data/SBGN/pathwayCommons.all.pathways.info.no.kegg.tsv"
                              ,sep="\t"
                              ,quote=""
                              ,header=T
                              )
head(all.pathway.info)
dim(all.pathway.info)

pre.generated.files = read.table("./data/SBGN/pathwayCommons.file.names.txt"
                              ,sep="\t"
                              ,quote=""
                              ,header=F
                              )
colnames(pre.generated.files) = "Name_of_pre.generated_SBGN_file"
head(pre.generated.files)
dim(pre.generated.files)

pathway.info.pre.generated.files = merge(pre.generated.files,all.pathway.info,all.x = T )
pathway.info.pre.generated.files = pathway.info.pre.generated.files[!is.na(pathway.info.pre.generated.files[,2]),]
head(pathway.info.pre.generated.files)
dim(pathway.info.pre.generated.files)
save(pathway.info.pre.generated.files,file="./data/SBGN/pathway.info.pre.generated.files.RData")
table(pathway.info.pre.generated.files[,"sub.database"])
write.table(pathway.info.pre.generated.files,"./data/SBGN/pathway.info.pre.generated.files.tsv",row.names=F,sep="\t",quote=FALSE)


system("ls ./data/SBGN/MetaCyc >./data/SBGN/MetaCyc.SBGN.file.names.tsv ")
system("ls ./data/SBGN/MetaCrop >./data/SBGN/MetaCrop.SBGN.file.names.tsv ")
# Manualy add MetaCyc and MetaCrop information in the file ./data/SBGN/pathway.info.pre.generated.files.tsv
pathway.info.pre.generated.files = read.table("./data/SBGN/pathway.info.pre.generated.files.tsv"
                              ,sep="\t"
                              ,quote=""
                              ,header=T
                              )

head(pathway.info.pre.generated.files)
dim(pathway.info.pre.generated.files)
save(pathway.info.pre.generated.files,file="./data/SBGN/pathway.info.pre.generated.files.RData")
table(pathway.info.pre.generated.files[,"sub.database"])

# add pathway description for metacyc
metacyc.id.to.name = read.table("C:/xiaoxi/work/input/metaCyc/metaCyc/pathways.col.id-to-description.tsv.txt",stringsAsFactors=F,header=T,sep="\t",quote='"')
head(metacyc.id.to.name)
row.names(metacyc.id.to.name) = metacyc.id.to.name[,1]
head(metacyc.id.to.name)
dim(metacyc.id.to.name)
library(SBGNview)
data(pathway.info.pre.generated.files)
head(pathway.info.pre.generated.files)

dim(pathway.info.pre.generated.files[pathway.info.pre.generated.files[,"pathway.database"] == "MetaCyc",])
head(pathway.info.pre.generated.files[pathway.info.pre.generated.files[,"pathway.database"] == "MetaCyc",])
metacyc.ids = pathway.info.pre.generated.files[pathway.info.pre.generated.files[,"pathway.database"] == "MetaCyc","pathway.id"]
head(metacyc.ids)
metacyc.names = metacyc.id.to.name[metacyc.ids,2]
head(metacyc.names)
length(metacyc.names)
pathway.info.pre.generated.files[pathway.info.pre.generated.files[,"pathway.database"] == "MetaCyc","DISPLAY_NAME"] = metacyc.names
head(pathway.info.pre.generated.files[pathway.info.pre.generated.files[,"pathway.database"] == "MetaCyc",])
save(pathway.info.pre.generated.files,file = "../SBGNview/data/pathway.info.pre.generated.files.RData")
head(pathway.info.pre.generated.files)
tail(pathway.info.pre.generated.files)
dim(pathway.info.pre.generated.files)












