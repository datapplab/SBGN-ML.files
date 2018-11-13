
download.from.kegg = function(ids,url,output.file,size.slice){
        all.species.entrez = paste(ids,collapse = "+")
        kegg.q.all.species.entrez.to.all.species.uniprot = paste(url,all.species.entrez,sep="")
        substr(kegg.q.all.species.entrez.to.all.species.uniprot ,1,50)
        # print(kegg.q.all.species.entrez.to.all.species.uniprot)
        # stop()
        if.found = try(download.file(kegg.q.all.species.entrez.to.all.species.uniprot ,destfile = output.file,quiet = T,mode = "a"))
        if(class(if.found) == "try-error"){
                n.genes = length(ids)
                slices = seq(1,n.genes,by=size.slice)
                for(j in 1:length(slices)){
                        if(j == length(slices)){
                                end.index = n.genes
                        }else{
                                end.index = slices[j+1] -1
                        }
                        slice.genes = slices[j]:end.index
                        all.species.entrez = paste(ids[slice.genes],collapse = "+")
                        cat("downloading to ", output.file,"\n",paste(min(slice.genes),max(slice.genes),sep="-"),"out of",n.genes,"genes\n")
                        kegg.q.all.species.entrez.to.all.species.uniprot = paste(url,all.species.entrez,sep="")
                        substr(kegg.q.all.species.entrez.to.all.species.uniprot ,1,50)
                        if.downloaded = try(download.file(kegg.q.all.species.entrez.to.all.species.uniprot ,destfile = output.file,quiet = T,mode = "a"))
                        if(class(if.downloaded) == "try-error"){
                                download.from.kegg(ids[slice.genes],url,output.file,(size.slice-1))
                        }
                }
        }
}


biopax.folder = "C:/xiaoxi/work/papers/SBGNview/package/data.package/pathway.commons/raw.data/v9"
protein.mapping =  readRDS(paste(biopax.folder,"/proteins.mapping.list.rds",sep=""))
names(protein.mapping)
protein.mapping$duplicated.ids
nrow(protein.mapping$mappings.no.database.id)
head(protein.mapping$mappings.database.id.name)
pathwayCommons_UNIPROT_human = protein.mapping$mappings.database.id.name[,2:1]
colnames(pathwayCommons_UNIPROT_human) = c("UNIPROT","pathwayCommons")
head(pathwayCommons_UNIPROT_human)
head(pathwayCommons_UNIPROT_human)
dim(pathwayCommons_UNIPROT_human)
human.uniprot.to.entrez.file = "kegg.human.uniprot.to.entrez.tsv"
human.entrez.to.ko.file = "kegg.human.entrez.to.ko.tsv"
all.species.entrez.to.uniprot.file = "kegg.ko.to.sll.species.entrez.tsv"
kegg.humanUniprot.to.NonHumanUniprot.orthologs.file = "C:/xiaoxi/work/input/id.mapping/kegg/kegg.humanUniprot.to.NonHumanUniprot.orthologs.tsv"
# file.remove(human.uniprot.to.entrez.file,human.entrez.to.ko.file,all.species.entrez.to.uniprot.file,kegg.humanUniprot.to.NonHumanUniprot.orthologs.file)
all.uniprots = unique(pathwayCommons_UNIPROT_human[,"UNIPROT"])
n.pathwayCommons.uniprot = length(all.uniprots)


# human uniprot to genes
# current.uniprot = all.uniprots[1:500]
current.uniprot = all.uniprots
output.folder = "C:/xiaoxi/work/input/id.mapping/kegg"
output.file = paste(output.folder,"/",human.uniprot.to.entrez.file,sep="")
download.from.kegg (paste("up:",current.uniprot,sep=""),"http://rest.kegg.jp/conv/hsa/",output.file,size.slice = 200)
kegg.q.uniprot.to.entrez.human = try(read.table(output.file,header = F,stringsAsFactors = F,sep="\t",quote=""))
colnames(kegg.q.uniprot.to.entrez.human) =c("uniprot","entrez")
head(kegg.q.uniprot.to.entrez.human)

# entrez to ko
output.file = paste(output.folder,"/",human.entrez.to.ko.file,sep="")
download.from.kegg (unique(kegg.q.uniprot.to.entrez.human[,"entrez"]),"http://rest.kegg.jp/link/ko/",output.file,size.slice = 200)
# download.file(kegg.q.human.entrez.to.ko,destfile = paste(output.folder,"/",human.entrez.to.ko.file,sep=""),mode="a")
kegg.human.entrez.to.ko = try(read.table(paste(output.folder,"/",human.entrez.to.ko.file,sep=""),header = F,stringsAsFactors = F,sep="\t",quote=""))
# if(!file.exists(paste(output.folder,"/",human.entrez.to.ko.file,sep=""))){
colnames(kegg.human.entrez.to.ko) =c("entrez.human","ko")
head(kegg.human.entrez.to.ko)

output.file = paste(output.folder,"/",all.species.entrez.to.uniprot.file,sep="")
download.from.kegg (unique(kegg.human.entrez.to.ko[,"ko"]),"http://rest.kegg.jp/link/genes/",output.file,size.slice = 20)
kegg.ko.to.all.species.entrez = try(read.table(output.file,header = F,stringsAsFactors = F,sep="\t",quote=""))
colnames(kegg.ko.to.all.species.entrez ) =c("ko","all.species.entrez")
head(kegg.ko.to.all.species.entrez )
dim(kegg.ko.to.all.species.entrez )
table(kegg.ko.to.all.species.entrez[,"ko"] )

# all species entrez to all species uniprot
dim(kegg.ko.to.all.species.entrez)
by(
        kegg.ko.to.all.species.entrez
        ,as.factor(kegg.ko.to.all.species.entrez[,"ko"])
        ,function(all.species.entrez.df){
                ko = all.species.entrez.df[1,"ko"]
                ko = gsub("ko:","",ko)
                ko.uniprot.file =  paste("C:/xiaoxi/work/input/id.mapping/kegg/ko/",ko,".tsv",sep="")
                
                output.file = ko.uniprot.file
                if(file.exists(output.file)){
                        print(ko)
                        return()
                }
                download.from.kegg (unique(all.species.entrez.df[,"all.species.entrez"]),"http://rest.kegg.jp/conv/uniprot/",output.file,size.slice = 200)
                
                all.species.entrez.to.uniprot = read.table(ko.uniprot.file,header = F,stringsAsFactors = F,sep="\t",quote="")
                colnames(all.species.entrez.to.uniprot) =c("gene","uniprot")
                all.species.entrez.to.uniprot[,"uniprot"] = gsub("up:","",all.species.entrez.to.uniprot[,"uniprot"])
                species = do.call(rbind,strsplit(all.species.entrez.to.uniprot[,"gene"],":"))[,1]
                uniprot.to.ko = all.species.entrez.to.uniprot[,c("uniprot","gene")]
                uniprot.to.ko = cbind(uniprot.to.ko,ko)
                uniprot.to.ko.file = "C:/xiaoxi/work/input/id.mapping/kegg/uniprot.to.ko.tsv"
                write.table(uniprot.to.ko,uniprot.to.ko.file,row.names = F,col.names =F,sep="\t",quote=F,append = T)
                # stop()
        }
)