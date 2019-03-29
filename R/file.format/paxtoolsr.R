# install.packages("rJava")
# need to install JRE and JDK
#  sudo apt-get install default-jdk
# then register jdk with R
#  sudo R CMD javareconf
# source("https://bioconductor.org/biocLite.R")
# biocLite("paxtoolsr")
# options(java.parameters = "-Xmx10g")
# Sys.setenv(JAVA_HOME="c:/Program Files (x86)/Java/jre1.8.0_191")
library(rJava)
library(paxtoolsr)
options(java.parameters = "-Xmx10240m")
# input.file = "/home/xiaoxi/0work/input/metacyc/biopax/pathway-biopax"
# input.folder = "/home/xiaoxi/0work/input/metacyc/biopax.level3"
# input.folder = "/media/sf_work/input/2018.10.3.linux/metacyc/biopax.level3/bile.acid"
# input.folder = "/media/sf_work/input/2018.10.3.linux/metacyc/biopax.level3/bile.acid"
input.folder = "/media/sf_work/input/plant.reactome/downloaded"
# input.folder = "/home/xiaoxi/0work/input/innatedb"
# input.file = "/media/sf_work/input/metaCyc/biopax-level3.owl"
input.files = list.files(path=input.folder,full.names=F,pattern="biopax")
head(input.files)

# output.folder = "/home/xiaoxi/0work/input/innatedb"
# output.folder = "/home/xiaoxi/0work/input/metacyc/biopax.test.sbgn"
output.folder = "/media/sf_work/input/plant.reactome/downloaded.sbgn"
system(paste("mkdir ",output.folder,sep=""))
for(i in 1:length(input.files)){
        input.file = input.files[i]
        input.file
        
        output.file.name = input.file
        output.file = gsub(" ","",output.file.name)
        output.file.name
        output.file = paste(output.folder,"/",output.file,".sbgn",sep="")
        output.file
        
        input.file = paste(input.folder,input.file,sep="/")
        input.file
        input.file.size = file.info(input.file)$size
        input.file.size
        cat("\n\n\n")
        print(input.file)
        
        input.content = readChar(input.file,nchars = input.file.size)
        substr(input.content,300,400)
        input.content = gsub("xml:base=\"[.]+\"",paste("xml:base=\"",output.file.name,"\"",sep=""),input.content)
        substr(input.content,300,400)
        
        if(file.exists(output.file) | input.file.size == 0){
                return()
        }
        print(input.file)
        print(output.file)
        toSBGN(input.file,output.file)
}
