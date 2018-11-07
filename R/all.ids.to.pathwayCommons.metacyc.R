
#' This is modified version of function geneannot.map() in pathview package. The difference is: 1.Uses all packages(pods) instead of one package for one species. 2. Always generate all-to-all mapping table(unique.map=F)
#' modified on 10/03/2018 previously the input ids are only human ids. In this version, the function can handle all species in bods and only map "eg" to target id

geneannot.map.all <- function( in.type, out.type, org="Hs", pkg.name=NULL, unique.map=TRUE, na.rm=TRUE, keep.order=F){
  library(pathview)
  if(is.null(pkg.name)) {#pkg.name=paste("org", org, "eg.db", sep=".")
    data(bods)
    ridx=grep(tolower(paste0(org, "[.]")), tolower(bods[,1]))
    if(length(ridx)==0) {
      ridx=grep(tolower(org), tolower(bods[,2:3])) %% nrow(bods)
      if(length(ridx)==0) stop("Wrong org value!")
      if(any(ridx==0)) ridx[ridx==0]=nrow(bods)
    }
    pkg.name=bods[ridx,1]
  }
  all.mappings = character(3)
  print("printing nrwos")
  print(nrow(bods))
  for (i in 1:nrow(bods)){
      cat("\n\n\n\n\n\n",i," \n\n")
      print(bods[i,"species"])
      pkg.name = bods[i,1]
      # generate input gene id list(only "eg") for this species
      if( bods[i,"id.type"] != "eg"){
          next()
      }
      gene.data <- try(sim.mol.data(mol.type = "gene", species = bods[i,"kegg code"],id.type = in.type, nmol = 99999999999,nexp=2))
      # gene.data <- try(sim.mol.data(mol.type = "gene", species = "eco",id.type = "ENSEMBL", nmol = 99999999999,nexp=2))
      if(class(gene.data) == "try-error"){
              next()
      }
      
      in.ids = row.names(gene.data)
      pkg.on=try(requireNamespace(pkg.name),silent = TRUE)
      if(!pkg.on) {
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg.name, suppressUpdates =TRUE)
        pkg.on=try(requireNamespace(pkg.name),silent = TRUE)
        if(!pkg.on) stop(paste("Fail to install/load gene annotation package ", pkg.name, "!",  sep=""))
      }
      if(! pkg.name   %in% rownames(installed.packages())){
        cat("\n\n\n",pkg.name," not installed!\n\n")
        stop()
      }
      print(pkg.name)
      db.obj <- eval(parse(text=paste0(pkg.name, "::", pkg.name)))
      id.types <- columns(db.obj) #columns(eval(as.name(pkg.name)))
      in.type=toupper(in.type)
      out.type=toupper(out.type)
      eii=in.type==toupper("entrez") | in.type==toupper("eg")
      if(any(eii)) in.type[eii]="ENTREZID"
      eio=out.type==toupper("entrez") | out.type==toupper("eg")
      if(any(eio)) out.type[eio]="ENTREZID"
      if(in.type==out.type) stop("in.type and out.type are the same, no need to map!")
      nin=length(in.type)
      if(nin!=1) stop("in.type must be of length 1!")
      out.type=out.type[!out.type %in% in.type]
      nout=length(out.type)
      msg <-  paste0("must from: ", paste(id.types, collapse=", "), "!")
      if (! in.type %in% id.types) stop("'in.type' ", msg)
      if (! all(out.type %in% id.types)){
         cat("\n\n\n out.type ", msg,"\n\n\n")
         next()
      }
      in.ids0=in.ids
      in.ids <- unique(as.character(in.ids))#unique necessary for select()# if(unique.map)
      out.ids=character(length(in.ids))
      cat("\n\n querying database \n\n")
      print(in.type)
      print(out.type)
      res <- try(suppressWarnings(select(db.obj,
                                         keys = in.ids,
                                         keytype = in.type,
                                         columns=c(in.type, out.type))))
      # cat("\n done database \n\n")
      res = cbind(res,bods[i,"species"])
      all.mappings = rbind(all.mappings,res)
      # print(head(res))
      # print(dim(res))
      # print(dim(all.mappings))
      # print(head(all.mappings))
      # print(table(all.mappings[,3]))
  }
  print("finished all species")
  # print(head(all.mappings))
  # print(table(all.mappings[,3]))
  # print(head(all.mappings))
  # print(dim(all.mappings))
  # stop()
  res = as.data.frame(all.mappings[-1,])
  colnames(res)[3]="species.pathview"
  # if(class(res)=="data.frame"){
  # print(in.type)
  # print(out.type)
  # print(head(res))
  res <- res[, c(in.type, out.type,"species.pathview")]
  # if(nout==1) na.idx <- is.na(res[,2])
  na.idx <- is.na(res[,2])
  # else na.idx <- apply(res[,-1],1,function(x) all(is.na(x)))
  if (sum(na.idx)>0) {
    n.na <- length(unique(res[na.idx, 1]))
    if (n.na>0) {
      print(paste("Note:", n.na, "of", length(in.ids), "unique input IDs unmapped."))
    }
    if (na.rm) res <- res[!na.idx, ]
  }
  cns=colnames(res)
  res=as.matrix(res)
  rownames(res)=NULL
  return(res)
}

#' generate pairwise mapping tables from all id types to SBGN ids

pairwiseMapping.for.sbgnIds = function(input.type,output.type,cpd.or.gene,id.to.chebi=c(),SBGNview.uniprot.chebi.to.SBGN.id.mapping){
  if(output.type == "pathwayCommons"){  # if output is pathwayCommons id, we can't use pathview.
    if(cpd.or.gene == "compound"){  # we need to convert other ids to chebi, then convert chebi to pathwaycommons
        chebi.to.pathway.commons.mapping = SBGNview.uniprot.chebi.to.SBGN.id.mapping$compound
        chebi.to.pathway.commons.mapping[,"ChEBI"] = gsub("CHEBI:","",chebi.to.pathway.commons.mapping[,"ChEBI"] )
      if(input.type == "chebi"){
        mapping = chebi.to.pathway.commons.mapping
      }else{
        # print(head(chebi.to.pathway.commons.mapping))
        # print(head(id.to.chebi))
        print("111111")
        mapping = merge(id.to.chebi,chebi.to.pathway.commons.mapping,by.x="chebi",by.y="ChEBI",all=F)
        # print(head(mapping))
        # stop()
        mapping = mapping[,c(input.type,"pathwayCommons.id")]
        # mapping = mapping[!is.na(mapping[,2]),]
        colnames(mapping)[2] = "pathwayCommons"
      }
    }else if(cpd.or.gene == "gene"){
          # data("SBGNview.uniprot.chebi.to.SBGN.id.mapping")
      uniprot.to.pathway.commons.mapping = SBGNview.uniprot.chebi.to.SBGN.id.mapping$gene
      if(input.type == "UNIPROT"){
        mapping = uniprot.to.pathway.commons.mapping
      }else{
        id.map <- geneannot.map.all(in.type = input.type, out.type = "UNIPROT")
        print("merging mappings")
        # print(head(id.map))
        # print(head(uniprot.to.pathway.commons.mapping))
        # mapping = merge(id.map,uniprot.to.pathway.commons.mapping[,1:2],by.x="UNIPROT",by.y="UNIPROT",all.y=T)
        mapping = merge(id.map,uniprot.to.pathway.commons.mapping,by.x="UNIPROT",by.y="UNIPROT",all.y=T)
          print("printing species count pathwayCommons")
        #   print(head(mapping))
        #   print(input.type)
        # print(table(id.map[,3]))
        # print(table(mapping[,c(3,5)]))
        # print(head(mapping))
        # stop()
        mapping = mapping[,c(input.type,"pathwayCommons","species")]
        colnames(mapping)[2] = "pathwayCommons"
      }
    }
  }else if(output.type == "metacyc.SBGN"){
    if(cpd.or.gene == "compound"){  # we need to convert other ids to chebi, then convert chebi to pathwaycommons
        if(input.type == "chebi"){
          mapping = compound.to.chebi.mapping
        }else{
          # print(head(chebi.to.pathway.commons.mapping))
          # print(head(id.to.chebi))
          print("2222222")
          mapping = merge(id.to.chebi,compound.to.chebi.mapping,by.x="chebi",by.y="chebi",all=F)
          # print(head(mapping))
          # stop()
          mapping = mapping[,c(input.type,"metacyc.SBGN")]
          # mapping = mapping[!is.na(mapping[,2]),]
          colnames(mapping)[2] = "metacyc.SBGN"
        }
    }else if(cpd.or.gene == "gene"){
      # load("/media/sf_work/papers/SBGNview/package/data.package/metacyc/merged.tables/uniprot.pathway.sbgnBiopaxProtein.RData")
      load("C:/xiaoxi/work/papers/SBGNview/package/data.package/metacyc/merged.tables/uniprot.pathway.sbgnBiopaxProtein.RData")
      uniprot.to.metacyc.mapping = as.matrix(uniprot.pathway.sbgnBiopaxProtein)
      uniprot.to.metacyc.mapping = uniprot.to.metacyc.mapping[!is.na(uniprot.to.metacyc.mapping[,"Uniprot"]),]
    # table(apply(uniprot.to.metacyc.mapping,1,function(line){sum(is.na(line))}))
      if(input.type == "UNIPROT"){
      # if(F){
          # mapping = uniprot.to.metacyc.mapping
          # mapping = mapping[!is.na(mapping[,"sbgnBiopaxProtein"]),]
          id.map <- geneannot.map.all(in.type = "ACCNUM", out.type = "UNIPROT")
          id.map = id.map[,c("UNIPROT","species.pathview")]
          # mapping = merge(id.map,uniprot.to.metacyc.mapping,by.x="UNIPROT",by.y="Uniprot",all=F)
          
      }else{
          # print(head(id.map))
          # print(head(uniprot.to.metacyc.mapping))
          id.map <- geneannot.map.all(in.type = input.type, out.type = "UNIPROT")
      }
          mapping = merge(id.map,uniprot.to.metacyc.mapping,by.x="UNIPROT",by.y="Uniprot",all=F)
          mapping = mapping[!is.na(mapping[,"sbgnBiopaxProtein"]) & !is.na(mapping[,"pathway"]),]
          # print(head(mapping))
          # stop()
          print("printing species count metacyc")
          print(table(id.map[,ncol(id.map)]))
          print(table(mapping[,"species.pathview"]))
          print(head(mapping))
          mapping = mapping[,c(input.type,"pathway","sbgnBiopaxProtein","species.pathview")]
          print ("generating new id with pathway names")
          # uniq.sbgn.id = paste(mapping[,"sbgnBiopaxProtein"],mapping[,"pathway"],sep=":@:")
          uniq.sbgn.id = as.character(mapping[,"sbgnBiopaxProtein"])
          cat("\n\n\n\n\n retrieving species \n\n\n\n")
          print(head(mapping))
          print(head(uniq.sbgn.id))
          mapping = cbind(as.character(mapping[,1]),as.character(uniq.sbgn.id),as.character(mapping[,"species.pathview"]))
          load("C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/id.mapping/mapping.each.pair/for.use.all/species.mapping.RData")
          mapping[,3] = species.mapping[as.character(mapping[,3]),1]
          print("cbind species done")
          # print(head(mapping))
          cat("\n\n\n\n\n")
          colnames(mapping) = c(input.type,"metacyc.SBGN","species")
    }
  }else {
    if(cpd.or.gene == "compound"){
      cat("\n\n\n using pathview to change compound id\n\n")
      if(output.type == "ChEBI"){  # if output is chebi, we can't use pathview.  need to convert other types to kegg, then convert kegg to chebi
        print("333333")
        in.data.target.id = id.to.chebi(data.input.id,input.type,fun.for.multi.id.node)
      }else if(output.type == "pathwayCommons.cpd.name"){ # output is cpd name, use chebi to cpd.name table to convert
        if(input.type == "ChEBI"){
          data.chebi.id = data.input.id
        }else{
          print("44444")
          data.chebi.id = id.to.chebi(data.input.id,input.type,fun.for.multi.id.node)
        }
        data("cpd.name.pathwayCommons.to.chebi")
        in.data.target.id <- mol.sum.multiple.mapping(mol.data = data.input.id, id.map = cpd.name.pathwayCommons.to.chebi,sum.method = fun.for.multi.id.node,input.sbgn.file = input.sbgn.file)
      }else{
        id.map <- cpdidmap(in.ids = row.names(data.input.id), in.type = input.type , out.type = output.type)
        in.data.target.id <- mol.sum.multiple.mapping(mol.data = data.input.id, id.map = id.map,sum.method = fun.for.multi.id.node,input.sbgn.file = input.sbgn.file)
      }
    }else if(cpd.or.gene == "gene"){
      cat("\n\n\n using pathview to change gene id\n\n")
      mapping <- geneannot.map.all(in.type = input.type, out.type = output.type)
    }
  }
  return(mapping)
}

# all.pairwise.mappings = function(){

library(pathview)
library(SBGNview)

other.proteinid.to.pathwayCommons = function(){
    data("gene.idtype.list")
    gene.idtype.list = gene.idtype.list[!gene.idtype.list %in% c("TAIR","ORF")]  # those type of ids can't be handled by pathview
    pathwayCommonsIdPairs = t(expand.grid(gene.idtype.list,"pathwayCommons"))
    all.pairs = combn(gene.idtype.list,2)  # generate pairwise id types other than pathwayCommons
    all.pairs = cbind(pathwayCommonsIdPairs ,all.pairs)
    # all.pairs = t(expand.grid("UNIPROT",gene.idtype.list[!gene.idtype.list %in% c("UNIPROT","TAIR","ORF")]))
    id.mapping.all.list = list()
    tp = apply(
      all.pairs
      ,2
      ,function(pair){
        cat("\n\n\n\n\n")
        # gene.data <- sim.mol.data(mol.type = "gene", id.type = pair[1], nmol = 5000000,nexp=2)
        # in.ids = row.names(gene.data)
        print(pair)
        # print(head(in.ids))
        mapping = pairwiseMapping.for.sbgnIds(input.type=pair[1],output.type=pair[2],cpd.or.gene="gene")
        print(head(mapping))
        type.pair.name = paste(sort(pair),collapse="_")
        id.mapping.all.list[["gene"]][[type.pair.name]] <<- mapping
      }
    )
    lapply(id.mapping.all.list$gene,function(mapping){dim(mapping)})
    names(id.mapping.all.list$gene)
    save(id.mapping.all.list,file="./R/id.mapping/id.mapping.all.list.RData")
    dim(all.pairs)
    # }
    # all.pairwise.mappings()
}

other.proteinid.to.metacyc.gene = function(){
    data("gene.idtype.list")
    gene.idtype.list = gene.idtype.list[!gene.idtype.list %in% c("TAIR","ORF")]  # those type of ids can't be handled by pathview
    pathwayCommonsIdPairs = t(expand.grid(gene.idtype.list,"metacyc.SBGN"))
    # all.pairs = combn(gene.idtype.list,2)  # generate pairwise id types other than pathwayCommons
    all.pairs = pathwayCommonsIdPairs
    # all.pairs = t(expand.grid("UNIPROT",gene.idtype.list[!gene.idtype.list %in% c("UNIPROT","TAIR","ORF")]))
    # load("./R/id.mapping/id.mapping.all.list.RData")
    tp = apply(
      all.pairs
      ,2
      ,function(pair){
        cat("\n\n\n\n\n")
        # gene.data <- sim.mol.data(mol.type = "gene", id.type = pair[1], nmol = 5000000,nexp=2)
        # in.ids = row.names(gene.data)
        mapping = pairwiseMapping.for.sbgnIds(input.type=pair[1],output.type=pair[2],cpd.or.gene="gene")
        # if(pair[1]=="UNIPROT"){
            print(pair)
            print(head(in.ids))
            print(head(mapping))
            # stop()
        # }
        type.pair.name = paste(sort(pair),collapse="_")
        id.mapping.all.list[["gene"]][[type.pair.name]] <<- mapping
      }
    )
    # refs = id.mapping.all.list$gene$metacyc.SBGN_UNIPROT
    # dim(refs)
    # tail(refs)
    # head(refs[is.na(refs[,2]),])
      lapply(id.mapping.all.list$gene,function(mapping){dim(mapping)})
      names(id.mapping.all.list$gene)
      head(id.mapping.all.list$gene$metacyc.SBGN_PROSITE)
      head(id.mapping.all.list$gene$metacyc.SBGN_UNIGENE)
      head(id.mapping.all.list$gene$metacyc.SBGN_UNIPROT)
      all.mappings = names(id.mapping.all.list$gene)
      pathwayCommons.mappings = all.mappings[grepl("pathwayCommons",all.mappings)]
      pathwayCommons.mappings.all = list(gene=list())
      pathwayCommons.mappings.all$gene = id.mapping.all.list$gene[pathwayCommons.mappings]
      names(pathwayCommons.mappings.all)
      head(pathwayCommons.mappings.all$gene$ENSEMBL_pathwayCommons)
      save(pathwayCommons.mappings.all,file="./R/id.mapping/pathwayCommons.mappings.all.RData")
      lapply(pathwayCommons.mappings.all$gene,function(mapping){dim(mapping)})

      metacyc.SBGN.mappings = all.mappings[grepl("metacyc.SBGN",all.mappings)]
      metacyc.SBGN.mappings.all = list(gene=list())
      metacyc.SBGN.mappings.all$gene = id.mapping.all.list$gene[metacyc.SBGN.mappings]
      names(metacyc.SBGN.mappings.all)
      head(metacyc.SBGN.mappings.all$gene$metacyc.SBGN_SYMBOL)
      save(metacyc.SBGN.mappings.all,file="./R/id.mapping/metacyc.SBGN.mappings.all.RData")
      lapply(metacyc.SBGN.mappings.all$gene,function(mapping){dim(mapping)})
      # save(id.mapping.all.list,file="./R/id.mapping/id.mapping.all.list.RData")
    # dim(all.pairs)
    # }
    # all.pairwise.mappings()
}


all.compounds.to.chebi.pathwayCommons.metacyc = function(){
        load("/home/xiaoxi/0work/rSBGN/R/id.mapping/metacyc/mapping.tables/compound.to.chebi.mapping.all.RData")
    mapping.files.folder = "/mdedia/sf_work/input/id.mapping/unichem/ChEBI.db.names"
    
    metacyc.cpd.id = paste(compound.to.chebi.mapping.all[,2],compound.to.chebi.mapping.all[,3],sep=":@:")
    compound.to.chebi.mapping = cbind(compound.to.chebi.mapping.all[,1],metacyc.cpd.id)
    colnames(compound.to.chebi.mapping) = c("chebi","metacyc.SBGN")
    head(compound.to.chebi.mapping)
    files = list.files(mapping.files.folder,full.names = F,pattern="txt")
    # id.mapping.all.list
      id.mapping.all.list$compound = list()
      metacyc.SBGN.mappings.all$compound = list()
      pathwayCommons.mappings.all$compound = list()
    for (i in 1: length(files)){
    # for (i in 1:3){
        file.name = files[i]
        file.name = gsub(".txt","",file.name)
        id.pair = strsplit(file.name,"_")[[1]]
        cat("\n\n\n")
        print(files[i])
        file.content = read.table(paste(mapping.files.folder,files[i],sep="/"),header=T,sep="\t",stringsAsFactors = F)
        id.pair = sort(id.pair)
        id.pair = paste(id.pair,collapse = "_")
        print(id.pair)
        print(head(file.content))
        id.mapping.all.list$compound[[id.pair]] = file.content
        input.type = setdiff(colnames(file.content),c("chebi"))

        pathwayCommons.mapping = pairwiseMapping.for.sbgnIds(file.content[,input.type],input.type=input.type,output.type="pathwayCommons",cpd.or.gene="compound",id.to.chebi=file.content)
        pathwayCommons.id.pair = paste(sort(c(input.type,"pathwayCommons")),collapse = "_")
        pathwayCommons.mappings.all$compound[[pathwayCommons.id.pair]] = pathwayCommons.mapping
        id.mapping.all.list$compound[[pathwayCommons.id.pair]] = pathwayCommons.mapping

        metacyc.SBGN.mapping = pairwiseMapping.for.sbgnIds(file.content[,input.type],input.type=input.type,output.type="metacyc.SBGN",cpd.or.gene="compound",id.to.chebi=file.content)
        metacyc.SBGN.id.pair = paste(sort(c(input.type,"metacyc.SBGN")),collapse = "_")
        metacyc.SBGN.mappings.all$compound[[metacyc.SBGN.id.pair]] = metacyc.SBGN.mapping
        id.mapping.all.list$compound[[metacyc.SBGN.id.pair]] = metacyc.SBGN.mapping
    }
      # add chebi to sbgn id mappings
      library(SBGNview)
      head(chebi.to.pathway.commons.mapping)
      chebi.to.pathway.commons.mapping = SBGNview.uniprot.chebi.to.SBGN.id.mapping$compound
      chebi.to.pathway.commons.mapping[,"ChEBI"] = gsub("CHEBI:","",chebi.to.pathway.commons.mapping[,"ChEBI"] )
      colnames(chebi.to.pathway.commons.mapping) = c("pathwayCommons","chebi")
      pathwayCommons.mappings.all$compound$chebi_pathwayCommons = chebi.to.pathway.commons.mapping
      head(pathwayCommons.mappings.all$compound$chebi_pathwayCommons)
      id.mapping.all.list$compound$chebi_pathwayCommons = chebi.to.pathway.commons.mapping
      head(id.mapping.all.list$compound$chebi_pathwayCommons)

      head(compound.to.chebi.mapping)
      metacyc.SBGN.mappings.all$compound$chebi_metacyc.SBGN = compound.to.chebi.mapping
      head(metacyc.SBGN.mappings.all$compound$chebi_metacyc.SBGN)
      id.mapping.all.list$compound$chebi_metacyc.SBGN = compound.to.chebi.mapping
      head(id.mapping.all.list$compound$chebi_metacyc.SBGN)

      head(pathwayCommons.mappings.all$compound$actor_pathwayCommons)
      head(metacyc.SBGN.mappings.all$compound$actor_metacyc.SBGN)
      head(id.mapping.all.list$compound$actor_pathwayCommons)
      head(id.mapping.all.list$compound$chebi)

      names(pathwayCommons.mappings.all$compound)
      names(metacyc.SBGN.mappings.all$compound)
      names(id.mapping.all.list$compound)
      str(pathwayCommons.mappings.all)
      str(metacyc.SBGN.mappings.all)
      str(id.mapping.all.list)
      save(metacyc.SBGN.mappings.all,file="./R/id.mapping/metacyc.SBGN.mappings.all.RData")
      save(pathwayCommons.mappings.all,file="./R/id.mapping/pathwayCommons.mappings.all.RData")
    save(id.mapping.all.list,file="./R/id.mapping/id.mapping.all.list.RData")
}

add.uniref.metacyc.mapping = function(){
        load("/home/xiaoxi/0work/rSBGN/R/id.mapping/metacyc/merged.tables/uniref.pathway.sbgnBiopaxProtein.RData")
        head(uniref.pathway.sbgnBiopaxProtein)
        metacyc.SBGN_uniref = paste(uniref.pathway.sbgnBiopaxProtein[,"sbgnBiopaxProtein"],uniref.pathway.sbgnBiopaxProtein[,"pathway"],sep=":@:")
        uniref.to.metacyc = cbind(uniref.pathway.sbgnBiopaxProtein[,1],metacyc.SBGN_uniref)
        head(uniref.to.metacyc)
        colnames(uniref.to.metacyc) = c("uniref","metacyc.SBGN")
        head(uniref.to.metacyc)
      metacyc.SBGN.mappings.all$gene$metacyc.SBGN_uniref = uniref.to.metacyc
      head(metacyc.SBGN.mappings.all$gene$metacyc.SBGN_uniref)
      id.mapping.all.list$gene$metacyc.SBGN_uniref = uniref.to.metacyc
      head(id.mapping.all.list$gene$metacyc.SBGN_uniref)
      # save(metacyc.SBGN.mappings.all,file="./R/id.mapping/metacyc.SBGN.mappings.all.RData")

      metacyc.gene.ids = names(metacyc.SBGN.mappings.all$gene)
      metacyc.gene.ids = metacyc.gene.ids[grepl("metacyc",metacyc.gene.ids)]
      metacyc.gene.ids =gsub("_metacyc.SBGN","",metacyc.gene.ids)
      metacyc.gene.ids =gsub("metacyc.SBGN_","",metacyc.gene.ids)
      metacyc.gene.ids
      metacyc.SBGN.mappings.all$metacyc.mappable.gene.ids = metacyc.gene.ids
      names(metacyc.SBGN.mappings.all)

      metacyc.compound.ids = names(metacyc.SBGN.mappings.all$compound)
      metacyc.compound.ids = metacyc.compound.ids[grepl("metacyc",metacyc.compound.ids)]
      metacyc.compound.ids =gsub("_metacyc.SBGN","",metacyc.compound.ids)
      metacyc.compound.ids =gsub("metacyc.SBGN_","",metacyc.compound.ids)
      metacyc.compound.ids
      metacyc.SBGN.mappings.all$metacyc.mappable.compound.ids = metacyc.compound.ids
      names(metacyc.SBGN.mappings.all)

      pathwayCommons.gene.ids = names(pathwayCommons.mappings.all$gene)
      pathwayCommons.gene.ids = pathwayCommons.gene.ids[grepl("pathwayCommons",pathwayCommons.gene.ids)]
      pathwayCommons.gene.ids =gsub("_pathwayCommons","",pathwayCommons.gene.ids)
      pathwayCommons.gene.ids =gsub("pathwayCommons_","",pathwayCommons.gene.ids)
      pathwayCommons.gene.ids
      pathwayCommons.mappings.all$pathwayCommons.mappable.gene.ids = pathwayCommons.gene.ids
      names(pathwayCommons.mappings.all)

      pathwayCommons.compound.ids = names(pathwayCommons.mappings.all$compound)
      pathwayCommons.compound.ids = pathwayCommons.compound.ids[grepl("pathwayCommons",pathwayCommons.compound.ids)]
      pathwayCommons.compound.ids =gsub("_pathwayCommons","",pathwayCommons.compound.ids)
      pathwayCommons.compound.ids =gsub("pathwayCommons_","",pathwayCommons.compound.ids)
      pathwayCommons.compound.ids
      pathwayCommons.mappings.all$pathwayCommons.mappable.compound.ids = pathwayCommons.compound.ids
      names(pathwayCommons.mappings.all)
      str(pathwayCommons.mappings.all)

      save(pathwayCommons.mappings.all,file="./R/id.mapping/pathwayCommons.mappings.all.RData")
      save(metacyc.SBGN.mappings.all,file="./R/id.mapping/metacyc.SBGN.mappings.all.RData")
    save(id.mapping.all.list,file="./R/id.mapping/id.mapping.all.list.RData")

}



# entrenz.to.pathwayCommons = function(uniprot.to.pathwayCommons.species.matrix,gene.to.pathwayCommons.each.pair.RData.folder){
all.gene.id.types.to.pathwayCommons = function(uniprot.to.pathwayCommons.species.matrix,gene.to.pathwayCommons.each.pair.RData.folder){
    data("gene.idtype.list")
    gene.idtype.list = gene.idtype.list[!gene.idtype.list %in% c("TAIR","ORF")]  # those type of ids can't be handled by pathview
    gene.idtype.list = c(gene.idtype.list,"ENTREZID")
    # all.pairs = as.matrix(c("ENTREZID","pathwayCommons"))
    all.pairs = t(expand.grid(gene.idtype.list,"pathwayCommons"))
    # all.pairs = combn(gene.idtype.list,2)  # generate pairwise id types other than pathwayCommons
    all.pairs
    id.mapping.all.list = list()
  # load("C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/R/id.mapping/pathway.commons/uniprot.to.pathwayCommons.all.species.list.RData")
  # SBGNview.uniprot.chebi.to.SBGN.id.mapping = uniprot.to.pathwayCommons.all.species.list
  # load("C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/id.mapping/mapping.each.pair/for.use.all/both/pathwayCommons_UNIPROT.RData")
    SBGNview.uniprot.chebi.to.SBGN.id.mapping = list()
  SBGNview.uniprot.chebi.to.SBGN.id.mapping$gene = uniprot.to.pathwayCommons.species.matrix
  head(SBGNview.uniprot.chebi.to.SBGN.id.mapping[[1]])
  
    tp = apply(
      all.pairs
      ,2
      ,function(pair){
        cat("\n\n\n\n\n")
        # gene.data <- sim.mol.data(mol.type = "gene", id.type = pair[1], nmol = 5000000,nexp=2)
        # in.ids = row.names(gene.data)
        print(pair)
        type.pair.name = paste(sort(pair),collapse="_")
        out.rd.file = paste(gene.to.pathwayCommons.each.pair.RData.folder,"/",type.pair.name,".RData",sep="")
        if(file.exists(out.rd.file)){
                return()
        }
        # print(head(in.ids))
        mapping = pairwiseMapping.for.sbgnIds(input.type=pair[1],output.type=pair[2],cpd.or.gene="gene",SBGNview.uniprot.chebi.to.SBGN.id.mapping = SBGNview.uniprot.chebi.to.SBGN.id.mapping)
        print(head(mapping))
        # stop()
        # print(mapping)
        mapping = mapping[!is.na(mapping[,1]),]
        mapping = mapping[!is.na(mapping[,2]),]
        mapping = mapping[!is.na(mapping[,3]),]
        mapping.list = list()
        mapping.list$gene = list()
        mapping.list$gene[[type.pair.name]] = mapping
        print(dim(mapping.list$gene[[1]]))
        print(type.pair.name)
        print(str(mapping.list))
        print(head(mapping.list[[1]][[1]]))
        save(mapping.list,file = out.rd.file)
        id.mapping.all.list[[type.pair.name]] <<- mapping
      }
    )
    
}

entrenz.to.metacyc.gene = function(){
    # data("gene.idtype.list")
    # gene.idtype.list = gene.idtype.list[!gene.idtype.list %in% c("TAIR","ORF")]  # those type of ids can't be handled by pathview
    # pathwayCommonsIdPairs = t(expand.grid(gene.idtype.list,"metacyc.SBGN"))
    # all.pairs = combn(gene.idtype.list,2)  # generate pairwise id types other than pathwayCommons
    # all.pairs = as.matrix(c("ENTREZID","metacyc.SBGN"))
    data("gene.idtype.list")
    gene.idtype.list = gene.idtype.list[!gene.idtype.list %in% c("TAIR","ORF")]  # those type of ids can't be handled by pathview
    # all.pairs = as.matrix(c("ENTREZID","pathwayCommons"))
    all.pairs = t(expand.grid(gene.idtype.list,"metacyc.SBGN"))
    # all.pairs = combn(gene.idtype.list,2)  # generate pairwise id types other than pathwayCommons
    all.pairs
    id.mapping.all.list = list()
    # all.pairs = cbind(pathwayCommonsIdPairs ,all.pairs)
    # all.pairs = t(expand.grid("UNIPROT",gene.idtype.list[!gene.idtype.list %in% c("UNIPROT","TAIR","ORF")]))
    # id.mapping.all.list = list()
    
    tp = apply(
      all.pairs
      ,2
      ,function(pair){
        cat("\n\n\n\n\n")
        # gene.data <- sim.mol.data(mol.type = "gene", id.type = pair[1], nmol = 5000000,nexp=2)
        # in.ids = row.names(gene.data)
        print(pair)
        type.pair.name = paste(sort(pair),collapse="_")
        output.file = paste("C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/id.mapping/mapping.each.pair/for.use.all/",type.pair.name,".RData",sep="")
        if(file.exists(output.file)){
                return()
        }
        # print(head(in.ids))
        mapping = pairwiseMapping.for.sbgnIds(input.type=pair[1],output.type=pair[2],cpd.or.gene="gene")
        print(head(mapping))
        mapping.list = list()
        mapping.list$gene = list()
        mapping.list$gene[[type.pair.name]] = mapping
        print(dim(mapping.list$gene[[1]]))
        print(type.pair.name)
        print(table(mapping[,"species"]))
        # stop()
        save(mapping.list,file = output.file)
        # print(mapping)
        # die()
        id.mapping.all.list[[type.pair.name]] <<- mapping
      }
    )
    # lapply(id.mapping.all.list$gene,function(mapping){dim(mapping)})
#     names(id.mapping.all.list)
#     table(id.mapping.all.list[[1]][,3])
#     table(id.mapping.all.list[[2]][,3])
#     head(id.mapping.all.list[[2]])
#     entrezid.to.metacyc.pathwayCommons.with.species.info = id.mapping.all.list
#     save(id.mapping.all.list,file="/media/sf_work/papers/SBGNview/package/data.package/entrezid.to.metacyc.pathwayCommons.with.species.info.RData")
#     # dim(all.pairs)
}




























