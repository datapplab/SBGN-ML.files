  #parse ports
# function:
# from pathwayCommons biopax files, generate mapping tables: 1. macromolecule id to uniprot. 2 simple chemical id to chebi
  library(xml2)
  biopax.to.id.mapping = function(biopax.xml,node.type){
        # function: from a xml file, for a type of node(node.type: small molecule/protein/protein reference), find their id mapping tables: c("id.q","display.name","id.found.entityReference","id.found.type").
        pattern = paste(".//",node.type,sep="")
        nodes.info = xml_find_all(biopax.xml,pattern)
        mapping = c("id.q","display.name","id.found.entityReference","id.found.type")
        # nodes.info = nodes.info[1:100000]
        # print(nodes.info[1:10])
        ids = lapply(
              nodes.info
              ,function(glyph){
                    glyph.id = xml_attr(glyph,"ID")
                    if(is.na(glyph.id)){
                        return()
                        # print(glyph)
                        # stop("glyph without id")
                    }
                    if(node.type %in% c("SmallMolecule","Protein")){
                        id.node.name = "entityReference"
                    }else if (node.type %in% c("ProteinReference","SmallMoleculeReference")){
                        id.node.name = "memberEntityReference"
                    }
                    displayName = ""
                    entityReference = ""
                   sapply(
                      xml_children(glyph)
                      ,function(child){
                          child.name = xml_name(child)
                          if(child.name == "displayName"){
                              displayName <<- xml_text(child)
                          }
                      }
                    )
                   sapply(
                      xml_children(glyph)
                      ,function(child){
                          child.name = xml_name(child)
                          if(child.name == id.node.name){
                              entityReference <<- xml_attr(child,"resource")
                              entityReference <<- gsub("#","",entityReference)
                              entityReference <<- gsub("http://identifiers.org/chebi/","",entityReference)
                              entityReference <<- gsub("http://identifiers.org/uniprot/","",entityReference)
                              if.ref = "mapped.id"
                              if(grepl("Reference_",entityReference)){
                                  if.ref = "ref"
                              }
                              mapping <<- rbind(mapping,c(glyph.id,displayName,entityReference,if.ref))
                          }
                      }
                    )
                   if(entityReference == ""){
                        # cat("\n\n\n  glyph with no database info\n\n")
                        # print(glyph)
                        mapping <<- rbind(mapping,c(glyph.id,displayName,"no.entity.ref","no.ref.for.ref"))
                   }
                    if(glyph.id %in% c("ProteinReference_986d21d017dc2de1b0ee9a5b806f3f7f","ProteinReference_5384d757cc969c5e4ea4ab4f51881448","SmallMoleculeReference_6a7acf777207de3b96410ac55de4666e")){
                        # print(glyph)
                    }
              }
        )
        colnames(mapping) = mapping[1,]
        mapping = mapping[-1,]
        return(mapping)
  }

find.id.mapping.each.biopax.file.pathwayCommons10 = function(biopax.file){
      # for each biopax file, find id mappings for different molecular types: small molecule or proteins
      # output table: c("id.q","database.ids","display.name_molecule.to.ref")
      biopax.xml = read_xml(biopax.file)
      xml_attrs(biopax.xml) = NULL # Remove root node attribute. This is necessary Otherwise xml2 won't find the nodes when using xml_find_all.
      print("molecule to ref: small molecules")
      small.molecule.to.ref = biopax.to.id.mapping(biopax.xml,"SmallMolecule")
      print("refs to database ids: small molecules")
      small.molecule.ref.to.chebi = biopax.to.id.mapping(biopax.xml,"SmallMoleculeReference")
      print("molecule to ref: proteins")
      protein.to.ref = biopax.to.id.mapping(biopax.xml,"Protein")
      print("refs to database ids: proteins")
      protein.ref.to.uniprot = biopax.to.id.mapping(biopax.xml,"ProteinReference")
      # for small chemicals, generate full table containing chebi ids
      small.molecules.all  = merge(small.molecule.to.ref,small.molecule.ref.to.chebi,by.x="id.found.entityReference",by.y="id.q",all=T,suffixes = c("_molecule.to.ref","_ref.to.id"))
      colnames(small.molecules.all)[1] = "id.found.entityReference_molecule.to.ref"
      database.ids = as.character(small.molecules.all[,"id.found.type_ref.to.id"])
      if.mapped.id = as.character(small.molecules.all[,"id.found.type_molecule.to.ref"]) == "mapped.id" & !is.na(as.character(small.molecules.all[,"id.found.type_molecule.to.ref"]))
      database.ids[if.mapped.id] = as.character(small.molecules.all[if.mapped.id,"id.found.entityReference_molecule.to.ref"])
      if.mapped.id = as.character(small.molecules.all[,"id.found.type_ref.to.id"]) == "mapped.id" & ! is.na(as.character(small.molecules.all[,"id.found.type_ref.to.id"]) )
      database.ids[if.mapped.id] = as.character(small.molecules.all[if.mapped.id,"id.found.entityReference"])
      small.molecules.all$database.ids = database.ids
      small.molecules.all = small.molecules.all[,c("id.q","database.ids","display.name_molecule.to.ref")]

      # for proteins ,generate full table containing uniprot ids
      proteins.all = merge(protein.to.ref,protein.ref.to.uniprot,by.x="id.found.entityReference",by.y="id.q",all=T,suffixes = c("_molecule.to.ref","_ref.to.id"))
      colnames(proteins.all)[1] = "id.found.entityReference_molecule.to.ref"
      database.ids = as.character(proteins.all[,"id.found.type_ref.to.id"])
      if.mapped.id = as.character(proteins.all[,"id.found.type_molecule.to.ref"]) == "mapped.id" & !is.na(as.character(proteins.all[,"id.found.type_molecule.to.ref"]))
      database.ids[if.mapped.id] = as.character(proteins.all[if.mapped.id,"id.found.entityReference_molecule.to.ref"])
      if.mapped.id = as.character(proteins.all[,"id.found.type_ref.to.id"]) == "mapped.id" & ! is.na(as.character(proteins.all[,"id.found.type_ref.to.id"]) )
      database.ids[if.mapped.id] = as.character(proteins.all[if.mapped.id,"id.found.entityReference"])
      proteins.all$database.ids = database.ids
      proteins.all = proteins.all[,c("id.q","database.ids","display.name_molecule.to.ref")]
      return(list(
            proteins.all = proteins.all
            ,small.molecules.all = small.molecules.all
      ))
}
mappings.all.files.to.single.column = function(mapping.all.files){
      # in the input table mapping.all.files, each database(biopax file) has its own two columns(e.g. panther has two columns, reactome has another two columns, resulting from the merge function in "biopax.files.to.id.mapping"), this function:
      #   1. check if any entity exists in multiple biopax files
      #   2. combine different databases into three columns: c(entity id,  database id, display name)
        duplicated.ids = list()
        mappings = apply(mapping.all.files
                    ,1
                    ,function(entity){
                        entity = entity[!is.na(entity)]
                        n.nas = sum(!is.na(entity))
                        if(n.nas>3){
                            duplicated.ids[[entity[1]]] <<- entity
                        }else{
                            return(c(entity))
                        }
                    }
                )
        mappings = t(mappings)
        colnames(mappings) = c("entity.id","database.id","display.name")
        mappings.no.database.id = mappings[mappings[,2]=="no.ref.for.ref",]
        mappings.database.id.name = mappings[mappings[,2]!="no.ref.for.ref",]
        return(list(
            duplicated.ids = duplicated.ids
            ,mappings.no.database.id = mappings.no.database.id
            ,mappings.database.id.name = mappings.database.id.name
        ))
}

# id.mapping.from.multiple.biopax.files = function()
biopax.files.to.id.mapping = function(biopax.folder){
    file.extension = ".owl"
    biopax.files = list.files(biopax.folder,pattern=file.extension,full.names = F)  # All pathways in each type of database has only one file(e.g. panther, reactome etc.)
    small.molecules.all = data.frame()
    proteins.all = data.frame()
    tp = sapply(
        biopax.files
        # sbgn.files[17:19]
        ,function(biopax.file){
            cat("\n\n\n\n biopax file")
            print(biopax.file)
            input.file = paste(biopax.folder,"/",biopax.file,sep="")
            result.list = find.id.mapping.each.biopax.file.pathwayCommons10(input.file)
            small.molecule.mapping = result.list$small.molecules.all
            proteins.mapping = result.list$proteins.all
            if(nrow(proteins.all)==0){
                small.molecules.all <<- small.molecule.mapping
                proteins.all <<- proteins.mapping
            }else{
                small.molecules.all <<- merge(small.molecule.mapping,small.molecules.all,by.x="id.q",by.y="id.q",all=T,suffixes = c("",paste("_",biopax.file,sep="")))
                proteins.all <<-  merge(proteins.all,proteins.mapping,by.x="id.q",by.y="id.q",all=T,suffixes = c("",paste("_",biopax.file,sep="")))
            }
            print(head(small.molecule.mapping))
            print(head(proteins.mapping))
        }
    )
    result.list = mappings.all.files.to.single.column(small.molecules.all)
    small.molecules.duplicated.ids = result.list$duplicated.ids
    small.molecules.to.chebi = result.list$mappings.database.id.name
    small.molecules.no.database.ids = result.list$mappings.no.database.id
    saveRDS(result.list,paste(biopax.folder,"/small.molecules.mapping.list.rds",sep=""))

    result.list = mappings.all.files.to.single.column(proteins.all)
    proteins.duplicated.ids = result.list$duplicated.ids
    proteins.to.uniprot = result.list$mappings.database.id.name
    proteins.no.database.ids = result.list$mappings.no.database.id
    saveRDS(result.list,paste(biopax.folder,"/proteins.mapping.list.rds",sep=""))
    dim(proteins.no.database.ids)
    dim(proteins.to.uniprot)
    return(biopax.folder)
}
