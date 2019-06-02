data.folder = "C:/xiaoxi/work/papers/SBGNview/package/backup/sbgnview.data/2019.5.22/original.variable.name.mapping.list.mapping.table"
output.folder = "C:/xiaoxi/work/papers/SBGNview/package/backup/sbgnview.data/2019.5.22/original.variable.name.mapping.list.mapping.table.unique.name"
data.files = list.files(data.folder,full.names = FALSE)
for(i in seq_along(data.files)){
        full.path = paste0(data.folder,"/",data.files[i])
        output.file.name = gsub(" ","",data.files[i])
        pair.name = gsub(".RData","",output.file.name)
        cat("\n\n\n\n")
        print(pair.name)
        print("loading data")
        loaded.name = tryCatch(load(full.path),warning=function(w) "no such data") 
        if(loaded.name == "mapping.list"){
                mapping.table = mapping.list[[1]][[1]]
        }else if(loaded.name == "mapping.table"){
        }else{
                print("data not found")
        }
        # mapping.table = mapping.table[, which(!colnames(mapping.table) %in% "species")]
        mapping.table[,1] = as.character(mapping.table[,1])
        mapping.table[,2] = as.character(mapping.table[,2])
        print("assigning data")
        assign(pair.name,mapping.table)
        print(head(get(pair.name)))
        output.full.name = paste0(output.folder,"/",output.file.name)
        print("saving data")
        save(list = pair.name, file = output.full.name, compress = "xz")
        # stop()
}
