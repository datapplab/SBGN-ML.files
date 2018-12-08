
# given two vectors of characters/words, this function find the best match between words in the two vectors, by spliting a word into sub-strings and use "jaccard" to calculate two words' similarity (similarity between thier sub-string vectors). It for one name, the function will output the pair with largest similarity score.
match.names = function(names.1,names.2,output.file){
    
    names.1 = tolower(names.1)
    names.2 = tolower(names.2)
    names.1 = unique(names.1)
    names.2 = unique(names.2)
    head(names.1)
    head(names.2)
    length(names.1)
    length(names.2)
    
    all.pairs = expand.grid(names.1,names.2)
    head(all.pairs[,1:2])
    
    ot = apply(all.pairs
          ,1
          ,function(pair){
              # print(pair)
              p1n = pair[1]
              p2n = pair[2]
                p1 = strsplit(p1n,"[^a-zA-Z0-9]+",perl=T)[[1]]
                p1 = unique(p1)
                p1
                p2 = strsplit(p2n,"[^a-zA-Z0-9]+",perl=T)[[1]]
                p2 = unique(p2)
                p2
                inter = length(intersect(p1,p2))
                uni = length(unique(union(p1,p2)))
                jc = inter/uni
                return(c(p1n,p2n,jc))
          }
          )
    ot = t(ot)
    head(ot)
    
    
    ot=as.data.frame(ot,stringsAsFactors = F)
    names(ot) = c("name.in.chebi","name.abundance","overlap")
    ot$overlap = as.numeric(ot$overlap)
    
    
    
    names.2[grepl("choe",names.2,ignore.case = T)]
    ot[ot$name.in.chebi == "choe(16:1)" & ot$name.abundance == "choe.16.1.",]
    
    ot1 = ot[1:2000,]
    ot1 = ot
    mapped.table = by(
        ot1
        ,as.factor(ot1$name.abundance)
        ,function(df){
            if.max = df$overlap == max(df$overlap)
            mx.df = df[if.max,]
            mx.df = mx.df[1,]
        }
        
    )
    mapped.table = do.call(rbind,mapped.table)
    row.names(mapped.table) = c()
    head(mapped.table)
    dim(mapped.table)
    mapped.table[mapped.table[,"name.abundance"] == "choe.16.1.",]
    write.table(mapped.table,output.file,sep="\t",row.names=F,quote=F)
    return(mapped.table)
}
