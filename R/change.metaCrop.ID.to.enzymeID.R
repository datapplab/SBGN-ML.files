library(xml2)
metaCrop.sbgn.folder = "C:/xiaoxi/work/input/metacrop"
metaCrop.sbgn.files = list.files(metaCrop.sbgn.folder,full.names = T)
for(i in 1:length(metaCrop.sbgn.folder)){
        metaCrop.sbgn.file = metaCrop.sbgn.files[i]
        sbgn.xml = read_xml(metaCrop.sbgn.file)
        xml_attrs(sbgn.xml) = NULL
        glyphs = xml_find_all(sbgn.xml,".//glyph")
        for (i in 1:length(glyphs)){
                glyph = glyphs[[i]]
                glyph.id = xml_attr(glyph,"id")
                label.node = xml_find_all(glyph,".//label")
                print(label.node)
                label = xml_attr(label.node,"text")
                print(c(glyph.id,label))
                stop()
                
        }
}