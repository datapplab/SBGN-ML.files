Each RData here is a list containing one pairwise ID mapping table. 
The name of the list is "mapping.list". Try this:
>variable.name = load("pathwayCommons_UNIPROT.RData")
>print(variable.name)
>str(mapping.list)

Data generation:

1. Compound:
    1.1 pathwayCommons:
        1.1.1 Map ChEBI to pathwayCommons Biopax protein ID.
                The mapping was extracted from BIOPAX files such as: http://www.pathwaycommons.org/archives/PC2/v10/PathwayCommons10.reactome.BIOPAX.owl.gz
        1.1.2 Use UniChem to map other ID types to ChEBI, then map to pathwayCommons Biopax protein IDs.
        
    1.2 MetaCyc
        1.2.1 Map ChEBI to MetaCyc Biopax protein ID.
                Extracted from downloaded BIOPAX files. For example : https://metacyc.org/META/pathway-biopax?type=3&object=PWY-7754
        1.2.2 Use UniChem to map other ID types to ChEBI, then map to MetaCyc Biopax protein IDs.
              

2. Protein/gene:
    2.1 pathwayCommons:
        2.1.1 Map UNIPROT to pathwayCommons Biopax protein ID.
              2.1.1.1 The mapping was extracted from BIOPAX files such as: http://www.pathwaycommons.org/archives/PC2/v10/PathwayCommons10.reactome.BIOPAX.owl.gz
                        pathwayCommons only curated huamn pathways, therefore all the UNIPROT IDs are from human.
              2.1.1.2 We used the ortholog mapping from Reactome to map UNIPROT of other species to UNIPROT of human:
                        The mapping was extracted from https://reactome.org/download/current/UniProt2Reactome_PE_Pathway.txt
              2.1.1.3  Then map UNIPROT from other species to pathwayCommons BIOPAX protein ID through human UNIPROT
                Result is in this RData "pathwayCommons_UNIPROT.RData"
        2.1.2 Map other IDs to pathwayCommons biopax protein IDs through UNIPROT, using pathview
                        Note: currently pathwview supports only several species (not the same as the species in Reactome). Therefore some species in UNIPROT-TO-pathwayCommons mapping may be lost.
    2.2 MetaCyc
        2.2.1 Map UNIPROT to MetaCyc Biopax protein ID.
                Extracted from downloaded BIOPAX files. For example : https://metacyc.org/META/pathway-biopax?type=3&object=PWY-7754
                It has UNIPROT from all supported species.
                Result is in this RData "metacyc.SBGN_UNIPROT.RData"
        2.2.2 Use pathview to map other ID types to UNIPROT, then map to MetaCyc Biopax protein IDs.
         Note: currently pathwview supports only several species (not the same as the species in Reactome). Therefore some species in UNIPROT-TO-metacyc mapping may be lost.
