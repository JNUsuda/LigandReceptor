## ligand-receptor circos neurotransmitter ------

#devtools::install_github("jinworks/CellChat")
#devtools::install_github("jokergoo/circlize")
library("CellChat")
library("circlize")
library("ComplexHeatmap")
library("dplyr")
library("tidyr")

# get CellChat list of ligand-receptors
db = CellChat::CellChatDB.human
db2 = db$interaction
#write.table(db2, file = "./CellChatdb.tsv", sep = "\t", row.names = F)

# if unable to load CellChat, use the table:
db2 = read.table(file = "./CellChatdb.tsv", sep = "\t", header = T)

# check for consistency in is.neurotransmitter annotation for pathways
df3 = db2 %>% 
  dplyr::group_by(pathway_name, is_neurotransmitter) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = is_neurotransmitter, values_from = n) %>% 
  dplyr::rename(true = 'TRUE', false = 'FALSE') %>% 
  filter(!is.na(true)) %>% # filter pathways which have a TRUE for 'is_neurotransmitter'
  filter(is.na(false)) # filter pathways which don't have a FALSE for 'is_neurotransmitter'

# filter CellChat db for neurotransmitter pathways 
db3 = db2 %>% dplyr::filter(pathway_name %in% df3$pathway_name)

# tables for circos plot/chord diagram
# load table with a list of genes of interest (one column named "Genes")
tab2 = read.table(file = "./LR_circos_example.tsv", header = T, sep = "\t")
colnames(tab2) = "Genes"
my_genes = tab2$Genes

# prepare empty tables
# table for CellChatdb filtered for neurotransmitter ligands or receptors
tab_allu = db3[0,] 
tab_allu = tab_allu %>% # make columns to rbind
  #mutate(disease = NA) %>% 
  mutate(LR = NA)
# table for genes from your list which were neurotransmitter ligands or receptors
pres = data.frame(pres = character())


# filter rows that have the genes of interest
  ## ligands
  dbL = db3 %>%
    separate_rows(ligand.symbol, sep = ", ") # separate ligand names 
  dbL = dbL %>% dplyr::filter(ligand.symbol %in% my_genes) # filter which are in the list of my_genes
  presL = dbL %>% dplyr::distinct(ligand.symbol) %>% # get list of which were in the list of my_genes
    dplyr::rename(pres = ligand.symbol) #rename column to 'pres' (presence, which were in the list of my_genes)
  dbL = dbL %>% distinct(interaction_name, .keep_all = T)
  tabL = db3 %>% dplyr::filter(interaction_name %in% dbL$interaction_name) # filter db 
  
  ## receptors
  dbR = db3 %>%
    separate_rows(receptor.symbol, sep = ", ") # separate receptor names 
  dbR = dbR %>% dplyr::filter(receptor.symbol %in% my_genes) # filter which are in the list of my_genes
  presR = dbR %>% dplyr::distinct(receptor.symbol) %>% 
    dplyr::rename(pres = receptor.symbol) # rename column to 'pres' (presence, which were in the list of my_genes)
  dbR = dbR %>% distinct(interaction_name, .keep_all = T)
  tabR = db3 %>% dplyr::filter(interaction_name %in% dbR$interaction_name) # filter db table for interactions where receptors were in the list
  
  # if there are results (has rows), add identification cols and append to tab_allu 
  if (nrow(tabL) != 0) {
    #tabL$disease = grupao
    tabL$LR = "LIGAND"
    tab_allu = rbind(tab_allu, tabL)
    pres = rbind(pres, presL)}
  
  if (nrow(tabR) != 0) {
    #tabR$disease = grupao
    tabR$LR = "RECEPTOR"
    tab_allu = rbind(tab_allu, tabR)
    pres = rbind(pres, presR)}

  
# reorder columns
tab_allu2 = tab_allu %>%
  relocate(pathway_name, ligand.symbol, receptor.symbol, 
           #disease, 
           LR, is_neurotransmitter,
           .before = 1)
# add count for each line
tab_allu2$num_rows = 1


tab_allu3 = tab_allu2 %>% 
  separate_rows(ligand.symbol, sep = ", ") %>% 
  separate_rows(receptor.symbol, sep = ", ") %>% 
  distinct(pathway_name, ligand.symbol, receptor.symbol, #disease, 
           LR, num_rows)


tab_allu4 = tab_allu3 %>%  
  #dplyr::filter(disease == 'PD') %>% 
  dplyr::rename(ligand = ligand.symbol, receptor = receptor.symbol) %>% 
  group_by(ligand, receptor, pathway_name) %>% 
  summarize(n = n()) %>% 
  mutate(cores = if_else(pathway_name == "GABA-A", "#74a892", false = 
                 if_else(pathway_name == "GABA-B", "#008585", false = ##008585
                 if_else(pathway_name == "Glutamate", "#e5c185", false =
                 if_else(pathway_name == "Glycine", "#c7522a", false = 
                 if_else(pathway_name == "TRH", "#db7d42", false =
                 if_else(pathway_name == "CALC", "#cee583", false =
                 if_else(pathway_name == "5-HT", "#195382", false =
                 if_else(pathway_name == "Dopamine", "#057dcd", false =
                 if_else(pathway_name == "SerotoninDopamin", "#74c1ed", false = NA
          )))))))))) %>% 
  mutate(pathway_name = factor(pathway_name, 
                               levels = c("GABA-A","GABA-B","Glutamate",
                                          "Glycine","TRH", "CALC", "5-HT",
                                          "Dopamine", "SerotoninDopamin")))


# cor da barrinha externa
grid.col = tab_allu4 %>% pivot_longer(cols = c(ligand, receptor))
grid.col2 = setNames(grid.col$cores, grid.col$value)

# cor dos tracks (conexoes)
col = tab_allu4$cores

# ordem dos setores
orderL = tab_allu4 %>% ungroup() %>% arrange(pathway_name) %>% distinct(ligand)
orderR = tab_allu4 %>% ungroup() %>% arrange(pathway_name) %>% distinct(receptor)
order1 = c(rev(orderL$ligand), (orderR$receptor))


png(filename = "../../../figures/circos/example1.png", res = 300, units = "px", width = 2000, height = 2000)
pdf(file = "../../../figures/circos/AD1.pdf", width = 6.5, height = 6.5)

#svglite(filename = "../../../figures/circos/PD1.svg",  width = 2000, height = 2000)

# plot two side by side
#par(mar = c(0, 0, 0, 0), mfrow = c(1, 2))

# Restart circular layout parameters
circos.clear()
# circos.par(start.degree = 90) # to rotate plot
chordDiagram(tab_allu4, 
             grid.col = grid.col2, # set colors
             preAllocateTracks = 1, 
             order = order1, # set sector order
             annotationTrackHeight = mm_h(2), # altura da barrinha externa
             annotationTrack =  "grid")
# add labels:
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  ## highlight specific ligands/receptors with bold font:
  if (sector.name %in% pres$pres) {
    circos.text(mean(xlim), ylim[1] + .1, # label position 
                sector.name, 
                cex = 0.9, # font size
                facing = "clockwise", 
                font = 2, # bold
                niceFacing = TRUE, adj = c(-0, 0.5))} 
  
  ## add the other labels in normal font:
  else{
    circos.text(mean(xlim), ylim[1] + .1, # label position 
                sector.name, 
                cex = 0.9, # font size
                facing = "clockwise", 
                niceFacing = TRUE, adj = c(-0, 0.5))}
}, bg.border = NA)

# add legend
tab_allu5 = tab_allu4 %>% ungroup() %>% distinct(pathway_name, .keep_all = T)
lgd_points = Legend(at = as.character(tab_allu5$pathway_name), 
                    type = "grid", #squares
                    legend_gp = gpar(fill = tab_allu5$cores),
                    title_position = "topleft", 
                    gap = unit(0.3, "cm"),  
                    title = "Pathway")

draw(lgd_points, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))

dev.off()

