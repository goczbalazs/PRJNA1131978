library(ggrepel)
library(tidyverse)
library("edgeR")
library('RColorBrewer')
library(writexl)
library(biomaRt)
library(ggThemeAssist)

files <- list.files(pattern = 'no_multi.txt')    
for(i in files) {
  x <-read_delim(i,                     delim = "\t", escape_double = FALSE, 
                 col_types = cols(Chr = col_skip(), Start = col_skip(), 
                                  End = col_skip(), Strand = col_skip(), 
                                  Length = col_skip()), trim_ws = TRUE, 
                 skip = 1)
  assign(i,x)  
}

########Coverage of ChINsIHC and ChINsZsG#####

CPU_1000 <- featureCounts_GRCm39.107_CPU_1000.txt %>% 
  rename_at(vars(contains('Final_GRCm39_107_')), list( ~ gsub('Final_GRCm39_107_', '', .))) %>%  
  rename_at(vars(contains('_trimmed_finalAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmed_finalAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ gsub('-', '_', .))) %>% 
  data.frame(row.names = 1)


CPU_300_data_original <- featureCounts_GRCm39.107_s2_no_multi.txt %>% 
  rename_at(vars(contains('GRCm39_107_')), list( ~ gsub('GRCm39_107_', '', .))) %>%  
  rename_at(vars(contains('_trimmed_finalAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmed_finalAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ gsub('-', '_', .))) %>% 
  data.frame(row.names = 1) 

CPU300_data_original <- CPU_300_data_original %>% 
  dplyr::select(smpl01, smpl03, smpl05, smpl07) %>% 
  rename(CHAT300_1 = smpl01) %>% 
  rename(CHAT300_3 = smpl03) %>% 
  rename(CHAT300_5 = smpl05) %>% 
  rename(CHAT300_7 = smpl07)

CPU_300_with_name <- CPU300_data_original %>%
  rownames_to_column(var = "Geneid")

CPU_1000_with_name <- CPU_1000 %>% 
  rownames_to_column(var = "Geneid")

data_1000vs300_107 <- CPU_1000_with_name %>% 
  inner_join(CPU_300_with_name) %>% 
  data.frame(row.names = 1)

ANN_file= getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                filters = c('ensembl_gene_id') ,
                values = CPU_1000_with_name$Geneid,
                mart = ensembl107)

CPU_1000_cpm = as.data.frame(cpm(CPU_1000))
colnames(CPU_1000_cpm) = paste0('cpm_', colnames(CPU_1000_cpm))
CPU_1000_cpm <- rownames_to_column(CPU_1000_cpm,var = "Geneid")

CPU_300_cpm = as.data.frame(cpm(CPU300_data_original))
colnames(CPU_300_cpm) = paste0('cpm_', colnames(CPU_300_cpm))
CPU_300_cpm <- rownames_to_column(CPU_300_cpm,var = "Geneid")

CPU_1000_all <- CPU_1000_cpm %>% 
  inner_join(ANN_file, by = c("Geneid" = "ensembl_gene_id")) %>% 
  dplyr::select(Geneid,external_gene_name, description, everything()) %>% 
  mutate(Mean_all_1000 = rowMeans(.[4:6])) %>% 
  mutate(round(.[, 4:7], digit = 2)) %>% 
  arrange(-Mean_all_1000)

CPU_300_all <- cpm_CPU300 %>% 
  rownames_to_column(var = "Geneid") %>% 
  inner_join(ANN_file, by = c("Geneid" = "ensembl_gene_id")) %>%
  dplyr::select(Geneid,external_gene_name, description, everything()) %>% 
  mutate(Mean_all_300 = rowMeans(.[4:7])) %>% 
  mutate(round(.[, 4:8], digit = 2)) %>% 
  arrange(-Mean_all_300)

CPU_1000vs300_all <- CPU_1000_all %>% 
  inner_join(CPU_300_all, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name.x, description.x, Mean_all_1000, Mean_all_300, starts_with("cpm_"))

data300vs1000_t <- data300vs1000 %>%
  t() %>% 
  data.frame()

data_1000vs300_107_t <- data_1000vs300_107 %>% 
  t() %>% 
  data.frame()

Count_1000vs300 <- as.data.frame(apply(data_1000vs300_107_t, 1, function(x) length(x[x>=5])))

coverage_1000v300 <- Count_1000vs300 %>% 
  rename(Gene_number = `apply(data_1000vs300_107_t, 1, function(x) length(x[x >= 5]))`) %>% 
  rownames_to_column(var='Name') %>% 
  arrange(Name) %>% 
  mutate(Type = c("CPU300","CPU300","CPU300","CPU300","CPU1000","CPU1000","CPU1000")) %>% 
  mutate(Mean = ifelse(.$Type == "CPU300", mean(.[c(1,2,3,4),2]), mean(.[c(5,6,7),2]))) %>% 
  mutate(SEM = ifelse(.$Type == "CPU300", sd(.[c(1,2,3,4),2])/sqrt(4), sd(.[c(5,6,7),2])/sqrt(3))) %>% 
  mutate(Mean_CPU300 = mean(.[c(1,2,3,4),2])) %>% 
  mutate(Mean_CPU1000 = mean(.[c(5,6,7),2]))

coverage_1000v300$Type = factor(coverage_1000v300$Type, 
                                levels = c("CPU1000", "CPU300"))

df_coverage_1000vs300 <- coverage_1000v300 %>% 
  dplyr::select(Type, Mean, SEM) %>% 
  .[c(1,5),]

ggplot(data = df_coverage_1000vs300, 
       aes(x=Type, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Type))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Detected transcripts")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 16100), n.breaks = 10)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = coverage_1000v300, aes(x= Type, y=Gene_number), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))
ggsave("coverage_300vs1000.pdf", width = 28 , height = 35 , units = "mm")

#####Categorization of ChINsIHC and ChINsZsG#####
All_receptor_with_4_cat <- G_Protein_Coupled_Receptors_mouse %>% 
  full_join(Cytokine_Receptors_mouse) %>% 
  full_join(Nuclear_Receptors_mouse) %>% 
  full_join(Pattern_Recognition_Receptors_mouse) %>% 
  unique()

Two_list_comparison <-  All_receptor_with_4_cat %>% 
  anti_join(All_receptor_with_4_cat)

#####Receptors#####
####TOP50 receptors#####

CPU1000_receptors <- CPU_1000_all %>% 
  inner_join(All_receptor_with_4_cat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_1000) %>% 
  mutate(Position_1000 = row_number())

CPU1000_receptors_posi <- CPU1000_receptors %>% 
  dplyr::select(Geneid, Position_1000)

CPU1000_receptors_top50 <- CPU_1000_all %>%
  inner_join(All_receptor_with_4_cat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_1000) %>%
  slice_head(., n= 50)

CPU300_receptors <- CPU_300_all %>% 
  inner_join(All_receptor_with_4_cat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_300)%>% 
  mutate(Position_300 = row_number())

CPU300_receptors_posi <- CPU300_receptors %>% 
  dplyr::select(Geneid, Position_300)

CPU300_receptors_top50 <- CPU_300_all %>%
  inner_join(All_receptor_with_4_cat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_300) %>%
  slice_head(., n= 50)

CPU_1000vs300_unio_receptors_top50 <- CPU1000_receptors_top50 %>% 
  inner_join(CPU300_receptors_top50, by = "Geneid") %>% 
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  dplyr::select(Geneid,external_gene_name.x, description.x, Mean_all_1000, Mean_all_300, logFC, starts_with("cpm_"))

CPU_1000_only_receptors_top50 <- CPU1000_receptors_top50 %>% 
  anti_join(CPU300_receptors_top50, by = "Geneid")

CPU_300_only_receptors_top50 <- CPU300_receptors_top50 %>% 
  anti_join(CPU1000_receptors_top50, by = "Geneid")

CPU_1000vs300_unio_receptors_top50_full <- CPU1000_receptors_top50 %>% 
  full_join(CPU300_receptors_top50, by = "Geneid")

CPU_1000vs300_unio_receptors_top50_full_names <- CPU_1000vs300_unio_receptors_top50_full %>% 
  dplyr::select(Geneid)

CPU_1000vs300_unio_All_receptor_with_4_cat_top50 <- CPU_1000vs300_unio_receptors_top50_full_names %>% 
  inner_join(CPU_1000vs300_all, by = "Geneid") %>%
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  inner_join(CPU1000_receptors_posi, by = "Geneid") %>%
  inner_join(CPU300_receptors_posi, by = "Geneid") %>%
  rename(Gene_symbol = external_gene_name.x) %>% 
  rename(Gene_description = description.x) %>% 
  dplyr::select(Geneid, Gene_symbol, Gene_description, Mean_all_1000, Mean_all_300, logFC, Position_1000, Position_300, starts_with("cpm_"))

####Top100 receptors#####

CPU1000_receptors_top100 <- CPU_1000_all %>%
  inner_join(All_receptor_with_4_cat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_1000) %>%
  slice_head(., n= 100)

CPU1000_receptors_top100_ <- CPU_1000_all %>%
  inner_join(All_receptors, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_1000) %>%
  slice_head(., n= 100)

CPU300_receptors_top100 <- CPU_300_all %>%
  inner_join(All_receptor_with_4_cat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_300) %>%
  slice_head(., n= 100)

CPU300_receptors_top100_ <- CPU_300_all %>%
  inner_join(All_receptors, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_300) %>%
  slice_head(., n= 100)

CPU_1000vs300_unio_receptors_top100 <- CPU1000_receptors_top100 %>% 
  inner_join(CPU300_receptors_top100, by = "Geneid") %>% 
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  dplyr::select(Geneid,external_gene_name.x, description.x, Mean_all_1000, Mean_all_300, logFC, starts_with("cpm_"))

CPU_1000_only_receptors_top100 <- CPU1000_receptors_top100 %>% 
  anti_join(CPU300_receptors_top100, by = "Geneid")

CPU_300_only_receptors_top100 <- CPU300_receptors_top100 %>% 
  anti_join(CPU1000_receptors_top100, by = "Geneid")

CPU_1000vs300_unio_receptors_top100_full <- CPU1000_receptors_top100 %>% 
  full_join(CPU300_receptors_top100, by = "Geneid")

CPU_1000vs300_unio_receptors_top100_full_names <- CPU_1000vs300_unio_receptors_top100_full %>% 
  dplyr::select(Geneid)

CPU_1000vs300_unio_All_receptor_with_4_cat_top100 <- CPU_1000vs300_unio_receptors_top100_full_names %>% 
  inner_join(CPU_1000vs300_all, by = "Geneid") %>%
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  inner_join(CPU1000_receptors_posi, by = "Geneid") %>%
  inner_join(CPU300_receptors_posi, by = "Geneid") %>%
  rename(Gene_symbol = external_gene_name.x) %>% 
  rename(Gene_description = description.x) %>% 
  dplyr::select(Geneid, Gene_symbol, Gene_description, Mean_all_1000, Mean_all_300, logFC, Position_1000, Position_300, starts_with("cpm_"))

VD_1000_rec <- CPU_1000vs300_unio_All_receptor_with_4_cat_top100 %>% 
  dplyr::filter(Position_1000<101) %>% 
  dplyr::select(Gene_symbol)

VD_1000_rec_only <- CPU_1000vs300_unio_All_receptor_with_4_cat_top100 %>% 
  dplyr::filter(Position_1000>100) %>% 
  dplyr::select(Gene_symbol)

VD_300_rec <- CPU_1000vs300_unio_All_receptor_with_4_cat_top100 %>% 
  dplyr::filter(Position_300<101) %>% 
  dplyr::select(Gene_symbol)

VD_300_rec_only <- CPU_1000vs300_unio_All_receptor_with_4_cat_top100 %>% 
  dplyr::filter(Position_300>100) %>% 
  dplyr::select(Gene_symbol)

VD_1000vs300_top4 <- CPU_1000vs300_unio_All_receptor_with_4_cat_top100 %>% 
  dplyr::filter(logFC<0.1 & logFC> -0.1) %>% 
  dplyr::arrange(-Mean_all_1000) %>%
  dplyr::filter(Mean_all_1000>=100 | Mean_all_300>=100) %>% 
  dplyr::slice_head(.,n = 4) %>%
  dplyr::mutate(SEM_all_1000 = rowSds(as.matrix(.[9:11]))/sqrt(3)) %>% 
  dplyr::mutate(SEM_all_300 = rowSds(as.matrix(.[12:15]))/sqrt(4)) %>% 
  dplyr::select(Gene_symbol, starts_with("cpm"))

VD_1000vs300_rec <- VD_1000vs300_top4 %>% 
  pivot_longer(
    cols = starts_with("cpm_"),
    names_to = "name_uni",
    values_to = "points"
  ) %>%
  mutate(Name = ifelse(grepl("cpm_CPU",name_uni), paste0("CPU300"),paste0("CPU1000")))

VD_1000vs300_rec_top1 <- VD_1000vs300_rec %>% 
  .[1:7,] %>% 
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))
VD_1000vs300_rec_top11 <- VD_1000vs300_rec_top1 %>%
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

VD_1000vs300_rec_top2 <- VD_1000vs300_rec %>% 
  .[8:14,] %>%
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))
VD_1000vs300_rec_top22 <- VD_1000vs300_rec_top2 %>% 
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

VD_1000vs300_rec_top3 <- VD_1000vs300_rec %>% 
  .[15:21,] %>%
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))
VD_1000vs300_rec_top33 <- VD_1000vs300_rec_top3 %>%
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

VD_1000vs300_rec_top4 <- VD_1000vs300_rec %>% 
  .[22:28,] %>%
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))
VD_1000vs300_rec_top44 <- VD_1000vs300_rec_top4 %>% 
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

#####top1 ggplot####
ggplot(data = VD_1000vs300_rec_top11, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 800), n.breaks = 5)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = VD_1000vs300_rec_top1, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = VD_1000vs300_rec_top1$Gene_symbol[1])
ggsave("VD_1000vs300_rec_top1.pdf", width = 28 , height = 35 , units = "mm")

######top2 ggplot#######
ggplot(data = VD_1000vs300_rec_top22, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 700), n.breaks = 7)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = VD_1000vs300_rec_top2, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = VD_1000vs300_rec_top2$Gene_symbol[1])
ggsave("VD_1000vs300_rec_top2.pdf", width = 28 , height = 35 , units = "mm")

######top3 ggplot#######
ggplot(data = VD_1000vs300_rec_top33, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 500), n.breaks = 7)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = VD_1000vs300_rec_top3, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = VD_1000vs300_rec_top3$Gene_symbol[1])
ggsave("VD_1000vs300_rec_top3.pdf", width = 28 , height = 35 , units = "mm")

######top4 ggplot#######
ggplot(data = VD_1000vs300_rec_top44, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 500), n.breaks = 7)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = VD_1000vs300_rec_top4, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = VD_1000vs300_rec_top4$Gene_symbol[1])
ggsave("VD_1000vs300_rec_top4.pdf", width = 28 , height = 35 , units = "mm")


VD_1000vs300_rec <- list(A = as.character(VD_1000_rec$Gene_symbol), B = as.character(VD_300_rec$Gene_symbol))

VD_1000vs300_dst <- VD_1000_rec %>% 
  inner_join(VD_300_rec)

VD_only_1000 <- VD_1000_rec %>% 
  anti_join(VD_300_rec)

VD_only_300 <- VD_300_rec %>% 
  anti_join(VD_1000_rec)

ggVennDiagram(VD_1000vs300_rec, category.names = c("CPU1000" , "CPU300 " ), set_color = "black")+
  coord_flip() +
  scale_fill_distiller(palette = "RdBu")

#####Ion_channels#####
####TOP50 Ion_channels#####

CPU1000_Ion_channels <- CPU_1000_all %>% 
  inner_join(ion_channel_mmu, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_1000) %>% 
  mutate(Position_1000 = row_number())

CPU1000_Ion_channels_posi <- CPU1000_Ion_channels %>% 
  dplyr::select(Geneid, Position_1000)

CPU1000_Ion_channels_top50 <- CPU_1000_all %>%
  inner_join(ion_channel_mmu, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_1000) %>%
  slice_head(., n= 50)

CPU300_Ion_channels <- CPU_300_all %>% 
  inner_join(ion_channel_mmu, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_300)%>% 
  mutate(Position_300 = row_number())

CPU300_Ion_channels_posi <- CPU300_Ion_channels %>% 
  dplyr::select(Geneid, Position_300)

CPU300_Ion_channels_top50 <- CPU_300_all %>%
  inner_join(ion_channel_mmu, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_300) %>%
  slice_head(., n= 50)

CPU_1000vs300_unio_Ion_channels_top50 <- CPU1000_Ion_channels_top50 %>% 
  inner_join(CPU300_Ion_channels_top50, by = "Geneid") %>% 
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  dplyr::select(Geneid,external_gene_name.x, description.x, Mean_all_1000, Mean_all_300, logFC, starts_with("cpm_"))

CPU_1000_only_Ion_channels_top50 <- CPU1000_Ion_channels_top50 %>% 
  anti_join(CPU300_Ion_channels_top50, by = "Geneid")

CPU_300_only_Ion_channels_top50 <- CPU300_Ion_channels_top50 %>% 
  anti_join(CPU1000_Ion_channels_top50, by = "Geneid")

CPU_1000vs300_unio_Ion_channels_top50_full <- CPU1000_Ion_channels_top50 %>% 
  full_join(CPU300_Ion_channels_top50, by = "Geneid")

CPU_1000vs300_unio_Ion_channels_top50_full_names <- CPU_1000vs300_unio_Ion_channels_top50_full %>% 
  dplyr::select(Geneid)

CPU_1000vs300_unio_all_ion_channel_top50 <- CPU_1000vs300_unio_Ion_channels_top50_full_names %>% 
  inner_join(CPU_1000vs300_all, by = "Geneid") %>%
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  inner_join(CPU1000_Ion_channels_posi, by = "Geneid") %>%
  inner_join(CPU300_Ion_channels_posi, by = "Geneid") %>%
  rename(Gene_symbol = external_gene_name.x) %>% 
  rename(Gene_description = description.x) %>% 
  dplyr::select(Geneid, Gene_symbol, Gene_description, Mean_all_1000, Mean_all_300, logFC, Position_1000, Position_300, starts_with("cpm_"))

####Top100 Ion_channels#####

CPU1000_Ion_channels_top100 <- CPU_1000_all %>%
  inner_join(ion_channel_mmu, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_1000) %>%
  slice_head(., n= 100)

CPU300_Ion_channels_top100 <- CPU_300_all %>%
  inner_join(ion_channel_mmu, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_300) %>%
  slice_head(., n= 100)

CPU_1000vs300_unio_Ion_channels_top100 <- CPU1000_Ion_channels_top100 %>% 
  inner_join(CPU300_Ion_channels_top100, by = "Geneid") %>% 
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  dplyr::select(Geneid,external_gene_name.x, description.x, Mean_all_1000, Mean_all_300, logFC, starts_with("cpm_"))

CPU_1000_only_Ion_channels_top100 <- CPU1000_Ion_channels_top100 %>% 
  anti_join(CPU300_Ion_channels_top100, by = "Geneid")

CPU_300_only_Ion_channels_top100 <- CPU300_Ion_channels_top100 %>% 
  anti_join(CPU1000_Ion_channels_top100, by = "Geneid")

CPU_1000vs300_unio_Ion_channels_top100_full <- CPU1000_Ion_channels_top100 %>% 
  full_join(CPU300_Ion_channels_top100, by = "Geneid")

CPU_1000vs300_unio_Ion_channels_top100_full_names <- CPU_1000vs300_unio_Ion_channels_top100_full %>% 
  dplyr::select(Geneid)

CPU_1000vs300_unio_all_ion_channel_top100 <- CPU_1000vs300_unio_Ion_channels_top100_full_names %>% 
  inner_join(CPU_1000vs300_all, by = "Geneid") %>%
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  inner_join(CPU1000_Ion_channels_posi, by = "Geneid") %>%
  inner_join(CPU300_Ion_channels_posi, by = "Geneid") %>%
  rename(Gene_symbol = external_gene_name.x) %>% 
  rename(Gene_description = description.x) %>% 
  dplyr::select(Geneid, Gene_symbol, Gene_description, Mean_all_1000, Mean_all_300, logFC, Position_1000, Position_300, starts_with("cpm_"))

VD_1000 <- CPU_1000vs300_unio_all_ion_channel_top100 %>% 
  dplyr::filter(Position_1000<101) %>% 
  dplyr::select(Gene_symbol)

VD_300 <- CPU_1000vs300_unio_all_ion_channel_top100 %>% 
  dplyr::filter(Position_300<101) %>% 
  dplyr::select(Gene_symbol)

VD_1000vs300 <- list(A = as.character(VD_1000$Gene_symbol), B = as.character(VD_300$Gene_symbol))

VD_1000vs300_dst <- VD_1000 %>% 
  inner_join(VD_300)

VD_only_1000 <- VD_1000 %>% 
  anti_join(VD_300)

VD_only_300 <- VD_300 %>% 
  anti_join(VD_1000)

ggVennDiagram(VD_1000vs300, category.names = c("CPU1000" , "CPU300 " ), set_color = "black")+
  coord_flip() +
  scale_fill_distiller(palette = "RdBu")

ICH_1000vs300_top4 <- CPU_1000vs300_unio_all_ion_channel_top100 %>% 
  dplyr::filter(logFC<0.1 & logFC> -0.1) %>% 
  dplyr::arrange(-Mean_all_1000) %>%
  dplyr::filter(Mean_all_1000>=100 | Mean_all_300>=100) %>% 
  dplyr::slice_head(.,n = 4) %>%
  dplyr::mutate(SEM_all_1000 = rowSds(as.matrix(.[9:11]))/sqrt(3)) %>% 
  dplyr::mutate(SEM_all_300 = rowSds(as.matrix(.[12:15]))/sqrt(4)) %>% 
  dplyr::select(Gene_symbol, starts_with("cpm"))


ICH_1000vs300 <- ICH_1000vs300_top4 %>% 
  pivot_longer(
    cols = starts_with("cpm_"),
    names_to = "name_uni",
    values_to = "points"
  ) %>%
  mutate(Name = ifelse(grepl("cpm_CPU",name_uni), paste0("CPU300"),paste0("CPU1000")))


ICH_1000vs300_ICH_top1 <- ICH_1000vs300 %>% 
  .[1:7,] %>% 
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))

ICH_1000vs300_ICH_top11 <- ICH_1000vs300_ICH_top1 %>%
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

ICH_1000vs300_ICH_top2 <- ICH_1000vs300 %>% 
  .[8:14,] %>%
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))

ICH_1000vs300_ICH_top22 <- ICH_1000vs300_ICH_top2 %>% 
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

ICH_1000vs300_ICH_top3 <- ICH_1000vs300 %>% 
  .[15:21,] %>%
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))

ICH_1000vs300_ICH_top33 <- ICH_1000vs300_ICH_top3 %>%
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

ICH_1000vs300_ICH_top4 <- ICH_1000vs300 %>% 
  .[22:28,] %>%
  mutate(Mean = ifelse(.$Name == "CPU1000", mean(.$points[1:3]), mean(.$points[4:7]))) %>% 
  mutate(SEM = ifelse(.$Name == "CPU1000",  sd(.$points[1:3])/sqrt(3), sd(.$points[4:7])/sqrt(4)))

ICH_1000vs300_ICH_top44 <- ICH_1000vs300_ICH_top4 %>% 
  dplyr::select(Name, Mean, SEM) %>% 
  .[c(1,4),]

#####top1 ggplot####
ggplot(data = ICH_1000vs300_ICH_top11, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 850), n.breaks = 6)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = ICH_1000vs300_ICH_top1, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = ICH_1000vs300_ICH_top1$Gene_symbol[1])
ggsave("ICH_1000vs300_ICH_top1.pdf", width = 28 , height = 35 , units = "mm")


#####top2 ggplot####
ggplot(data = ICH_1000vs300_ICH_top22, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 800), n.breaks = 6)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = ICH_1000vs300_ICH_top2, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = ICH_1000vs300_ICH_top2$Gene_symbol[1])
ggsave("ICH_1000vs300_ICH_top2.pdf", width = 28 , height = 35 , units = "mm")

#####top3 ggplot####
ggplot(data = ICH_1000vs300_ICH_top33, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 650), n.breaks = 6)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = ICH_1000vs300_ICH_top3, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = ICH_1000vs300_ICH_top3$Gene_symbol[1])
ggsave("ICH_1000vs300_ICH_top3.pdf", width = 28 , height = 35 , units = "mm")


#####top4 ggplot####
ggplot(data = ICH_1000vs300_ICH_top44, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", width = 0.8, size = 0.1)+
  labs(x=NULL, y="Gene expression (CPM)")+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 600), n.breaks = 6)+
  scale_fill_manual(values = c( "#00539CFF", "#EEA47FFF"))+
  geom_errorbar(width=0.2, size = 0.1)+
  geom_point(data = ICH_1000vs300_ICH_top4, aes(x= Name, y=points), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.2)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold"))+ 
  theme(text = element_text(size= 5)) + theme(axis.text = element_text(size= 5),
                                              axis.text.x = element_text(colour = "black"),
                                              axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size= 5, face = "bold"))+labs(title = ICH_1000vs300_ICH_top4$Gene_symbol[1])
ggsave("ICH_1000vs300_ICH_top4.pdf", width = 28 , height = 35 , units = "mm")



#####supplementary table 1 Transcription profile#####

#####1.Excel Tab ChINsIHC and ChINsZsG transcriptome profiles#####

CPU_1000vs300_tp <- CPU_1000_all %>% 
  inner_join(CPU_300_all) %>% 
  dplyr::select(Geneid, external_gene_name, description, Mean_all_1000, Mean_all_300) %>% 
  rename(Gene_symbols = external_gene_name) %>% 
  rename(Gene_description = description) %>% 
  filter(Mean_all_1000>=1 | Mean_all_300>=1)


CPU_1000vs300_all_TP <- CPU_1000vs300_all %>% 
  dplyr::select(Geneid, external_gene_name.x, description.x, Mean_all_1000, Mean_all_300) %>% 
  rename(Gene_symbol = external_gene_name.x) %>% 
  rename(Gene_description = description.x) %>% 
  filter(Mean_all_1000>=1 | Mean_all_300>=1)

#####2.Excel Tab Top receptors of ChINsIHC and ChINsZsG#####

CPU_1000vs300_unio_All_receptor_with_4_cat_top100 <- CPU_1000vs300_unio_receptors_top100_full_names %>% 
  inner_join(CPU_1000vs300_all, by = "Geneid") %>%
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  inner_join(CPU1000_receptors_posi, by = "Geneid") %>%
  inner_join(CPU300_receptors_posi, by = "Geneid") %>%
  rename(Gene_symbol = external_gene_name.x) %>% 
  rename(Gene_description = description.x) %>% 
  dplyr::select(Geneid, Gene_symbol, Gene_description, Mean_all_1000, Mean_all_300, logFC, Position_1000, Position_300, starts_with("cpm_"))

#####3.Excel Tab Top ion channel of ChINsIHC and ChINsZsG######

CPU_1000vs300_unio_all_ion_channel_top100 <- CPU_1000vs300_unio_Ion_channels_top100_full_names %>% 
  inner_join(CPU_1000vs300_all, by = "Geneid") %>%
  mutate(logFC = log(Mean_all_1000/Mean_all_300,2)) %>% 
  inner_join(CPU1000_Ion_channels_posi, by = "Geneid") %>%
  inner_join(CPU300_Ion_channels_posi, by = "Geneid") %>%
  rename(Gene_symbol = external_gene_name.x) %>% 
  rename(Gene_description = description.x) %>% 
  dplyr::select(Geneid, Gene_symbol, Gene_description, Mean_all_1000, Mean_all_300, logFC, Position_1000, Position_300, starts_with("cpm_"))


CPU_1000vs300_all_marker_genes <- CPU_1000vs300_all %>% 
  dplyr::filter(external_gene_name.x %in% "Chat" |
                  external_gene_name.x %in% "Slc5a7"|
                  external_gene_name.x %in% "Slc18a3"|
                  external_gene_name.x %in% "Ache"|
                  external_gene_name.x %in% "Chrm2"|
                  external_gene_name.x %in% "Lhx8"|
                  external_gene_name.x %in% "Isl1"|
                  external_gene_name.x %in% "Gbx2"|
                  external_gene_name.x %in% "Slc17a8") %>%
  dplyr::rename("Gene_symbol" = external_gene_name.x) %>% 
  dplyr::select(Gene_symbol, Mean_all_300, Mean_all_1000)


###Load count matrix of Rat GnRH neuron transcriptome####

files <- list.files(pattern = '.txt')    
for(i in files) {
  x <-read_delim(i, delim = "\t", escape_double = FALSE, 
                 col_types = cols(Chr = col_skip(), Start = col_skip(), 
                                  End = col_skip(), Strand = col_skip(), 
                                  Length = col_skip()), trim_ws = TRUE, 
                 skip = 1)
  assign(i,x)  
}


Rat_iLCM<- featureCounts_GRCr39_111__GnRH.txt %>% 
  rename_at(vars(contains('Final_GRCr39_111_')), list( ~ gsub('Final_GRCr39_111_', '', .))) %>%  
  rename_at(vars(contains('_trimmed_finalAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmed_finalAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ gsub('-', '_', .)))

Rat_iLCM_wgn <- Rat_iLCM %>% 
  data.frame(row.names = 1)

Rat_iLCM_ann= getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                    filters = c('ensembl_gene_id') ,
                    values = Rat_iLCM$Geneid,
                    mart = ensembl111_rat)

Rat_iLCM_cpm = as.data.frame(cpm(Rat_iLCM_wgn))
colnames(Rat_iLCM_cpm) = paste0('cpm_', colnames(Rat_iLCM_cpm))
Rat_iLCM_cpm <- rownames_to_column(Rat_iLCM_cpm,var = "Geneid")

Rat_iLCM_all <- Rat_iLCM%>% 
  inner_join(Rat_iLCM_ann, by = c("Geneid" = "ensembl_gene_id")) %>%
  inner_join(Rat_iLCM_cpm, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name, description, TSeq10_S10, TSeq8_S8, TSeq9_S9, cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9) %>% 
  mutate(Mean_all_Rat = rowMeans(.[7:9]))

Rat_iLCM_all_p <- Rat_iLCM%>% 
  inner_join(Rat_iLCM_ann, by = c("Geneid" = "ensembl_gene_id")) %>%
  inner_join(Rat_iLCM_cpm, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name, description, cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9) %>% 
  mutate(Mean_all_Rat = rowMeans(.[4:6]))


filter3__rat_cpm <- apply(Rat_iLCM_p, 1, function(x) length(x[x>1])==3)
Rat_iLCM_p_3_3_cpm <- Rat_iLCM_p[filter3__rat_cpm,]

filter0_rat_cpm <- apply(Rat_iLCM_p, 1, function(x) length(x[x==0])==3)
Rat_iLCM_p_3_0_cpm <- Rat_iLCM_p[filter0_rat_cpm,]

Rat_iLCM_p_3_3_cpm <- rownames_to_column(Rat_iLCM_p_3_3_cpm, var='Geneid')
Rat_iLCM_p_3_0_cpm <- rownames_to_column(Rat_iLCM_p_3_0_cpm, var='Geneid')

Rat_iLCM_p_full_filtered <- Rat_iLCM_p_3_3_cpm %>% 
  bind_rows(., Rat_iLCM_p_3_0_cpm) %>%
  unique() %>% 
  mutate(Mean_all_Rat = rowMeans(.[2:4])) %>%
  arrange(-Mean_all_Rat) %>% 
  inner_join(Rat_ann) %>% 
  dplyr::select(Geneid, external_gene_name, description, everything()) 

#####Load count matrix of Mouse GnRH neuron transcriptome#####

files <- list.files(pattern = '.txt')    
for(i in files) {
  x <-read_delim(i, delim = "\t", escape_double = FALSE, 
                 col_types = cols(Chr = col_skip(), Start = col_skip(), 
                                  End = col_skip(), Strand = col_skip(), 
                                  Length = col_skip()), trim_ws = TRUE, 
                 skip = 1)
  assign(i,x)  
}

Mouse_iLCM<- featureCounts_GRCm39.111_GnRH.txt %>% 
  rename_at(vars(contains('Final_GRCm39_111_')), list( ~ gsub('Final_GRCm39_111_', '', .))) %>%  
  rename_at(vars(contains('_trimmed_finalAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmed_finalAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ gsub('-', '_', .)))

Mouse_iLCM_wgn <- Mouse_iLCM %>% 
  data.frame(row.names = 1)

Mouse_iLCM_ann= getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                      filters = c('ensembl_gene_id') ,
                      values = Mouse_iLCM$Geneid,
                      mart = ensembl111_mouse)

Mouse_iLCM_cpm = as.data.frame(cpm(Mouse_iLCM_wgn))
colnames(Mouse_iLCM_cpm) = paste0('cpm_', colnames(Mouse_iLCM_cpm))
Mouse_iLCM_cpm <- rownames_to_column(Mouse_iLCM_cpm,var = "Geneid")

Mouse_iLCM_all <- Mouse_iLCM%>% 
  inner_join(Mouse_iLCM_ann, by = c("Geneid" = "ensembl_gene_id")) %>%
  inner_join(Mouse_iLCM_cpm, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name, description, TSeq13_S13, TSeq14_S14, TSeq1_S1, cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1) %>% 
  mutate(Mean_all_mouse = rowMeans(.[7:9]))

Mouse_iLCM_all_p <- Mouse_iLCM%>% 
  inner_join(Mouse_iLCM_ann, by = c("Geneid" = "ensembl_gene_id")) %>%
  inner_join(Mouse_iLCM_cpm, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name, description, cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1) %>% 
  mutate(Mean_all_mouse = rowMeans(.[4:6]))

########Coverage of Mouse vs Rat GnRH transcriptomes######

Mouse_coverage <- Mouse_iLCM_all %>% 
  dplyr::select(Geneid, TSeq1_S1,TSeq13_S13,TSeq14_S14) %>%
  data.frame(row.names = 1) %>% 
  t() %>% 
  data.frame()

Mouse_coverage_cpm <- Mouse_iLCM_all %>% 
  dplyr::select(Geneid, cpm_TSeq1_S1,cpm_TSeq13_S13,cpm_TSeq14_S14) %>%
  data.frame(row.names = 1) %>% 
  t() %>% 
  data.frame()

Count_Mouse <- as.data.frame(apply(Mouse_coverage, 1, function(x) length(x[x>=5])))

Count_Mouse_1cpm <- as.data.frame(apply(Mouse_coverage_cpm, 1, function(x) length(x[x>=1])))

coverage_Mouse <- Count_Mouse %>% 
  rename(Gene_number = `apply(Mouse_coverage, 1, function(x) length(x[x >= 5]))`) %>% 
  rownames_to_column(var='Name') %>% 
  arrange(Name) %>% 
  mutate(Mean = mean(.[1:3,2])) %>% 
  mutate(SEM = sd(.[1:3,2]/sqrt(3)))

coverage_Mouse_1cpm <- Count_Mouse_1cpm %>% 
  rename(Gene_number = `apply(Mouse_coverage_cpm, 1, function(x) length(x[x >= 1]))`) %>% 
  rownames_to_column(var='Name') %>% 
  arrange(Name) %>% 
  mutate(Mean = mean(.[1:3,2])) %>% 
  mutate(SEM = sd(.[1:3,2]/sqrt(3)))

Rat_coverage <- Rat_iLCM_all %>% 
  dplyr::select(Geneid, TSeq8_S8, TSeq9_S9, TSeq10_S10) %>%
  data.frame(row.names = 1) %>%
  t() %>% 
  data.frame()

Rat_coverage_1cpm <- Rat_iLCM_all %>% 
  dplyr::select(Geneid, cpm_TSeq8_S8, cpm_TSeq9_S9, cpm_TSeq10_S10) %>%
  data.frame(row.names = 1) %>%
  t() %>% 
  data.frame()

Count_Rat <- as.data.frame(apply(Rat_coverage, 1, function(x) length(x[x>=5])))
Count_Rat_1cpm <- as.data.frame(apply(Rat_coverage_1cpm, 1, function(x) length(x[x>=1])))

coverage_Rat <- Count_Rat %>% 
  rename(Gene_number = `apply(Rat_coverage, 1, function(x) length(x[x >= 5]))`) %>% 
  rownames_to_column(var='Name') %>% 
  arrange(Name) %>% 
  mutate(Mean = mean(.[1:3,2])) %>% 
  mutate(SEM = sd(.[1:3,2]/sqrt(3)))

coverage_Rat_1cpm <- Count_Rat_1cpm %>% 
  rename(Gene_number = `apply(Rat_coverage_1cpm, 1, function(x) length(x[x >= 1]))`) %>% 
  rownames_to_column(var='Name') %>% 
  arrange(Name) %>% 
  mutate(Mean = mean(.[1:3,2])) %>% 
  mutate(SEM = sd(.[1:3,2]/sqrt(3)))

#####Transcriptome profile of rodent GnRH neurons####

#####Mouse GnRH neuron transcriptome#####

Mouse_ann <- Mouse_iLCM_all %>% 
  dplyr::select(Geneid, external_gene_name, description)

Mouse_iLCM_p <- Mouse_iLCM_all %>% 
  dplyr::select(Geneid, cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1) %>% 
  data.frame(row.names = 1)

filter3 <- apply(Mouse_iLCM_p, 1, function(x) length(x[x>1])==3)
Mouse_iLCM_p_3_3 <- Mouse_iLCM_p[filter3,]

filter2 <- apply(Mouse_iLCM_p, 1, function(x) length(x[x>1])==2)
Mouse_iLCM_p_3_2 <- Mouse_iLCM_p[filter2,]

filter1 <- apply(Mouse_iLCM_p, 1, function(x) length(x[x>1])==1)
Mouse_iLCM_p_3_1 <- Mouse_iLCM_p[filter1,]


Mouse_iLCM_p_3_1 <- rownames_to_column(Mouse_iLCM_p_3_1, var='Geneid')
Mouse_iLCM_p_3_2 <- rownames_to_column(Mouse_iLCM_p_3_2, var='Geneid')
Mouse_iLCM_p_3_3 <- rownames_to_column(Mouse_iLCM_p_3_3, var='Geneid')

Mouse_iLCM_p_3_3 <- Mouse_iLCM_p_3_3 %>% 
  mutate(Mean3 = rowMeans(.[2:4])) %>% 
  arrange(-Mean3)%>% 
  dplyr::select(-Mean3)

Mouse_iLCM_p_3_1 <- Mouse_iLCM_p_3_1 %>% 
  mutate(Mean1 = rowMeans(.[2:4])) %>% 
  arrange(-Mean1)%>% 
  dplyr::select(-Mean1)

Mouse_iLCM_p_3_2 <- Mouse_iLCM_p_3_2 %>% 
  mutate(Mean2 = rowMeans(.[2:4])) %>% 
  arrange(-Mean2) %>% 
  dplyr::select(-Mean2)

Mouse_iLCM_p_full <- Mouse_iLCM_p_3_3 %>% 
  bind_rows(., Mouse_iLCM_p_3_2) %>% 
  bind_rows(., Mouse_iLCM_p_3_1) %>%
  unique() %>% 
  mutate(Mean_all_mouse = rowMeans(.[2:4]))

Mouse_iLCM_p_full <- Mouse_iLCM_p_full %>% 
  inner_join(Mouse_ann, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name, description, everything() ) %>% 
  unique()

Mouse_iLCM_p_full_jav <- Mouse_iLCM_p_full %>%
  mutate(round(.[, 4:7], digit = 2)) %>% 
  rename("Ensembl ID" = Geneid) %>% 
  rename("Gene symbol" = external_gene_name) %>% 
  rename("Gene description" = description) %>% 
  rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .)))

####RAT GnRH neuron transcriptome####

Rat_ann <- Rat_iLCM_all %>% 
  dplyr::select(Geneid, external_gene_name, description)

Rat_iLCM_p <- Rat_iLCM_all %>% 
  dplyr::select(Geneid, cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9 ) %>% 
  data.frame(row.names = 1)

filter3 <- apply(Rat_iLCM_p, 1, function(x) length(x[x>1])==3)
Rat_iLCM_p_3_3 <- Rat_iLCM_p[filter3,]

filter2 <- apply(Rat_iLCM_p, 1, function(x) length(x[x>1])==2)
Rat_iLCM_p_3_2 <- Rat_iLCM_p[filter2,]

filter1 <- apply(Rat_iLCM_p, 1, function(x) length(x[x>1])==1)
Rat_iLCM_p_3_1 <- Rat_iLCM_p[filter1,]


Rat_iLCM_p_3_1 <- rownames_to_column(Rat_iLCM_p_3_1, var='Geneid')
Rat_iLCM_p_3_2 <- rownames_to_column(Rat_iLCM_p_3_2, var='Geneid')
Rat_iLCM_p_3_3 <- rownames_to_column(Rat_iLCM_p_3_3, var='Geneid')

Rat_iLCM_p_3_3 <- Rat_iLCM_p_3_3 %>% 
  mutate(Mean3 = rowMeans(.[2:4])) %>% 
  arrange(-Mean3)%>% 
  dplyr::select(-Mean3)

Rat_iLCM_p_3_1 <- Rat_iLCM_p_3_1 %>% 
  mutate(Mean1 = rowMeans(.[2:4])) %>% 
  arrange(-Mean1)%>% 
  dplyr::select(-Mean1)

Rat_iLCM_p_3_2 <- Rat_iLCM_p_3_2 %>% 
  mutate(Mean2 = rowMeans(.[2:4])) %>% 
  arrange(-Mean2) %>% 
  dplyr::select(-Mean2)

Rat_iLCM_p_full <- Rat_iLCM_p_3_3 %>% 
  bind_rows(., Rat_iLCM_p_3_2) %>% 
  bind_rows(., Rat_iLCM_p_3_1) %>%
  unique() %>% 
  mutate(Mean_all_Rat = rowMeans(.[2:4]))

Rat_iLCM_p_full <- Rat_iLCM_p_full %>% 
  inner_join(Rat_ann, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name, description, everything() ) %>% 
  unique()

Rat_iLCM_p_full_jav <- Rat_iLCM_p_full %>%
  mutate(round(.[, 4:7], digit = 2)) %>% 
  rename("Ensembl ID" = Geneid) %>% 
  rename("Gene symbol" = external_gene_name) %>% 
  rename("Gene description" = description) %>% 
  rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .)))


#####Categorization of the transcriptome of mouse GnRH neurons####

Mouse_iLCM_all_cpm_cat <- Mouse_iLCM_p_full %>%
  mutate(round(.[,4:7], digit = 2)) %>% 
  dplyr::select(Geneid, external_gene_name, description, starts_with("cpm"), Mean_all_mouse)

Neuropeptide_Mouse_gnrh <- Neuropeptides_mouse %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

TF_Mouse_gnrh <- TF_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

Transporter_Mouse_gnrh <- Transporters_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

Nuclear_Receptors_Mouse_gnrh <- Nuclear_Receptors_mouse %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

GPCR_Receptors_Mouse_gnrh <- G_Protein_Coupled_Receptors_mouse %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

Pattern_Recognition_Receptors_Mouse_gnrh <- Pattern_Recognition_Receptors_mouse %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

Cytokine_Receptors_Mouse_gnrh <- Cytokine_Receptors_mouse %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

Cell_Adhesion_Molecules_Mouse_gnrh <- Cell_Adhesion_Molecules_mouse %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

Ion_channel_Mouse_gnrh <- Ion_Channels_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_mouse >= 5)

#####Categorization of the transcriptome of Rat GnRH neurons####

Rat_iLCM_all_cpm_cat <- Rat_iLCM_p_full %>%
  mutate(round(.[,4:7], digit = 2)) %>%
  filter(!grepl("ENSRNOG00000013433", Geneid))%>% 
  dplyr::select(Geneid, external_gene_name, description, starts_with("cpm"), Mean_all_Rat)

Neuropeptide_Rat_gnrh <- Neuropeptides_Rat %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

TF_Rat_gnrh <- Transcription_Factors_Rat  %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

Transporter_Rat_gnrh <- Transporters_Rat  %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

Nuclear_Receptors_Rat_gnrh <- Nuclear_Receptors_Rat %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

GPCR_Receptors_Rat_gnrh <- G_Protein_Coupled_Receptors_Rat %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

Pattern_Recognition_Receptors_Rat_gnrh <- Pattern_Recognition_Receptors_Rat %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

Cytokine_Receptors_Rat_gnrh <- Cytokine_Receptors_Rat %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

Cell_Adhesion_Molecules_Rat_gnrh <- Cell_Adhesion_Molecules_Rat %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

Ion_channel_Rat_gnrh <- Ion_Channels_Rattus_norvegicus_rat_  %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>% 
  dplyr::filter(Mean_all_Rat >= 5)

MmuvsRno_all <- Mouse_iLCM_all_ %>% 
  inner_join(Rat_iLCM_all_, by = "external_gene_name") %>% 
  dplyr::select(Geneid.x,Geneid.y, external_gene_name, description.x, Mean_all_mouse, Mean_all_Rat, starts_with("cpm_"))

#####Neuropeptides#####
####top20 Neuropeptides#####

Neuropeptides_mmu <- Neuropeptides_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(!c(Geneid,description)) %>% 
  dplyr::select(symbol, everything()) %>%
  mutate(Position_mmu = row_number())%>% 
  unique()

Neuropeptides_mmu_posi <- Neuropeptides_mmu %>% 
  dplyr::select(symbol, Position_mmu)

Neuropeptides_mmu_top20 <- Mouse_iLCM_all_cpm_cat %>%
  inner_join(Neuropeptides_mouse, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_mouse) %>%
  slice_head(., n= 20)

Neuropeptides_rno <- Neuropeptides_Rat  %>%
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(!c(Geneid,description)) %>%
  dplyr::select( symbol, everything()) %>%
  mutate(Position_rno = row_number()) %>%
  unique()

Neuropeptides_rno_posi <- Neuropeptides_rno %>% 
  dplyr::select(symbol, Position_rno)

Neuropeptides_rno_top20 <- Rat_iLCM_all_cpm_cat %>%
  inner_join(Neuropeptides_Rat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_Rat) %>%
  slice_head(., n= 20)

# MmuvsRno_unio_Neuropeptides <- Neuropeptides_mmu %>%
#   inner_join(Neuropeptides_rno, by = "symbol") %>% 
#   mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
#   dplyr::select(symbol,Mean_all_mouse, Mean_all_Rat, logFC, starts_with("cpm_"))


MmuvsRno_unio_Neuropeptides_top20 <- Neuropeptides_mmu_top20 %>% 
  inner_join(Neuropeptides_rno_top20, by = "external_gene_name") %>% 
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  dplyr::select(Geneid.x,external_gene_name, description.x, Mean_all_mouse, Mean_all_Rat, logFC, starts_with("cpm_"))

Neuropeptides_only_mmu_top20 <- Neuropeptides_mmu_top20 %>% 
  anti_join(Neuropeptides_rno_top20, by = "external_gene_name")

Neuropeptides_only_rno_top20 <- Neuropeptides_rno_top20 %>% 
  anti_join(Neuropeptides_mmu_top20, by = "external_gene_name")

MmuvsRno_unio_Neuropeptides_top20_full <- Neuropeptides_mmu_top20 %>% 
  full_join(Neuropeptides_rno_top20, by = "external_gene_name")

MmuvsRno_unio_Neuropeptides_top20_full_names <- MmuvsRno_unio_Neuropeptides_top20_full %>% 
  dplyr::select(external_gene_name)

MmuvsRno_unio_all_Neuropeptides_top20 <- MmuvsRno_unio_Neuropeptides_top20_full_names %>% 
  inner_join(MmuvsRno_all, by = "external_gene_name") %>%
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  inner_join(Neuropeptides_mmu_posi, by = c("external_gene_name"="symbol")) %>%
  inner_join(Neuropeptides_rno_posi, by = c("external_gene_name"="symbol")) %>%
  rename(Gene_symbol = external_gene_name) %>% 
  rename(Gene_description = description.x) %>% 
  rename(Geneid_Mouse = Geneid.x) %>% 
  rename(Geneid_Rat = Geneid.y) %>% 
  dplyr::select(Geneid_Mouse,Geneid_Rat, Gene_symbol, Gene_description, Mean_all_mouse, Mean_all_Rat, logFC, 
                Position_mmu, Position_rno, starts_with("cpm_"))


####top40 Neuropeptides#####

Neuropeptides_mmu <- Neuropeptides_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_mmu = row_number())

Neuropeptides_mmu_posi <- Neuropeptides_mmu %>% 
  dplyr::select(symbol, Position_mmu)

Neuropeptides_mmu_top40 <- Mouse_iLCM_all_cpm_cat %>%
  inner_join(Neuropeptides_mouse, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_mouse) %>%
  slice_head(., n= 40)

Neuropeptides_rno <- Neuropeptides_Rat %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_rno = row_number())

Neuropeptides_rno_posi <- Neuropeptides_rno %>% 
  dplyr::select(symbol, Position_rno)

Neuropeptides_rno_top40 <- Rat_iLCM_all_cpm_cat %>%
  inner_join(Neuropeptides_Rat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_Rat) %>%
  slice_head(., n= 40)

MmuvsRno_unio_Neuropeptides_top40 <- Neuropeptides_mmu_top40 %>% 
  inner_join(Neuropeptides_rno_top40, by = "external_gene_name") %>% 
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  dplyr::select(Geneid.x,external_gene_name, description.x, Mean_all_mouse, Mean_all_Rat, logFC, starts_with("cpm_"))

Neuropeptides_only_mmu_top40 <- Neuropeptides_mmu_top40 %>% 
  anti_join(Neuropeptides_rno_top40, by = "external_gene_name")

Neuropeptides_only_rno_top40 <- Neuropeptides_rno_top40 %>% 
  anti_join(Neuropeptides_mmu_top40, by = "external_gene_name")

MmuvsRno_unio_Neuropeptides_top40_full <- Neuropeptides_mmu_top40 %>% 
  full_join(Neuropeptides_rno_top40, by = "external_gene_name")

MmuvsRno_unio_Neuropeptides_top40_full_names <- MmuvsRno_unio_Neuropeptides_top40_full %>% 
  dplyr::select(external_gene_name)

MmuvsRno_unio_all_Neuropeptides_top40 <- MmuvsRno_unio_Neuropeptides_top40_full_names %>% 
  inner_join(MmuvsRno_all, by = "external_gene_name") %>%
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  inner_join(Neuropeptides_mmu_posi, by = c("external_gene_name"="symbol")) %>%
  inner_join(Neuropeptides_rno_posi, by = c("external_gene_name"="symbol")) %>%
  rename(Gene_symbol = external_gene_name) %>% 
  rename(Gene_description = description.x) %>% 
  rename(Geneid_Mouse = Geneid.x) %>% 
  rename(Geneid_Rat = Geneid.y) %>% 
  dplyr::select(Geneid_Mouse,Geneid_Rat, Gene_symbol, Gene_description, Mean_all_mouse, Mean_all_Rat, logFC, 
                Position_mmu, Position_rno, starts_with("cpm_"))
#####Nuclear_Receptors#####
####top20 Nuclear_Receptors#####

Nuclear_Receptors_mmu <- Nuclear_Receptors_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_mmu = row_number())

Nuclear_Receptors_mmu_posi <- Nuclear_Receptors_mmu %>% 
  dplyr::select(symbol, Position_mmu)

Nuclear_Receptors_mmu_top20 <- Mouse_iLCM_all_cpm_cat %>%
  inner_join(Nuclear_Receptors_mouse, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_mouse) %>%
  slice_head(., n= 20)

Nuclear_Receptors_rno <- Nuclear_Receptors_Rat  %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_rno = row_number())

Nuclear_Receptors_rno_posi <- Nuclear_Receptors_rno %>% 
  dplyr::select(symbol, Position_rno)

Nuclear_Receptors_rno_top20 <- Rat_iLCM_all_cpm_cat %>%
  inner_join(Nuclear_Receptors_Rat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_Rat) %>%
  slice_head(., n= 20)

MmuvsRno_unio_Nuclear_Receptors_top20 <- Nuclear_Receptors_mmu_top20 %>% 
  inner_join(Nuclear_Receptors_rno_top20, by = "external_gene_name") %>% 
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  dplyr::select(Geneid.x,external_gene_name, description.x, Mean_all_mouse, Mean_all_Rat, logFC, starts_with("cpm_"))

Nuclear_Receptors_only_mmu_top20 <- Nuclear_Receptors_mmu_top20 %>% 
  anti_join(Nuclear_Receptors_rno_top20, by = "external_gene_name")

Nuclear_Receptors_only_rno_top20 <- Nuclear_Receptors_rno_top20 %>% 
  anti_join(Nuclear_Receptors_mmu_top20, by = "external_gene_name")

MmuvsRno_unio_Nuclear_Receptors_top20_full <- Nuclear_Receptors_mmu_top20 %>% 
  full_join(Nuclear_Receptors_rno_top20, by = "external_gene_name")

MmuvsRno_unio_Nuclear_Receptors_top20_full_names <- MmuvsRno_unio_Nuclear_Receptors_top20_full %>% 
  dplyr::select(external_gene_name)

MmuvsRno_unio_all_Nuclear_Receptors_top20 <- MmuvsRno_unio_Nuclear_Receptors_top20_full_names %>% 
  inner_join(MmuvsRno_all, by = "external_gene_name") %>%
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  inner_join(Nuclear_Receptors_mmu_posi, by = c("external_gene_name"="symbol")) %>%
  inner_join(Nuclear_Receptors_rno_posi, by = c("external_gene_name"="symbol")) %>%
  rename(Gene_symbol = external_gene_name) %>% 
  rename(Gene_description = description.x) %>% 
  rename(Geneid_Mouse = Geneid.x) %>% 
  rename(Geneid_Rat = Geneid.y) %>% 
  dplyr::select(Geneid_Mouse,Geneid_Rat, Gene_symbol, Gene_description, Mean_all_mouse, Mean_all_Rat, logFC, 
                Position_mmu, Position_rno, starts_with("cpm_"))


####top40 Nuclear_Receptors#####

Nuclear_Receptors_mmu <- Nuclear_Receptors_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_mmu = row_number())

Nuclear_Receptors_mmu_posi <- Nuclear_Receptors_mmu %>% 
  dplyr::select(symbol, Position_mmu)

Nuclear_Receptors_mmu_top40 <- Mouse_iLCM_all_cpm_cat %>%
  inner_join(Nuclear_Receptors_mouse, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_mouse) %>%
  slice_head(., n= 40)

Nuclear_Receptors_rno <- Nuclear_Receptors_Rat  %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_rno = row_number())

Nuclear_Receptors_rno_posi <- Nuclear_Receptors_rno %>% 
  dplyr::select(symbol, Position_rno)

Nuclear_Receptors_rno_top40 <- Rat_iLCM_all_cpm_cat %>%
  inner_join(Nuclear_Receptors_Rat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_Rat) %>%
  slice_head(., n= 40)

MmuvsRno_unio_Nuclear_Receptors_top40 <- Nuclear_Receptors_mmu_top40 %>% 
  inner_join(Nuclear_Receptors_rno_top40, by = "external_gene_name") %>% 
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  dplyr::select(Geneid.x,external_gene_name, description.x, Mean_all_mouse, Mean_all_Rat, logFC, starts_with("cpm_"))

Nuclear_Receptors_only_mmu_top40 <- Nuclear_Receptors_mmu_top40 %>% 
  anti_join(Nuclear_Receptors_rno_top40, by = "external_gene_name")

Nuclear_Receptors_only_rno_top40 <- Nuclear_Receptors_rno_top40 %>% 
  anti_join(Nuclear_Receptors_mmu_top40, by = "external_gene_name")

MmuvsRno_unio_Nuclear_Receptors_top40_full <- Nuclear_Receptors_mmu_top40 %>% 
  full_join(Nuclear_Receptors_rno_top40, by = "external_gene_name")

MmuvsRno_unio_Nuclear_Receptors_top40_full_names <- MmuvsRno_unio_Nuclear_Receptors_top40_full %>% 
  dplyr::select(external_gene_name)

MmuvsRno_unio_all_Nuclear_Receptors_top40 <- MmuvsRno_unio_Nuclear_Receptors_top40_full_names %>% 
  inner_join(MmuvsRno_all, by = "external_gene_name") %>%
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  inner_join(Nuclear_Receptors_mmu_posi, by = c("external_gene_name"="symbol")) %>%
  inner_join(Nuclear_Receptors_rno_posi, by = c("external_gene_name"="symbol")) %>%
  rename(Gene_symbol = external_gene_name) %>% 
  rename(Gene_description = description.x) %>% 
  rename(Geneid_Mouse = Geneid.x) %>% 
  rename(Geneid_Rat = Geneid.y) %>% 
  dplyr::select(Geneid_Mouse,Geneid_Rat, Gene_symbol, Gene_description, Mean_all_mouse, Mean_all_Rat, logFC, 
                Position_mmu, Position_rno, starts_with("cpm_"))


#####G_Protein_Coupled_Receptors#####
####TOP50 G_Protein_Coupled_Receptors#####

G_Protein_Coupled_Receptors_mmu <- G_Protein_Coupled_Receptors_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_mmu = row_number())

G_Protein_Coupled_Receptors_mmu_posi <- G_Protein_Coupled_Receptors_mmu %>% 
  dplyr::select(symbol, Position_mmu)

G_Protein_Coupled_Receptors_mmu_top50 <- Mouse_iLCM_all_cpm_cat %>%
  inner_join(G_Protein_Coupled_Receptors_mouse, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_mouse) %>%
  slice_head(., n= 50)

G_Protein_Coupled_Receptors_rno <- G_Protein_Coupled_Receptors_Rat  %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_rno = row_number())

G_Protein_Coupled_Receptors_rno_posi <- G_Protein_Coupled_Receptors_rno %>% 
  dplyr::select(symbol, Position_rno)

G_Protein_Coupled_Receptors_rno_top50 <- Rat_iLCM_all_cpm_cat %>%
  inner_join(G_Protein_Coupled_Receptors_Rat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_Rat) %>%
  slice_head(., n= 50)

MmuvsRno_unio_G_Protein_Coupled_Receptors_top50 <- G_Protein_Coupled_Receptors_mmu_top50 %>% 
  inner_join(G_Protein_Coupled_Receptors_rno_top50, by = "external_gene_name") %>% 
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  dplyr::select(Geneid.x,external_gene_name, description.x, Mean_all_mouse, Mean_all_Rat, logFC, starts_with("cpm_"))

G_Protein_Coupled_Receptors_only_mmu_top50 <- G_Protein_Coupled_Receptors_mmu_top50 %>% 
  anti_join(G_Protein_Coupled_Receptors_rno_top50, by = "external_gene_name")

G_Protein_Coupled_Receptors_only_rno_top50 <- G_Protein_Coupled_Receptors_rno_top50 %>% 
  anti_join(G_Protein_Coupled_Receptors_mmu_top50, by = "external_gene_name")

MmuvsRno_unio_G_Protein_Coupled_Receptors_top50_full <- G_Protein_Coupled_Receptors_mmu_top50 %>% 
  full_join(G_Protein_Coupled_Receptors_rno_top50, by = "external_gene_name")

MmuvsRno_unio_G_Protein_Coupled_Receptors_top50_full_names <- MmuvsRno_unio_G_Protein_Coupled_Receptors_top50_full %>% 
  dplyr::select(external_gene_name)

MmuvsRno_unio_all_G_Protein_Coupled_Receptors_top50 <- MmuvsRno_unio_G_Protein_Coupled_Receptors_top50_full_names %>% 
  inner_join(MmuvsRno_all, by = "external_gene_name") %>%
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  inner_join(G_Protein_Coupled_Receptors_mmu_posi, by = c("external_gene_name"="symbol")) %>%
  inner_join(G_Protein_Coupled_Receptors_rno_posi, by = c("external_gene_name"="symbol")) %>%
  rename(Gene_symbol = external_gene_name) %>% 
  rename(Gene_description = description.x) %>% 
  rename(Geneid_Mouse = Geneid.x) %>% 
  rename(Geneid_Rat = Geneid.y) %>% 
  dplyr::select(Geneid_Mouse,Geneid_Rat, Gene_symbol, Gene_description, Mean_all_mouse, Mean_all_Rat, logFC, 
                Position_mmu, Position_rno, starts_with("cpm_"))


####Top100 G_Protein_Coupled_Receptors#####

G_Protein_Coupled_Receptors_mmu <- G_Protein_Coupled_Receptors_mouse  %>% 
  inner_join(Mouse_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_mouse)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_mmu = row_number())

G_Protein_Coupled_Receptors_mmu_posi <- G_Protein_Coupled_Receptors_mmu %>% 
  dplyr::select(symbol, Position_mmu)

G_Protein_Coupled_Receptors_mmu_top100 <- Mouse_iLCM_all_cpm_cat %>%
  inner_join(G_Protein_Coupled_Receptors_mouse, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_mouse) %>%
  slice_head(., n= 100)

G_Protein_Coupled_Receptors_rno <- G_Protein_Coupled_Receptors_Rat  %>% 
  inner_join(Rat_iLCM_all_cpm_cat, by = c("symbol"="external_gene_name")) %>%
  dplyr::arrange(desc(Mean_all_Rat)) %>%
  dplyr::select(Geneid, symbol, description, everything()) %>%
  mutate(Position_rno = row_number())

G_Protein_Coupled_Receptors_rno_posi <- G_Protein_Coupled_Receptors_rno %>% 
  dplyr::select(symbol, Position_rno)

G_Protein_Coupled_Receptors_rno_top100 <- Rat_iLCM_all_cpm_cat %>%
  inner_join(G_Protein_Coupled_Receptors_Rat, by = c("external_gene_name" = "symbol")) %>%
  arrange(-Mean_all_Rat) %>%
  slice_head(., n= 100)

MmuvsRno_unio_G_Protein_Coupled_Receptors_top100 <- G_Protein_Coupled_Receptors_mmu_top100 %>% 
  inner_join(G_Protein_Coupled_Receptors_rno_top100, by = "external_gene_name") %>% 
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  dplyr::select(Geneid.x,external_gene_name, description.x, Mean_all_mouse, Mean_all_Rat, logFC, starts_with("cpm_"))

G_Protein_Coupled_Receptors_only_mmu_top100 <- G_Protein_Coupled_Receptors_mmu_top100 %>% 
  anti_join(G_Protein_Coupled_Receptors_rno_top100, by = "external_gene_name")

G_Protein_Coupled_Receptors_only_rno_top100 <- G_Protein_Coupled_Receptors_rno_top100 %>% 
  anti_join(G_Protein_Coupled_Receptors_mmu_top100, by = "external_gene_name")

MmuvsRno_unio_G_Protein_Coupled_Receptors_top100_full <- G_Protein_Coupled_Receptors_mmu_top100 %>% 
  full_join(G_Protein_Coupled_Receptors_rno_top100, by = "external_gene_name")

MmuvsRno_unio_G_Protein_Coupled_Receptors_top100_full_names <- MmuvsRno_unio_G_Protein_Coupled_Receptors_top100_full %>% 
  dplyr::select(external_gene_name)

MmuvsRno_unio_all_G_Protein_Coupled_Receptors_top100 <- MmuvsRno_unio_G_Protein_Coupled_Receptors_top100_full_names %>% 
  inner_join(MmuvsRno_all, by = "external_gene_name") %>%
  mutate(logFC = log(Mean_all_mouse/Mean_all_Rat,2)) %>% 
  inner_join(G_Protein_Coupled_Receptors_mmu_posi, by = c("external_gene_name"="symbol")) %>%
  inner_join(G_Protein_Coupled_Receptors_rno_posi, by = c("external_gene_name"="symbol")) %>%
  rename(Gene_symbol = external_gene_name) %>% 
  rename(Gene_description = description.x) %>% 
  rename(Geneid_Mouse = Geneid.x) %>% 
  rename(Geneid_Rat = Geneid.y) %>% 
  dplyr::select(Geneid_Mouse,Geneid_Rat, Gene_symbol, Gene_description, Mean_all_mouse, Mean_all_Rat, logFC, 
                Position_mmu, Position_rno, starts_with("cpm_"))



#####radar plot of Nuclear Receptors#####
scale1 = 2 
scale2 = 5 
filter = 0
filter2 = Inf 
categ = "Nuclear Receptors"
type = "abun"
cim = paste("RBC_GnRH_rat_mouse_", categ, "_sc1_", scale1, "_sc2_", scale2,
            "_filter_", filter, "_","_filter2_", filter2, "_", type, "__.pdf", sep = "")

a <- data.frame(Transcription_profiles_categorization_Mouse)
b <- data.frame(Transcription_profiles_categorization_Rat)
c <- full_join(a,b, by = c("symbol"))
c <- c %>% select(2,1,8,4:6,10:12)
c$symbol <- ifelse(is.na(c$symbol), ifelse(is.na(c$Geneid.y), c$Geneid.x, c$Geneid.y), c$symbol)
rownames(c) <- c[,1]

data <- c %>%
  select(1,4:9)
data[is.na(data)] <- 0
data <- unique(data)

data <- data %>%
  mutate(srmouse = rowMeans(data[,2:4], na.rm=TRUE)) %>% #sqrt()
  mutate(srrat = rowMeans(data[,5:7], na.rm=TRUE)) %>% #sqrt()
  mutate(log2FoldChange = log(srrat/srmouse,2))  
data[data == Inf] <- 100
data[data == -Inf] <- -100
data  <- data[rowMeans(data[2:4]) >= filter | rowMeans(data[5:7]) >= filter, ]
data <- data[rowMeans(data[2:4]) <= filter2 & rowMeans(data[5:7]) <= filter2,]

dataA <- data %>% 
  arrange(desc(srmouse)) %>% head(40)                     

dataB <- data %>%
  arrange(desc(srrat)) %>% head(40)

dataC <- semi_join(dataA,dataB)
dataA <- anti_join(dataA,dataB)
dataE <- dataA %>% filter(log2FoldChange > 0)
dataA <- anti_join(dataA,dataE)
dataB <- anti_join(dataB,dataC)
dataF <- dataB %>% filter(log2FoldChange < 0)
dataB <- anti_join(dataB,dataF)
dataD <- dataC %>% filter(log2FoldChange > 0)
dataC <- dataC %>% filter(log2FoldChange < 0)

dataA <- rbind(dataA,dataC,dataF)
dataB <- rbind(dataB,dataD,dataE)

dataA <- arrange(dataA, desc(dataA$srmouse))#desc()
dataB <- arrange(dataB, desc(dataB$srrat))

dataC <- dataA %>%
  mutate(GnRH_Mouse = srmouse-srrat) %>% 
  mutate(Overlap = srrat) %>%
  select(symbol, log2FoldChange, GnRH_Mouse, Overlap,
         cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1,
         cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9)

dataD <- dataB %>%
  mutate(GnRH_Rat = srrat-srmouse) %>% 
  mutate(Overlap = srmouse) %>%
  select(symbol, log2FoldChange, GnRH_Rat, Overlap,
         cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1,
         cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9)

dataE <- full_join(dataD,dataC) 
dataE[is.na(dataE)] <- 0
dataE <- dataE %>% mutate(CPM = GnRH_Rat+Overlap+GnRH_Mouse)

dataE$cast <- ifelse(dataE$CPM > max(dataE$CPM)/(scale1+0.1), 1,ifelse(dataE$CPM > max(dataE$CPM)/(scale2+1),2,3))

dataG <- dataE %>% filter(cast == 1)
dataH <- dataE %>% filter(cast == 2)
dataI <- dataE %>% filter(cast == 3)

dat <- dataG[1:2,]
dat[] <- 0
dat[1:2,1] <-c('z','zz')

dataG = rbind(dataG,dat)
dat[1:2,1] <-c('zzz','zzzz')
dataH = rbind(dataH,dat)
if (dataH[1,2] != 0){
    dat[1:2,1] <-c('zzzzz','zzzzzz')
  dataI = rbind(dataI,dat)
  if (dataI[1,2] != 0){
    dataE = full_join(dataG,full_join(dataH,dataI))
  } else {dataE = full_join(dataG,dataH)}
} else {dataE = dataG}

dataJ <- dataE %>% filter(grepl("Pgr|Ar|Esr1|Esr2", symbol))
dataE <- dataE %>% filter(!grepl("Pgr|Ar|Esr1|Esr2", symbol))
dataJ <- dataJ %>% arrange(factor(symbol, c('Esr2', 'Pgr', 'Ar', 'Esr1')))

dat <- dataJ[1:2,]
dat[] <- 0
dat[1:2,1] <-c('zzzzzzz','zzzzzzzz')
dataJ = rbind(dataJ,dat)
dataJ$cast <- 3
if (dataJ[1,2] != 0){
  dataE = full_join(dataJ, dataE)
}

dataE <- dataE %>% mutate(Position = row_number())

dat <- dataE %>%
  pivot_longer(cols = c(cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1,
                        cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9),
               names_to = "Cell_type",
               values_to = "sqrt_CPM")

dat$Cell_type = gsub('13|14|1' ,'' , matrix(unlist(strsplit(dat$Cell_type, '_')), nc=3, byrow=T)[,2])
dat$Cell_type[dat$Cell_type == "TSeq"] <- "GnRH_Mouse"
dat$Cell_type = gsub('0|8|9' ,'' , dat$Cell_type) 
dat$Cell_type[dat$Cell_type == "TSeq"] <- "GnRH_Rat"

dataE$GnRH_Mouse <- ifelse(dataE$cast == 2, dataE$GnRH_Mouse * scale1, dataE$GnRH_Mouse)
dataE$GnRH_Rat <- ifelse(dataE$cast == 2, dataE$GnRH_Rat * scale1, dataE$GnRH_Rat)
dataE$Overlap <- ifelse(dataE$cast == 2, dataE$Overlap * scale1, dataE$Overlap)
dataE$GnRH_Mouse <- ifelse(dataE$cast == 3, dataE$GnRH_Mouse * scale2, dataE$GnRH_Mouse)
dataE$GnRH_Rat <- ifelse(dataE$cast == 3, dataE$GnRH_Rat * scale2, dataE$GnRH_Rat)
dataE$Overlap <- ifelse(dataE$cast == 3, dataE$Overlap * scale2, dataE$Overlap)

dat$sqrt_CPM <- ifelse(dat$cast == 2, dat$sqrt_CPM * scale1, dat$sqrt_CPM)
dat$sqrt_CPM <- ifelse(dat$cast == 3, dat$sqrt_CPM * scale2, dat$sqrt_CPM)

dataE <- dataE %>%
  pivot_longer(cols = c(GnRH_Rat, GnRH_Mouse, Overlap),
               names_to = "Cell_type",
               values_to = "sqrt_CPM")

dataE$Cell_type <- factor(dataE$Cell_type, levels = c("GnRH_Rat", "GnRH_Mouse", "Overlap"))

label_data <- dataE %>% group_by(Position, symbol) %>% summarize(tot=sum(sqrt_CPM))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$Position-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)


maX <- max(ceiling(label_data$tot/10)*10)
miN <- -max(ceiling(label_data$tot/10)*4)
dat <-dat %>% mutate(cast = paste(Cell_type, cast))



ggplot(data = dataE, aes(reorder(symbol,Position), sqrt_CPM, fill = Cell_type))+
  geom_bar(stat="identity", color="black", position="stack",width=0.8,size=0.1)+
  coord_polar(theta = "x",start=0) +#start=4.42) +
  ylim(c(miN,maX))+
  
  geom_jitter(data = dat, aes(reorder(symbol,Position), sqrt_CPM, shape = Cell_type, colour = cast), position = position_dodge2(width = 0.5),
              show.legend = FALSE,  size = 0.5, stroke = 0.1)+
  theme_classic()+
  scale_fill_manual(values = c("#f1f3ff","#ef0079", "#f086c3"))+
  theme(panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray80"),
        axis.line = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 5, angle = label_data$angle, colour="black"), 
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA),
        legend.text = element_text(size = 6))+
  labs(fill = NULL)

ggsave(cim, width = 100 , height = 80 , units = "mm")

#####Radar plot of Neuropeptides#####
scale1 = 10
scale2 = 50
scale3 = 250
filter = 0
filter2 = Inf # Inf
categ = "Neuropeptide"
type = "abun"
cim = paste("RBC_GnRH_rat_mouse_", categ, "_sc1_", scale1, "_sc2_", scale2,
            "_filter_", filter, "_","_filter2_", filter2, "_", type, "__.pdf", sep = "")
a <- data.frame(Transcription_profiles_categorization_Mouse_NP)
b <- data.frame(Transcription_profiles_categorization_Rat_NP)
c <- full_join(a,b, by = c("symbol"))
c <- c %>% select(2,1,8,4:6,10:12)
c$symbol <- ifelse(is.na(c$symbol), ifelse(is.na(c$Geneid.y), c$Geneid.x, c$Geneid.y), c$symbol)
rownames(c) <- c[,1]

data <- c %>%
  select(1,4:9)
data[is.na(data)] <- 0
data <- unique(data)
data <- data %>%
  mutate(srmouse = rowMeans(data[,2:4], na.rm=TRUE)) %>% #sqrt()
  mutate(srrat = rowMeans(data[,5:7], na.rm=TRUE)) %>% #sqrt()
  mutate(log2FoldChange = log(srrat/srmouse,2))  
data[data == Inf] <- 100
data[data == -Inf] <- -100

data  <- data[rowMeans(data[2:4]) >= filter | rowMeans(data[5:7]) >= filter, ]
data <- data[rowMeans(data[2:4]) <= filter2 & rowMeans(data[5:7]) <= filter2,]

dataA <- data %>% 
  arrange(desc(srmouse)) %>% head(40)

dataB <- data %>%
  arrange(desc(srrat)) %>% head(40)
dataC <- semi_join(dataA,dataB)
dataA <- anti_join(dataA,dataB)
dataE <- dataA %>% filter(log2FoldChange > 0)
dataA <- anti_join(dataA,dataE)
dataB <- anti_join(dataB,dataC)
dataF <- dataB %>% filter(log2FoldChange < 0)
dataB <- anti_join(dataB,dataF)
dataD <- dataC %>% filter(log2FoldChange > 0)
dataC <- dataC %>% filter(log2FoldChange < 0)

dataA <- rbind(dataA,dataC,dataF)
dataB <- rbind(dataB,dataD,dataE)

dataA <- arrange(dataA, desc(dataA$srmouse))#desc()
dataB <- arrange(dataB, desc(dataB$srrat))

dataC <- dataA %>%
  mutate(GnRH_Mouse = srmouse-srrat) %>% 
  mutate(Overlap = srrat) %>%
  select(symbol, log2FoldChange, GnRH_Mouse, Overlap,
         cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1,
         cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9)

dataD <- dataB %>%
  mutate(GnRH_Rat = srrat-srmouse) %>% 
  mutate(Overlap = srmouse) %>%
  select(symbol, log2FoldChange, GnRH_Rat, Overlap,
         cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1,
         cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9)

dataE <- full_join(dataD,dataC) 
dataE[is.na(dataE)] <- 0
dataE <- dataE %>% mutate(CPM = GnRH_Rat+Overlap+GnRH_Mouse)

dataE$cast <- ifelse(dataE$CPM > max(dataE$CPM)/(scale1+0.1), 1,
                     ifelse(dataE$CPM > max(dataE$CPM)/(scale2+1),2,
                            ifelse(dataE$CPM > max(dataE$CPM)/(scale3+10),3,4)))

dataG <- dataE %>% filter(cast == 1)
dataH <- dataE %>% filter(cast == 2)
dataI <- dataE %>% filter(cast == 3)
dataK <- dataE %>% filter(cast == 4)

dat <- dataG[1:2,]
dat[] <- 0
dat[1:2,1] <-c('z','zz')

dataG = rbind(dataG,dat)
dat[1:2,1] <-c('zzz','zzzz')
dataH = rbind(dataH,dat)
if (dataH[1,2] != 0){
  dat[1:2,1] <-c('zzzzz','zzzzzz')
  dataI = rbind(dataI,dat)
  if (dataI[1,2] != 0){
    dat[1:2,1] <-c('zzzzzzzzz','zzzzzzzzzz')
    dataK = rbind(dataK,dat)
    if (dataK[1,2] !=0){
      dataE = full_join(dataG,full_join(dataH,full_join(dataI,dataK)))
    } else{ dataE = full_join(dataG,full_join(dataH,dataI))}
  } else {dataE = full_join(dataG,dataH)}
} else {dataE = dataG}

dataE <- dataE %>% mutate(Position = row_number())

dat <- dataE %>%
  pivot_longer(cols = c(cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1,
                        cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9),
               names_to = "Cell_type",
               values_to = "sqrt_CPM")

dat$Cell_type = gsub('13|14|1' ,'' , matrix(unlist(strsplit(dat$Cell_type, '_')), nc=3, byrow=T)[,2])
dat$Cell_type[dat$Cell_type == "TSeq"] <- "GnRH_Mouse"
dat$Cell_type = gsub('0|8|9' ,'' , dat$Cell_type)
dat$Cell_type[dat$Cell_type == "TSeq"] <- "GnRH_Rat"

dataE$GnRH_Mouse <- ifelse(dataE$cast == 2, dataE$GnRH_Mouse * scale1, dataE$GnRH_Mouse)
dataE$GnRH_Rat <- ifelse(dataE$cast == 2, dataE$GnRH_Rat * scale1, dataE$GnRH_Rat)
dataE$Overlap <- ifelse(dataE$cast == 2, dataE$Overlap * scale1, dataE$Overlap)
dataE$GnRH_Mouse <- ifelse(dataE$cast == 3, dataE$GnRH_Mouse * scale2, dataE$GnRH_Mouse)
dataE$GnRH_Rat <- ifelse(dataE$cast == 3, dataE$GnRH_Rat * scale2, dataE$GnRH_Rat)
dataE$Overlap <- ifelse(dataE$cast == 3, dataE$Overlap * scale2, dataE$Overlap)
dataE$GnRH_Mouse <- ifelse(dataE$cast == 4, dataE$GnRH_Mouse * scale3, dataE$GnRH_Mouse)
dataE$GnRH_Rat <- ifelse(dataE$cast == 4, dataE$GnRH_Rat * scale3, dataE$GnRH_Rat)
dataE$Overlap <- ifelse(dataE$cast == 4, dataE$Overlap * scale3, dataE$Overlap)

dat$sqrt_CPM <- ifelse(dat$cast == 2, dat$sqrt_CPM * scale1, dat$sqrt_CPM)
dat$sqrt_CPM <- ifelse(dat$cast == 3, dat$sqrt_CPM * scale2, dat$sqrt_CPM)
dat$sqrt_CPM <- ifelse(dat$cast == 4, dat$sqrt_CPM * scale3, dat$sqrt_CPM)
dataE <- dataE %>%
  pivot_longer(cols = c(GnRH_Rat, GnRH_Mouse, Overlap),
               names_to = "Cell_type",
               values_to = "sqrt_CPM")

dataE$Cell_type <- factor(dataE$Cell_type, levels = c("GnRH_Rat", "GnRH_Mouse", "Overlap"))

label_data <- dataE %>% group_by(Position, symbol) %>% summarize(tot=sum(sqrt_CPM))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$Position-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
maX <- max(ceiling(dat$sqrt_CPM/100)*100)
miN <- -max(ceiling(dat$sqrt_CPM/100)*40)
dat <-dat %>% mutate(cast = paste(Cell_type, cast))

ggplot(data = dataE, aes(reorder(symbol,Position), sqrt_CPM, fill = Cell_type))+
  geom_bar(stat="identity", color="black", position="stack",width=0.8,size=0.1)+
  coord_polar(theta = "x",start=0) +
  ylim(c(miN,maX))+
  geom_jitter(data = dat, aes(reorder(symbol,Position), sqrt_CPM, shape = Cell_type, colour = cast), position = position_dodge2(width = 0.5),
              show.legend = FALSE,  size = 0.5, stroke = 0.1)+
  theme_classic()+
  scale_fill_manual(values = c("#f1f3ff","#ef0079", "#f086c3"))+
  theme(panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray80"),
        axis.line = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 5, angle = label_data$angle, colour="black"),#size = 16, 
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA),
        legend.text = element_text(size = 6))+
  labs(fill = NULL)
ggsave(cim, width = 100 , height = 80 , units = "mm")

#######Analysis of disease databases####
reproductive_system_disease_data_mouse <- Mouse_iLCM_all_cpm_cat %>% 
  inner_join(reproductive_system_disease, by = c( "external_gene_name" = "symbol")) %>% 
  filter(Mean_all_mouse>5)

endocrine_system_disease_data_mouse <- Mouse_iLCM_all_cpm_cat %>% 
  inner_join(endocrine_system_disease, by = c( "external_gene_name" = "symbol")) %>% 
  filter(Mean_all_mouse>5)

HHG_genes_OMIM_mouse <- Mouse_iLCM_all_cpm_cat %>% 
  inner_join(OMIM_HHG_genes_database, by = c( "external_gene_name" = "Symbol")) %>% 
  filter(Mean_all_mouse>5)

HHG_genes_OMIM_Rat <- Rat_iLCM_all_cpm_cat %>% 
  inner_join(OMIM_HHG_genes_database, by = c( "external_gene_name" = "Symbol")) %>% 
  filter(Mean_all_Rat>5)

endocrine_system_disease_data_filtered_mouse <- endocrine_system_disease_data %>%
  rename(Disease_Term = 'Disease Term') %>%
  arrange(-Mean_all_mouse) %>% 
  filter(!grepl("mellitus",Disease_Term),
         !grepl("insulin",Disease_Term),
         !grepl("thyroi",Disease_Term),
         !grepl("pancre",Disease_Term),
         !grepl("tumor",Disease_Term),
         !grepl("glucoco",Disease_Term),
         !grepl("adenoma",Disease_Term),
         !grepl("goiter",Disease_Term),
         !grepl("adrenocortical",Disease_Term),
         !grepl("cortisone",Disease_Term) )

hypogonadotropic_hypogonadism_genes_in_GnRH <- endocrine_system_disease_data_filtered_mouse %>% 
  filter(grepl("hypogonadism",Disease_Term) |
           grepl("Kallmann",Disease_Term))

Pub_disease_genes_in_GnRH <- endocrine_system_disease_data_filtered_mouse %>% 
  filter(grepl("puberty",Disease_Term))

Pub_disease_genes_in_GnRH <- endocrine_system_disease_data_filtered_mouse %>% 
  filter(grepl("puberty",Disease_Term))

All_diseases_without_HHG <- endocrine_system_disease_data_filtered_mouse %>%
  filter(!grepl("hypogonadism",Disease_Term),
         !grepl("Kallmann",Disease_Term))

Graves_disease_genes_in_GnRH <- endocrine_system_disease_data_filtered_mouse %>%
  filter(grepl("Graves",Disease_Term))

All_diseases_without_HHG_Graves <- endocrine_system_disease_data_filtered_mouse %>%
  filter(!grepl("hypogonadism",Disease_Term),
         !grepl("Kallmann",Disease_Term),
         !grepl("Graves",Disease_Term))

reproductive_system_disease_data_filtered_mouse <- reproductive_system_disease_data %>%
  rename(Disease_Term = 'Disease Term') %>%
  arrange(-Mean_all_mouse) %>% 
  filter(!grepl("mellitus",Disease_Term),
         !grepl("insulin",Disease_Term),
         !grepl("thyroi",Disease_Term),
         !grepl("pancre",Disease_Term),
         !grepl("tumor",Disease_Term),
         !grepl("glucoco",Disease_Term),
         !grepl("adenoma",Disease_Term),
         !grepl("goiter",Disease_Term),
         !grepl("adrenocortical",Disease_Term))

reproductive_system_disease_data_Rat <- Rat_iLCM_all_cpm_cat %>% 
  inner_join(reproductive_system_disease, by = c( "external_gene_name" = "symbol")) %>% 
  filter(Mean_all_Rat>5)

endocrine_system_disease_data_Rat <- Rat_iLCM_all_cpm_cat %>% 
  inner_join(endocrine_system_disease, by = c( "external_gene_name" = "symbol")) %>% 
  filter(Mean_all_Rat>5)

hypogonadotropic_hypogonadism_genes_in_GnRH_Rat <- endocrine_system_disease_data_Rat %>%
  rename(Disease_Term = 'Disease Term') %>%
  arrange(-Mean_all_Rat) %>% 
  filter(grepl("hypogonadism",Disease_Term) |
           grepl("Kallmann",Disease_Term))

Just_Mouse_HHG <- HHG_genes_OMIM_mouse %>% 
  anti_join(HHG_genes_OMIM_Rat, by = "external_gene_name")

Just_Rat_HHG <- HHG_genes_OMIM_Rat %>% 
  anti_join(HHG_genes_OMIM_mouse, by = "external_gene_name")


####Transcriptome of the GnRH neuron in mixed cycle female mice ####
Mouse_iLCM_female<- featureCounts_GRCm39.111_test_06_03.txt %>% 
  rename_at(vars(contains('Final_GRCm39_111_')), list( ~ gsub('Final_GRCm39_111_', '', .))) %>%  
  rename_at(vars(contains('_trimmed_finalAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmed_finalAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ gsub('-', '_', .)))

Mouse_iLCM_female_ann= getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                             filters = c('ensembl_gene_id') ,
                             values = Mouse_iLCM_female$Geneid,
                             mart = ensembl111_mouse)

Mouse_iLCM_female_wgn <- Mouse_iLCM_female %>% 
  data.frame(row.names = 1)

Mouse_iLCM_female_cpm = as.data.frame(cpm(Mouse_iLCM_female_wgn))
colnames(Mouse_iLCM_female_cpm) = paste0('cpm_', colnames(Mouse_iLCM_female_cpm))
Mouse_iLCM_female_cpm <- rownames_to_column(Mouse_iLCM_female_cpm,var = "Geneid")

Mouse_iLCM_female_all <- Mouse_iLCM_female%>% 
  inner_join(Mouse_iLCM_female_ann, by = c("Geneid" = "ensembl_gene_id")) %>%
  inner_join(Mouse_iLCM_female_cpm, by = "Geneid") %>% 
  dplyr::select(Geneid, external_gene_name, description, everything())

######Supp. table 2######

######1.Excel Tab Transcriptome of Mouse GnRH neuron ######

Mouse_iLCM_p_full_1cpm <- Mouse_iLCM_p_full %>%
  mutate(round(.[, 4:7], digit = 2)) %>% 
  rename("Ensembl ID" = Geneid) %>% 
  rename("Gene symbol" = external_gene_name) %>% 
  rename("Gene description" = description) %>% 
  rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .))) %>% 
  filter(Mean_all_mouse>=1)

######2.Excel Tab Transcriptome of Rat GnRH neuron######

Rat_iLCM_p_full_1cpm <- Rat_iLCM_p_full %>%
  mutate(round(.[, 4:7], digit = 2)) %>% 
  rename("Ensembl ID" = Geneid) %>% 
  rename("Gene symbol" = external_gene_name) %>% 
  rename("Gene description" = description) %>% 
  rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .))) %>% 
  filter(Mean_all_Rat>=1)

######3.Excel Tab Neuropeptides #########

Neuropeptide_MvsR_symb <- Neuropeptide_Mouse_gnrh %>% 
  full_join(Neuropeptide_Rat_gnrh, by = "symbol") %>% 
  dplyr::select(symbol)

Neuropeptide_MvsR <- Neuropeptide_MvsR_symb %>% 
  inner_join(Mouse_iLCM_all, by = c("symbol" = "external_gene_name")) %>% 
  inner_join(Rat_iLCM_all, by = c("symbol" = "external_gene_name")) %>%
  filter(!Geneid.y == "ENSRNOG00000013433") %>% 
  dplyr::select(Geneid.x, Geneid.y, symbol, description.x, starts_with("cpm"), Mean_all_mouse, Mean_all_Rat) %>% 
  rename(Ensembl_ID_Mouse = Geneid.x) %>% 
  rename(Ensembl_ID_Rat = Geneid.y) %>%
  rename(Gene_symbol = symbol) %>%
  rename(Mean_Mouse = Mean_all_mouse) %>%
  rename(Mean_Rat = Mean_all_Rat) %>%
  rename(Gene_description = description.x) %>%
  dplyr::select(Ensembl_ID_Mouse, Ensembl_ID_Rat, Gene_symbol, Gene_description, cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1, 
                Mean_Mouse, cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9, Mean_Rat) %>% 
  arrange(-Mean_Mouse)

######4.Excel Tab Nuclear receptors #####

Nuclear_Receptors_MvsR_symb <- Nuclear_Receptors_Mouse_gnrh %>% 
  full_join(Nuclear_Receptors_Rat_gnrh, by = "symbol") %>% 
  dplyr::select(symbol)

Nuclear_Receptors_MvsR <- Nuclear_Receptors_MvsR_symb %>% 
  inner_join(Mouse_iLCM_all, by = c("symbol" = "external_gene_name")) %>% 
  inner_join(Rat_iLCM_all, by = c("symbol" = "external_gene_name")) %>% 
  dplyr::select(Geneid.x, Geneid.y, symbol, description.x, starts_with("cpm"), Mean_all_mouse, Mean_all_Rat) %>% 
  rename(Ensembl_ID_Mouse = Geneid.x) %>% 
  rename(Ensembl_ID_Rat = Geneid.y) %>%
  rename(Gene_symbol = symbol) %>%
  rename(Mean_Mouse = Mean_all_mouse) %>%
  rename(Mean_Rat = Mean_all_Rat) %>%
  rename(Gene_description = description.x) %>%
  dplyr::select(Ensembl_ID_Mouse, Ensembl_ID_Rat, Gene_symbol, Gene_description, cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1, 
                Mean_Mouse, cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9, Mean_Rat) %>% 
  arrange(-Mean_Mouse)

####### 5. Excel Tab GPCR receptors ######

GPCR_Receptors_MvsR_symb <- GPCR_Receptors_Mouse_gnrh %>% 
  full_join(GPCR_Receptors_Rat_gnrh, by = "symbol") %>% 
  dplyr::select(symbol)

GPCR_Receptors_M <- GPCR_Receptors_MvsR_symb %>% 
  inner_join(Mouse_iLCM_all, by = c("symbol" = "external_gene_name"))

Rat_iLCM_all_1 <- Rat_iLCM_all %>% 
  dplyr::filter(!Geneid == "ENSRNOG00000063253") %>% 
  dplyr::filter(!Geneid == "ENSRNOG00000063141")

GPCR_Receptors_R <- GPCR_Receptors_MvsR_symb %>% 
  inner_join(Rat_iLCM_all_1, by = c("symbol" = "external_gene_name"))

GPCR_Receptors_MvsR <- GPCR_Receptors_M %>% 
  full_join(GPCR_Receptors_R, by = "symbol") %>%
  dplyr::select(Geneid.x, Geneid.y, symbol, description.x, starts_with("cpm"), Mean_all_mouse, Mean_all_Rat) %>%
  rename(Ensembl_ID_Mouse = Geneid.x) %>% 
  rename(Ensembl_ID_Rat = Geneid.y) %>%
  rename(Gene_symbol = symbol) %>%
  rename(Mean_Mouse = Mean_all_mouse) %>%
  rename(Mean_Rat = Mean_all_Rat) %>%
  rename(Gene_description = description.x) %>%
  dplyr::select(Ensembl_ID_Mouse, Ensembl_ID_Rat, Gene_symbol, Gene_description, cpm_TSeq13_S13, cpm_TSeq14_S14, cpm_TSeq1_S1, 
                Mean_Mouse, cpm_TSeq10_S10, cpm_TSeq8_S8, cpm_TSeq9_S9, Mean_Rat) %>% 
  arrange(-Mean_Mouse)

mivan <- GPCR_Receptors_MvsR_symb %>% 
  full_join(GPCR_Receptors_MvsR, by = c("symbol" = "Gene_symbol"))

####### 6. Excel Tab Transcriptome of the GnRH neuron in mixed cycle female mice #######

Mouse_female_gnrh <- Mouse_iLCM_female_all %>% 
  dplyr::select(Geneid, external_gene_name, description, cpm_GnRH_S2) %>% 
  arrange(-cpm_GnRH_S2) %>%
  filter(cpm_GnRH_S2>=1) %>% 
  rename(Ensembl_ID = Geneid) %>%
  rename(Gene_symbol = external_gene_name) %>%
  rename(Gene_description = description)
