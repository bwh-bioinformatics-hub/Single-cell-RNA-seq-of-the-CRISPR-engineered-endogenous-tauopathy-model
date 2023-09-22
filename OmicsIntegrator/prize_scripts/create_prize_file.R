library(dplyr)
# Glia
Glia <- read.csv('/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/csv/Glia.csv')

colnames(Glia)
#idx <- is.na(Glia$geneSymbol)
#grep("NA",Glia$geneSymbol)
#allgenes <- Glia$geneID
#Glia$geneID[2]
#allgenes[idx]
#Glia$geneSymbol[idx]
prize <- Glia$avg_log2FC
prize <- as.data.frame(prize)
colnames(prize) <- "prize"
rownames(prize) <- Glia$geneSymbol
prize_df <- as.data.frame(cbind(prize,Glia$p_val,Glia$avg_log2FC,Glia$pct.1,Glia$pct.2,Glia$p_val_adj))

colnames(prize_df) <- c("prize","p_val","avg_log2FC","pct.1","pct.2","p_val_adj")

Glia_prize_df <- prize_df
Glia_prize_df$prize <- abs(Glia_prize_df$prize)
Glia_prize_df <- Glia_prize_df %>% filter(p_val_adj <= 0.05)
rownames(Glia_prize_df) <- Glia$geneSymbol
write.table(prize_df, file="/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/prize_scripts/Glia_prize.tsv", 
            append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)

# Optic Lobe 
OL <- read.csv('/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/csv/Optic_lobe.csv')

colnames(OL)
#idx <- is.na(Glia$geneSymbol)
#grep("NA",Glia$geneSymbol)
#allgenes <- Glia$geneID
#Glia$geneID[2]
#allgenes[idx]
#Glia$geneSymbol[idx]
prize <- abs(OL$avg_log2FC)
prize <- as.data.frame(prize)
colnames(prize) <- "prize"
OL$geneSymbol[17] <- OL$geneID[17]
OL$geneSymbol[18] <- OL$geneID[18]
OL$geneSymbol[19] <- OL$geneID[19]
rownames(prize) <- OL$geneSymbol

prize_df <- as.data.frame(cbind(prize,OL$p_val,OL$avg_log2FC,OL$pct.1,OL$pct.2,OL$p_val_adj))

colnames(prize_df) <- c("prize","p_val","avg_log2FC","pct.1","pct.2","p_val_adj")
#rownames(prize_df) <- Glia$geneSymbol
write.table(prize_df, file="/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/scripts/prize_scripts/OL_prize.tsv", 
            append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)

# Central Body 
CB <- read.csv('/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/csv/Central_Body.csv')
CB <- CB %>% filter(p_val_adj <= 0.05)
colnames(CB)
#idx <- is.na(Glia$geneSymbol)
#grep("NA",Glia$geneSymbol)
#allgenes <- Glia$geneID
#Glia$geneID[2]
#allgenes[idx]
#Glia$geneSymbol[idx]
prize <- abs(CB$avg_log2FC)
prize <- as.data.frame(prize)
colnames(prize) <- "prize"
CB$geneSymbol[1] <- CB$geneID[1]
CB$geneSymbol[21] <- CB$geneID[21]
CB$geneSymbol[606] <- CB$geneID[606]
rownames(prize) <- CB$geneSymbol

prize_df <- as.data.frame(cbind(prize,CB$p_val,CB$avg_log2FC,CB$pct.1,CB$pct.2,CB$p_val_adj))

colnames(prize_df) <- c("prize","p_val","avg_log2FC","pct.1","pct.2","p_val_adj")
#rownames(prize_df) <- Glia$geneSymbol
write.table(prize_df, file="/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/scripts/prize_scripts/CB_prize.tsv", 
            append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)

# Uannota 
CB <- read.csv('/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/scripts/prize_scripts/Central_Body.csv')
CB <- CB %>% filter(p_val_adj <= 0.05)
colnames(CB)
#idx <- is.na(Glia$geneSymbol)
#grep("NA",Glia$geneSymbol)
#allgenes <- Glia$geneID
#Glia$geneID[2]
#allgenes[idx]
#Glia$geneSymbol[idx]
prize <- abs(CB$avg_log2FC)
prize <- as.data.frame(prize)
colnames(prize) <- "prize"
CB$geneSymbol[1] <- CB$geneID[1]
CB$geneSymbol[21] <- CB$geneID[21]
CB$geneSymbol[606] <- CB$geneID[606]
rownames(prize) <- CB$geneSymbol

prize_df <- as.data.frame(cbind(prize,CB$p_val,CB$avg_log2FC,CB$pct.1,CB$pct.2,CB$p_val_adj))

colnames(prize_df) <- c("prize","p_val","avg_log2FC","pct.1","pct.2","p_val_adj")
#rownames(prize_df) <- Glia$geneSymbol
write.table(prize_df, file="/mnt/data0/projects/biohub/hassan2022/output/omicsintegrator/scripts/prize_scripts/CB_prize.tsv", 
            append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
