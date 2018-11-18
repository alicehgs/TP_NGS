path = "/Users/ahugues/desktop/tp_ngs_alice/part2/results/"
setwd(path)

library(ggplot2)
library(dplyr)
library(ggpubr)

# Load data
ReadData = function(Fichier){
  data = read.csv(Fichier, h=T, sep=";") # Read the file
  # Extract metadata from the file name 
  dataCutOff = as.numeric(strsplit(strsplit(Fichier, split="_")[[1]][[3]], ".csv")[[1]])
  
  # Add metadata to the data frame
  data$CutOff = dataCutOff
  return(data)
}

ListFiles = list.files(path = ".", pattern=".csv")
data = lapply(ListFiles, ReadData)
data = do.call(rbind, data)

threshold = 0.05

data = data %>%
  mutate(f = CutOff/(2*NbrIndividuals))

data$Epistasis = ifelse(data$Pvalue > 1-threshold/2, "Antagonistic",
                     ifelse(data$Pvalue < threshold/2, "Synergistic",
                            "None"))

nb_fp_null = length(data$Epistasis)*threshold
nb_positifs = length(data$Epistasis)-sum(ifelse(data$Epistasis == "None", 1, 0))
taux_antagoniste = sum(ifelse(data$Epistasis == "Antagonistic", 1, 0))/length(data$Epistasis)
taux_synergique = sum(ifelse(data$Epistasis == "Synergistic", 1, 0))/length(data$Epistasis)
taux_fp = nb_fp_null/nb_positifs

col = c("navyblue", "darkgoldenrod1", "springgreen3")

raw = ggplot(data, aes(x=as.factor(CutOff), y=Population, color = Epistasis)) + 
  geom_point(fill = NA, size = 3) + theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 0), legend.title=element_text(face="bold")) +
  xlab("Cut-off") + scale_size_discrete(name = "Cut-off") + scale_color_manual(values = col)
ggsave("raw.pdf", plot = raw, width=9, height=10, units="cm", dpi=500)
ggsave("raw.png", plot = raw, width=9, height=10, units="cm", dpi=500)


ggplot(data) + geom_bar(aes(x = as.factor(CutOff), fill = as.factor(Epistasis))) + xlab("Cut-off") +
  scale_fill_manual(values=col, name = "Epistasis") + theme_linedraw() + theme(strip.text.x = element_text(size=10, face="bold", color = "black"),
                                                                               strip.background = element_rect(colour="white", fill="white"), 
                                                                               legend.title=element_text(face="bold")) 
ggsave("prop.png", plot = last_plot(), width=8, height=6, units="cm", dpi=500)


f_cutoff = ggplot(data, aes(y = f, x = Population, color = Epistasis, size = as.factor(CutOff))) + geom_point(shape = 21, fill = NA, stroke = 1.5) + 
  theme_linedraw() + theme(strip.text.x = element_text(size=10, face="bold", color = "black"),
                           strip.background = element_rect(colour="white", fill="white"), 
                           legend.title=element_text(face="bold"), axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept = 0.01, color = "red") + geom_hline(yintercept = 1, color = "black", size = 2) +
  ylab("Frequency") + coord_trans(y="log2") + scale_size_discrete(name = "Cut-off") + scale_color_manual(values = col)
ggsave("f_cutoff.pdf", plot = f_cutoff, width=17, height=10, units="cm", dpi=500)
ggsave("f_cutoff.png", plot = f_cutoff, width=17, height=10, units="cm", dpi=500)


ggplot(data, aes(x = CutOff, y = NbrStops, color = Population)) + geom_line() + 
  theme_linedraw() + scale_x_continuous(trans='log2') + xlab("Cut-off (log2)") + ylab("Number of LoF SNPs (log2)") +
  theme(legend.position = "None") + scale_y_continuous(trans='log2')
ggsave("fLoF_2.png", plot = last_plot(), width=8, height=5, units="cm", dpi=500)

test = data %>% filter(CutOff == 2504)
test = rep(test$NbrStops, times = 4)
data$fLoF = data$NbrStops/test 
ggplot(data, aes(x = CutOff, y = fLoF, color = Population)) + geom_line() + theme_linedraw() + 
  theme(legend.position = "None") + 
  scale_x_continuous(trans='log2') + xlab("Cut-off (log2)") + ylab("% of LoF SNPs")
ggsave("fLoF.png", plot = last_plot(), width=8, height=5, units="cm", dpi=500)




co2504 = data %>% filter(CutOff == 2504)
co50 = data %>% filter(CutOff == 50)
co2504$FrequentLoF = (co2504$NbrStops - co50$NbrStops)/co2504$NbrStops # Proportion de LoF > 50
ggplot(co2504, aes(x = Epistasis, y = FrequentLoF, color = Epistasis)) + 
  geom_boxplot() + theme_linedraw() +
  stat_compare_means(comparisons = list(c("Antagonistic", "None")), method = "wilcox.test", label.y = 0.09) +
  theme(strip.text.x = element_text(size=10, face="bold", color = "black"),
        strip.background = element_rect(colour="white", fill="white"), 
        legend.title=element_text(face="bold")) +
  scale_color_manual(values=col) + ylab("% of frequent LoF \n(>50 occurences)") +
  geom_jitter() + ylim(0, 0.095)
ggsave("boxplot.png", plot = last_plot(), width=10, height=6, units="cm", dpi=500)
