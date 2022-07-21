library(readxl)
library(viridis)
library(hrbrthemes)
library(ggplot2)

# load clinical data
sheets <- excel_sheets("41588_2020_696_MOESM3_ESM.xlsx")
clinical <- data.frame(read_excel("41588_2020_696_MOESM3_ESM.xlsx",sheet = sheets[1],col_names = TRUE))
colnames(clinical) <- clinical[2,]
clinical <- clinical[-c(1,2),]
# Patient 1 not in RNA analysis
clinical <- clinical[-1,]

#table gender in study group
table(clinical$`Study Group`,clinical$Gender)

# median of age per group
median(as.numeric(clinical$`Age at death`)[clinical$`Study Group`=='Young'])
median(as.numeric(clinical$`Age at death`)[clinical$`Study Group`=='Old'])
median(as.numeric(clinical$`Age at death`)[clinical$`Study Group`=='AD'])

# standard error of age per group
sd(as.numeric(clinical$`Age at death`)[clinical$`Study Group`=='Young'])
sd(as.numeric(clinical$`Age at death`)[clinical$`Study Group`=='Old'])
sd(as.numeric(clinical$`Age at death`)[clinical$`Study Group`=='AD'])

colnames(clinical)[c(2,4)] <- c('condition','age')
clinical$age <- as.numeric(clinical$age)


ggplot(clinical, aes(x=condition, y=age, fill=condition)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.9, option="C") +
  theme_ipsum() +
  scale_fill_brewer(palette="PiYG")+
  labs(x="Groups", y = "Age")

# description: boxplot of age distribution in young, old and AD. No overlap of age 
# distribution in young patients vs old/AD patients

#Cause of Death from young and old patients
table(clinical$`Cause of Death`[1:18])
