library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(TCGAbiolinksGUI.data)
library("TCGAbiolinks")
library(dplyr)

########### Download data from TCGA ##############
##TCGA-GBM
query_GBM = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal"))

GDCdownload(query = query_GBM)
GBM_data = GDCprepare(query_GBM)


##TCGA-LGG

query_LGG = GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal"))

GDCdownload(query = query_LGG)
LGG_data = GDCprepare(query_LGG)

##TCGA-ACC

query_ACC = GDCquery(
  project = "TCGA-ACC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor"))

GDCdownload(query = query_ACC)
ACC_data <- GDCprepare(query_ACC)


## TCGA-BLCA

query_BLCA = GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query = query_BLCA)
BLCA_data <- GDCprepare(query_BLCA)


## TCGA-BRCA

query_BRCA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_BRCA)
BRCA_data <- GDCprepare(query_BRCA)


## TCGA-CESC

query_CESC = GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_CESC)
CESC_data <- GDCprepare(query_CESC)

## TCGA-CHOL

query_CHOL = GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal"))

GDCdownload(query = query_CHOL)
CHOL_data <- GDCprepare(query_CHOL)

## TCGA-COAD

query_COAD = GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_COAD)
COAD_data <- GDCprepare(query_COAD)

## TCGA-DLBC

query_DLBC = GDCquery(
  project = "TCGA-DLBC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor"))

GDCdownload(query = query_DLBC)
DLBC_data <- GDCprepare(query_DLBC)

## TCGA-ESCA

query_ESCA = GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_ESCA)
ESCA_data <- GDCprepare(query_ESCA)

## TCGA-HNSC
query_HNSC = GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_HNSC)
HNSC_data <- GDCprepare(query_HNSC)

## TCGA-KICH

query_KICH = GDCquery(
  project = "TCGA-KICH",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query = query_KICH)
KICH_data <- GDCprepare(query_KICH)

## TCGA-KIRC

query_KIRC = GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query = query_KIRC)
KIRC_data <- GDCprepare(query_KIRC)

## TCGA-KIRP

query_KIRP = GDCquery(
  project = "TCGA-KIRP",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query = query_KIRP)
KIRP_data <- GDCprepare(query_KIRP)

## TCGA-LAML

query_LAML = GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Blood Derived Cancer - Peripheral Blood"))

GDCdownload(query = query_LAML)
LAML_data <- GDCprepare(query_LAML)

## TCGA-LIHC

query_LIHC = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal"))

GDCdownload(query = query_LIHC)
LIHC_data <- GDCprepare(query_LIHC)


## TCGA-LUAD

query_LUAD = GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal"))

GDCdownload(query = query_LUAD)
LUAD_data <- GDCprepare(query_LUAD)


## TCGA-LUSC

query_LUSC = GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query = query_LUSC)
LUSC_data <- GDCprepare(query_LUSC)


## TCGA-MESO

query_MESO = GDCquery(
  project = "TCGA-MESO",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor"))

GDCdownload(query = query_MESO)
MESO_data <- GDCprepare(query_MESO)


## TCGA-OV

query_OV = GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor"))

GDCdownload(query = query_OV)
OV_data <- GDCprepare(query_OV)


## TCGA-PAAD

query_PAAD = GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_PAAD)
PAAD_data <- GDCprepare(query_PAAD)


## TCGA-PCPG

query_PCPG = GDCquery(
  project = "TCGA-PCPG",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal","Additonal - New Primary"))

GDCdownload(query = query_PCPG)
PCPG_data <- GDCprepare(query_PCPG)


## TCGA-PRAD

query_PRAD = GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_PRAD)
PRAD_data <- GDCprepare(query_PRAD)


## TCGA-READ

query_READ = GDCquery(
  project = "TCGA-READ",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal"))

GDCdownload(query = query_READ)
READ_data <- GDCprepare(query_READ)


## TCGA-SARC

query_SARC = GDCquery(
  project = "TCGA-SARC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_SARC)
SARC_data <- GDCprepare(query_SARC)



## TCGA-SKCM

query_SKCM = GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_SKCM)
SKCM_data <- GDCprepare(query_SKCM)


## TCGA-STAD


query_STAD = GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query = query_STAD)
STAD_data <- GDCprepare(query_STAD)



## TCGA-TGCT


query_TGCT = GDCquery(
  project = "TCGA-TGCT",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Additional - New Primary"))

GDCdownload(query = query_TGCT)
TGCT_data <- GDCprepare(query_TGCT)


## TCGA-THCA

query_THCA = GDCquery(
  project = "TCGA-THCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Metastatic","Solid Tissue Normal"))

GDCdownload(query = query_THCA)
THCA_data <- GDCprepare(query_THCA)


## TCGA-THYM


query_THYM = GDCquery(
  project = "TCGA-THYM",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal"))

GDCdownload(query = query_THYM)
THYM_data <- GDCprepare(query_THYM)


## TCGA-UCEC

query_UCEC = GDCquery(
  project = "TCGA-UCEC",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal"))

GDCdownload(query = query_UCEC)
UCEC_data <- GDCprepare(query_UCEC)


## TCGA-UCS


query_UCS = GDCquery(
  project = "TCGA-UCS",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor"))

GDCdownload(query = query_UCS)
UCS_data <- GDCprepare(query_UCS)



## TCGA-UVM

query_UVM = GDCquery(
  project = "TCGA-UVM",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type="Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor"))

GDCdownload(query = query_UVM)
UVM_data <- GDCprepare(query_UVM)

########### Extract TPM data ##############

selected_data<-c('ACC_data','BRCA_data','CESC_data','CHOL_data',
                 'COAD_data','DLBC_data','ESCA_data','HNSC_data',
                 'KIRC_data','KIRP_data','LIHC_data','LUAD_data',
                 'LUSC_data','PAAD_data','PRAD_data','READ_data',
                 'SKCM_data','STAD_data','TGCT_data','THYM_data',
                 'UCS_data')

tpm_ACC<-assays(ACC_data)$tpm_unstrand
write.csv(tpm_ACC,file="tpm_ACC.csv")

tpm_BRCA<-assays(BRCA_data)$tpm_unstrand
write.csv(tpm_BRCA,file="tpm_BRCA.csv")

tpm_CESC<-assays(CESC_data)$tpm_unstrand
write.csv(tpm_CESC,file="tpm_CESC.csv")

tpm_CHOL<-assays(CHOL_data)$tpm_unstrand
write.csv(tpm_CHOL,file="tpm_CHOL.csv")

tpm_COAD<-assays(COAD_data)$tpm_unstrand
write.csv(tpm_COAD,file="tpm_COAD.csv")

tpm_DLBC<-assays(DLBC_data)$tpm_unstrand
write.csv(tpm_DLBC,file="tpm_DLBC.csv")

tpm_ESCA<-assays(ESCA_data)$tpm_unstrand
write.csv(tpm_ESCA,file="tpm_ESCA.csv")

tpm_HNSC<-assays(HNSC_data)$tpm_unstrand
write.csv(tpm_HNSC,file="tpm_HNSC.csv")

tpm_KIRC<-assays(KIRC_data)$tpm_unstrand
write.csv(tpm_KIRC,file="tpm_KIRC.csv")

tpm_KIRP<-assays(KIRP_data)$tpm_unstrand
write.csv(tpm_KIRP,file="tpm_KIRP.csv")

tpm_LIHC<-assays(LIHC_data)$tpm_unstrand
write.csv(tpm_LIHC,file="tpm_LIHC.csv")

tpm_LUAD<-assays(LUAD_data)$tpm_unstrand
write.csv(tpm_LUAD,file="tpm_LUAD.csv")

tpm_LUSC<-assays(LUSC_data)$tpm_unstrand
write.csv(tpm_LUSC,file="tpm_LUSC.csv")

tpm_PAAD<-assays(PAAD_data)$tpm_unstrand
write.csv(tpm_PAAD,file="tpm_PAAD.csv")

tpm_PRAD<-assays(PRAD_data)$tpm_unstrand
write.csv(tpm_PRAD,file="tpm_PRAD.csv")

tpm_READ<-assays(READ_data)$tpm_unstrand
write.csv(tpm_READ,file="tpm_READ.csv")

tpm_SKCM<-assays(SKCM_data)$tpm_unstrand
write.csv(tpm_SKCM,file="tpm_SKCM.csv")

tpm_STAD<-assays(STAD_data)$tpm_unstrand
write.csv(tpm_STAD,file="tpm_STAD.csv")

tpm_TGCT<-assays(TGCT_data)$tpm_unstrand
write.csv(tpm_TGCT,file="tpm_TGCT.csv")

tpm_THYM<-assays(THYM_data)$tpm_unstrand
write.csv(tpm_THYM,file="tpm_THYM.csv")

tpm_UCS<-assays(UCS_data)$tpm_unstrand
write.csv(tpm_UCS,file="tpm_UCS.csv")


########### Plot correlation ##############

tmp_csv<-list.files(pattern="tpm_")

tpm_list <- lapply(tmp_csv, function(x) read.csv(x, check.names = F,row.names = 1))
selected_genes<-c("ENSG00000170027.7",'ENSG00000080824.19','ENSG00000116898.12') #YWHAG, HSP90AA1, MRPS15


tpm_list_sub <- lapply(tpm_list, function(x) {x[selected_genes,]})
tpm_list_sub <- lapply(tpm_list_sub, function(x) as.data.frame(t(x)))


#rename list in tpm_list_sub
dataname<-str_split_i(selected_data,"_",1)

names(tpm_list_sub)<-dataname

#rbind 
tpm_df<-tpm_list_sub %>% 
  bind_rows(.id = "tumor")

colnames(tpm_df)<-c("tumor",'YWHAG','HSP90AA1','MRPS15')
tpm_df$log2YWHAG<-log2(tpm_df$YWHAG)
tpm_df$log2HSP90AA1<-log2(tpm_df$HSP90AA1)
tpm_df$log2MRPS15<-log2(tpm_df$MRPS15)


#plot graph

HSPvYWHAG_corplot<-ggplot(tpm_df,aes(x=log2YWHAG,y=log2HSP90AA1))+
  geom_point()+
  facet_wrap(~tumor,ncol=3)+
  geom_smooth(method="lm",color="red",se=F)+
  stat_cor(aes(label = ..r.label..), color = "red", geom = "label",
           p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab(bquote(log[2](YWHAG~TPM)))+
  ylab(bquote(log[2](HSP90AA1~TPM)))
#coord_fixed()

dev.new()
pdf("HSPvYWHAG_corplot.pdf",width = 10,height = 10)
HSPvYWHAG_corplot
dev.off()




HSPvMRPS15_corplot<-ggplot(tpm_df,aes(x=log2MRPS15,y=log2HSP90AA1))+
  geom_point()+
  facet_wrap(~tumor,ncol=3)+
  geom_smooth(method="lm",color="red",se=F)+
  stat_cor(aes(label = ..r.label..), color = "red", geom = "label",
           p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab(bquote(log[2](MRPS15~TPM)))+
  ylab(bquote(log[2](HSP90AA1~TPM)))


dev.new()
pdf("HSPvMRPS15_corplot.pdf",width = 10,height = 10)
HSPvMRPS15_corplot
dev.off()



MRPS15vYWHAG_corplot<-ggplot(tpm_df,aes(x=log2YWHAG,y=log2MRPS15))+
  geom_point()+
  facet_wrap(~tumor,ncol=3)+
  geom_smooth(method="lm",color="red",se=F)+
  stat_cor(aes(label = ..r.label..), color = "red", geom = "label",
           p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab(bquote(log[2](YWHAG~TPM)))+
  ylab(bquote(log[2](MRPS15~TPM)))
#coord_fixed()

dev.new()
pdf("MRPS15vYWHAG_corplot.pdf",width = 10,height = 10)
MRPS15vYWHAG_corplot
dev.off()
