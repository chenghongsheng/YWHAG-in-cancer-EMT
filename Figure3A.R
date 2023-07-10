# Figure 3A & S7C

library("Hmisc")
library('gplots')
library('RColorBrewer')
library("genefilter")
library('pheatmap')
library("tidyverse")
library("ggplot2")
library('ggrepel')


##################### DMOG ##################################
DMOG_wide <- read.delim("./DMOG_wide.txt",
                        check.names = F,row.names=2)

DMOG_wide_trim<-DMOG_wide[,c("Inhibitor","DMOG_TotalAUC",'siDMOG_TotalAUC')]
summary(DMOG_wide_trim$DMOG_TotalAUC)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-8.4215 -2.7778 -1.4321 -1.3417 -0.0619  5.3823 
summary(DMOG_wide_trim$siDMOG_TotalAUC)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-4.95465 -1.33220 -0.09692  0.20872  1.54324  6.71424 

DMOG_wide <- DMOG_wide %>%  # classify the YWHAG status
  mutate(status = case_when(DMOG_wide$DMOG_TotalAUC<= -2.7778 & DMOG_wide$siDMOG_TotalAUC > -1.3322 ~ "YHWAG-dependent EMT Inhibitor",
                             DMOG_wide$DMOG_TotalAUC<= -2.7778 & DMOG_wide$siDMOG_TotalAUC <= -1.3322 ~ "YHWAG-independent EMT Inhibitor",
                             DMOG_wide$DMOG_TotalAUC> -2.7778 & DMOG_wide$siDMOG_TotalAUC >= 1.54324 ~ "YHWAG-dependent EMT Activator",
                            TRUE~"No role"))
DMOG_wide$status<-factor(DMOG_wide$status,levels = c("YHWAG-dependent EMT Inhibitor","YHWAG-dependent EMT Activator",
                                                         'YHWAG-independent EMT Inhibitor','No role'))

dotplot_col<-c("red", "Blue","#EF9739",'grey')

DMOG_inhibitor<-ggplot(DMOG_wide,aes(x=DMOG_TotalAUC,y=siDMOG_TotalAUC,color=status))+
  geom_point(size=3)+
  scale_color_manual(values=dotplot_col)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme(panel.background = element_blank(),
        panel.border = element_blank())+
  ggtitle("DMOG")+
  xlab("DMOG")+
  ylab("DMOG + siYHWAG")

tiff('DMOG_inhibitor.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
DMOG_inhibitor
dev.off()

#heatmap
##siDMOG_DMOG

DMOG_siDMOG<-DMOG_wide[,c(1,2:14,41:53)]

colnames(DMOG_siDMOG)<-c( "Inhibitor","0h.DMOG","4h.DMOG","8h.DMOG","12h.DMOG","16h.DMOG","20h.DMOG","24h.DMOG","28h.DMOG",
                           "32h.DMOG","36h.DMOG","40h.DMOG","44h.DMOG","48h.DMOG",
                           "0h.siDMOG","4h.siDMOG","8h.siDMOG","12h.siDMOG","16h.siDMOG","20h.siDMOG","24h.siDMOG","28h.siDMOG","32h.siDMOG","36h.siDMOG",
                           "40h.siDMOG","44h.siDMOG","48h.siDMOG" )
DMOG_siDMOG$status<-DMOG_wide$status
DMOG_siDMOG_shortlist<-DMOG_siDMOG[!DMOG_siDMOG$status=="No role",]
DMOG_siDMOG_row<-DMOG_siDMOG_shortlist[,28,drop=F]
DMOG_siDMOG_mat<-as.matrix(DMOG_siDMOG_shortlist[,-c(1,28)])

break_LL = seq(min(DMOG_siDMOG_mat), 1,length.out=50)
break_HL = seq(1, max(DMOG_siDMOG_mat),length.out=50)[-1]
breaks<-c(break_LL,break_HL)
col = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(length(breaks))

my_color_annotation<-list(status = c('YHWAG-dependent EMT Inhibitor' = "red", 'YHWAG-dependent EMT Activator' = "Blue",'YHWAG-independent EMT Inhibitor'="#EF9739"))

heatmap_DMOG_siDMOG<-pheatmap(DMOG_siDMOG_mat,scale="none",border_color = NA,cluster_cols = F,cluster_rows = T,
                              show_rownames = F,show_colnames = T,breaks=breaks,color = col,
                              annotation_row=DMOG_siDMOG_row,annotation_colors = my_color_annotation,
                              clustering_distance_rows = "euclidean",fontsize_row = 3,gaps_col = 13)

tiff('heatmap_DMOG_siDMOG.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
heatmap_DMOG_siDMOG
dev.off()

##################### TGFb ##################################
TGFB_wide <- read.delim("./TGFB_wide.txt",
                        check.names = F,row.names=2)

TGFB_wide_trim<-TGFB_wide[,c("Inhibitor","TGFB_TotalAUC",'siTGFB_TotalAUC')]

summary(c(TGFB_wide$TGFB_TotalAUC))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-6.2098 -2.7314 -1.3941 -1.5090 -0.5253  5.7008

summary(c(TGFB_wide$siTGFB_TotalAUC))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-5.8433 -2.2076 -0.9039 -0.8362  0.4320  4.9810 

TGFB_wide <- TGFB_wide %>% 
  mutate(status = case_when(TGFB_wide$TGFB_TotalAUC<= -2.7318 & TGFB_wide$siTGFB_TotalAUC > -2.2076 ~ "YHWAG-dependent EMT Inhibitor",
                             TGFB_wide$TGFB_TotalAUC<= -2.7318 & TGFB_wide$siTGFB_TotalAUC <= -2.2076 ~ "YHWAG-independent EMT Inhibitor",
                             TGFB_wide$TGFB_TotalAUC> -2.7318 & TGFB_wide$siTGFB_TotalAUC >= 0.4320 ~ "YHWAG-dependent EMT Activator"))
TGFB_wide$status<-replace_na(TGFB_wide$status,"No role")
TGFB_wide$status<-factor(TGFB_wide$status,levels = c("YHWAG-dependent EMT Inhibitor","YHWAG-dependent EMT Activator",
                                                         'YHWAG-independent EMT Inhibitor','No role'))

TGFB_inhibitor<-ggplot(TGFB_wide,aes(x=TGFB_TotalAUC,y=siTGFB_TotalAUC,color=status))+
  geom_point(size=3)+
  scale_color_manual(values=dotplot_col)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  theme(panel.background = element_blank(),
        panel.border = element_blank())+
  ggtitle("TGFB")+
  xlab("TGFB")+
  ylab("TGFB + siYHWAG")

tiff('TGFB_inhibitor.tiff',width=2500,height=2000,units='px',res=300,compression='lzw')
TGFB_inhibitor
dev.off()

#heatmap
##siTGFB_TGFB

TGFB_siTGFB<-TGFB_wide[,c(1,2:14,41:53)]

colnames(TGFB_siTGFB)<-c( "Inhibitor","0h.TGFB","4h.TGFB","8h.TGFB","12h.TGFB","16h.TGFB","20h.TGFB","24h.TGFB","28h.TGFB",
                          "32h.TGFB","36h.TGFB","40h.TGFB","44h.TGFB","48h.TGFB",
                          "0h.siTGFB","4h.siTGFB","8h.siTGFB","12h.siTGFB","16h.siTGFB","20h.siTGFB","24h.siTGFB","28h.siTGFB","32h.siTGFB","36h.siTGFB",
                          "40h.siTGFB","44h.siTGFB","48h.siTGFB" )
TGFB_siTGFB$status<-TGFB_wide$status
TGFB_siTGFB_shortlist<-TGFB_siTGFB[!TGFB_siTGFB$status=="No role",]
TGFB_siTGFB_row<-TGFB_siTGFB_shortlist[,28,drop=F]
TGFB_siTGFB_mat<-as.matrix(TGFB_siTGFB_shortlist[,-c(1,28)])

heatmap_TGFB_siTGFB<-pheatmap(TGFB_siTGFB_mat,scale="none",border_color = NA,cluster_cols = F,cluster_rows = T,
                              show_rownames = F,show_colnames = T,breaks=breaks,color = col,
                              annotation_row=TGFB_siTGFB_row,annotation_colors = my_color_annotation,
                              clustering_distance_rows = "euclidean",fontsize_row = 3,gaps_col = 13)

tiff('heatmap_TGFB_siTGFB.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
heatmap_TGFB_siTGFB
dev.off()


############# End ################