
library(limma)
library(Biobase)
library(readxl)
library(plyr)
library(impute)
library(imputeLCMD)

path="C:\\Users\\laslo\\Documents\\Protein\\"


#phenotipic data

pData<-read_excel(paste0(path,"2020-10-22-TP-Inventory.xlsx"),
                  sheet = 1, range = "B16:N48")
pData<-as.data.frame(pData)
colnames(pData)<-make.names(colnames(pData))
names(pData)

unique(as.data.frame(pData)$Age)
newAGES=mapvalues(as.data.frame(pData)$Age, from= c("Newborn", "5 YR",    
                                                    "1 YR",    "12 YR",   "3 MO" ,   "2 MO" ,  
                                                    "9 YR" ,   "3 DA"),
                  to=c(0,5,1,12,0,0,9,0))
as.numeric(as.character(newAGES))
ind=which(colnames(pData)=="Age")
ind
colnames(pData)
pData[,ind]<-as.numeric(as.character(newAGES))
head(pData)

pData_t<-pData %>% 
  mutate(Group=
  case_when(Treatment.Description == "Untreated" 
            & grepl("Healthy", Sample.Description, fixed = TRUE)
            ~ "Ctrl UT",
            Treatment.Description != "Untreated" 
            & grepl("Healthy", Sample.Description, fixed = TRUE)
            ~ "Ctrl AOAA",
            Treatment.Description == "Untreated" 
            & grepl("Trisomy", Sample.Description, fixed = TRUE)
            ~ "DS UT",
            Treatment.Description != "Untreated" 
            & grepl("Trisomy", Sample.Description, fixed = TRUE)
            ~ "DS AOAA"))
pData_t


xData<-read_excel(paste0(path,"Batch_corrected_log_intensities_data.xlsx"),
           sheet = 1, range = "A1:AG5539")

xData_m<-as.matrix(xData[2:ncol(xData)])
rownames(xData_m)<-xData[[1]]
colnames(xData_m)<-make.names(colnames(xData_m))

rownames(pData_t)<-colnames(xData_m)

head(xData_m)

eset <- ExpressionSet(assayData = xData_m,
                      phenoData = AnnotatedDataFrame(pData_t))
                                                     
                      
head(exprs(eset))
boxplot(exprs(eset))

saveRDS(eset, file = paste0(path,"Protein_eset.rds"))

###############################################################################

##############################################################################
library(dplyr)
library(readxl)
library(Biobase)
library(KEGG.db)
library(org.Hs.eg.db)
library(gplots)
library(HGNChelper)
library(Biobase)
library(writexl)
library(org.Hs.eg.db)


path="C:\\Users\\laslo\\Documents\\Protein\\"

eset <- readRDS(paste0(path,"Protein_eset.rds"))

names(pData(eset))
pData(eset)$Group <-factor(pData(eset)$Group)
pData(eset)$Gender <-factor(pData(eset)$Gender)

summary(pData(eset))

levels(pData(eset)$Group)
pData(eset)$Group <- factor(pData(eset)$Group, levels = c( "Ctrl UT"  , "Ctrl AOAA","DS UT" ,    "DS AOAA" ))  

rownames(exprs(eset))

names(pData(eset))



#Contrast5...............................
#Positive interaction- AOAA modifies slightly the metabolite level on Ctrl cell lines,
#but it has a stronger effect on DS
#Negative interaction- AOAA modifies the metabolite level on Ctrl and DS cells,
#but with an opposite tendency




df_interaction<-read_excel(paste0(path,"005_interaction_genelist.xlsx"),
                                sheet = 1)

df_interaction<- as.data.frame(df_interaction)

names(df_interaction)

df_interaction_sel<-df_interaction %>% 
  filter(P.Value<0.05 & FC<0) %>% 
  dplyr::select(-Reactome.Pathways) %>% 
  arrange(P.Value)

writexl::write_xlsx(df_interaction_sel, paste0(path,"PosInteraction.xlsx"))



#selected_metabolite="P05386 "

################################################################
make_boxplot<-function(x, 
                       selected_metabolite="",
                       subtitle="",
                       titleadd="",
                       group, pair) {
  
require(ggplot2)
  stopifnot(is(x, 'ExpressionSet'))
  if (!missing(selected_metabolite)) {
    stopifnot(selected_metabolite %in% rownames(exprs(x)))
  }
  if (!missing(group)) {
    stopifnot(group %in% colnames(pData(x)))
  }
  if (!missing(pair)) {
    stopifnot(pair %in% colnames(pData(x)))
  }
  
  Values <-exprs(x)[rownames(exprs(x))==selected_metabolite,]
  Group<-pData(x)[[group]]
  
  names(pData(x))
  
  df <-data.frame(Values=exprs(x)[rownames(exprs(x))==selected_metabolite,],
                  Group=pData(x)[[group]])
  
  r <- runif(nrow(df), -0.1, 0.1)
  
  names(df)
  ggp<- ggplot(df)+
    aes(x=Group, y=Values, fill=Group) +
    geom_boxplot(outlier.colour=NA,outlier.shape = NA, coef=0, 
                  colour="grey20", alpha=0.5) + 
    {if (!missing(pair)) {
                     geom_line(aes(x  = as.numeric(Group)+r, y = Values,
                  group = pData(x)[[pair]]), linetype = "dashed")}} +
    geom_point(aes(x  = as.numeric(Group)+r, fill=Group),
               size=4, shape=21, colour="grey20") +
    labs(title=paste(selected_metabolite,titleadd),
         subtitle= subtitle,
         y="Log2 normalised intensity")+
    theme_bw() +
    theme(axis.line=element_line(colour="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_blank())+
    theme(axis.text.x = element_text( color="black", 
                                      size=10),
          axis.text.y = element_text( color="black", 
                                      size=10))+
    theme(legend.position = "none")+
    theme(axis.title.x=element_blank())
  
  ggp
  
}

library(HGNChelper)
new.hgnc.table <- getCurrentHumanMap()

head(new.hgnc.table)

#Create df_ana using df_interaction_sel$Id Uniprot names
df_ana<-AnnotationDbi::select(org.Hs.eg.db,keys=df_interaction_sel$Id, keytype = "UNIPROT", c("SYMBOL", "GENENAME"))
df_ana
df_ana<-na.omit(df_ana)

names(df_ana)



#Genes
ch21genes<-readRDS(paste0(path,"ch21genes.Rds"))
ch21genes<-checkGeneSymbols(ch21genes, map=new.hgnc.table)[[3]]
ch21genes<-ch21genes[!is.na(ch21genes)]
ch21genes

df_ana$ischr21 <-ifelse(df_ana$SYMBOL %in% ch21genes, "chr21", "")
df_ana

make_boxplot(eset,df_ana$UNIPROT[6],subtitle=df_ana$GENENAME[6], titleadd = df_ana$ischr21[6],
             group = "Group", pair = "Cell.Line")
unlink(paste0(paste0(path,"Figures\\"),
       list.files(paste0(path,"Figures\\"))))

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Figures\\"), 
         height = 0.75/20, width=1/20, scale=4*20)  
}


##################Difference in DS vs. Ctrl
list.files(paste0(path))

df_typeDS<-read_excel(paste0(path,"004_type_DS_genelist.xlsx"),
                           sheet = 1)
names(df_typeDS)[2:ncol(df_typeDS)]<-
  paste0(names(df_typeDS)[2:ncol(df_typeDS)], "_DS")

df_typeUT<-read_excel(paste0(path,"001_type_UT_genelist.xlsx"),
                      sheet = 1)

names(df_typeUT)[2:ncol(df_typeUT)]<-
  paste0(names(df_typeUT)[2:ncol(df_typeUT)], "_UT")
  

df_tot<-merge(df_typeUT, df_typeDS, by="Id")
names(df_tot)
df_tot$Prod.P.Value <-df_tot$P.Value_DS*df_tot$P.Value_UT


df_tot<- as.data.frame(df_tot)

names(df_tot)

df_tot_sel<-df_tot %>% 
  filter(Prod.P.Value<0.05) %>%  
    filter((FC_UT< (-1.05) & FC_DS> 1.05)|(FC_UT>1.05 & FC_DS<(-1.05)))%>% 
  dplyr::select(-Reactome.Pathways_UT,-Reactome.Pathways_DS) %>% 
  arrange(Prod.P.Value)


dim(df_tot_sel)
head(df_tot_sel)


writexl::write_xlsx(df_tot_sel, paste0(path,"CTR_DS_DSAOAA.xlsx"))


#Create df_ana using df_tot_sel$Id Uniprot names
df_ana<-AnnotationDbi::select(org.Hs.eg.db,
                              keys=df_tot_sel$Id, 
                              keytype = "UNIPROT",
                              c("SYMBOL", "GENENAME"))
df_ana
df_ana<-na.omit(df_ana)

names(df_ana)



#Genes
ch21genes<-readRDS(paste0(path,"ch21genes.Rds"))
ch21genes<-checkGeneSymbols(ch21genes, map=new.hgnc.table)[[3]]
ch21genes<-ch21genes[!is.na(ch21genes)]
ch21genes

df_ana$ischr21 <-ifelse(df_ana$SYMBOL %in% ch21genes, "chr21", "")

dim(df_ana)
df_ana
make_boxplot(eset,df_ana$UNIPROT[1],subtitle=df_ana$GENENAME[6], titleadd = df_ana$ischr21[6],
             group = "Group", pair = "Cell.Line")
unlink(paste0(paste0(path,"CTRLvsDSvsDS_AOAA\\"),
list.files(paste0(path,"CTRLvsDSvsDS_AOAA\\"))))

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"CTRLvsDSvsDS_AOAA\\"), 
         height = 0.75/20, width=1/20, scale=4*20)  
}


h#Energymetabolites

#User defined function#############################################
Collect_GenesfromPathway<-function(pName){
  res3 <-NULL
  for (i in 1:length(pName)){
    pId <- mget(pName[i], KEGGPATHNAME2ID)[[1]]
    pId
    Geneid<-mget(pId, org.Hs.egPATH2EG)
    res<-mapIds(org.Hs.eg.db, Geneid[[1]], 'SYMBOL', 'ENTREZID')
    res2<-unname(res)
    res3 <- c(res2, res3)
  }
  return(unique(res3))}
##############################################



pName<-"Glycolysis / Gluconeogenesis"
genes_glyc<-Collect_GenesfromPathway(pName)
genes_glyc <-genes_glyc[!genes_glyc %in% c("ADH1A",   "ADH1B" ,  "ADH1C",   "ADH4" ,
                                           "ADH5" ,   "ADH6" ,   "ADH7", "ALDH2" ,
                                           "ALDH3A1", "ALDH1B1", "ALDH1A3","ALDH3B1",
                                           "ALDH3B2", "ALDH9A1","ALDH3A2" ,"ALDH7A1")]
#"PDK1", "PGM2"
#genes_glyc <-c(genes_glyc, "SLC2A1" ,  "SLC2A2" ,  "SLC2A3" ,  "SLC2A4",
#               "SLC16A3",  "SLC16A4")
genes_glyc
genes_glyc<-checkGeneSymbols(genes_glyc, map=new.hgnc.table)[[3]]
genes_glyc<-genes_glyc[order(genes_glyc)]

#Transporters
genes_transp <-c(
  #Glucose transporters
  "SLC2A1" ,  "SLC2A2" ,  "SLC2A3" ,  "SLC2A4", "SLC2A14",
  #Lactate transporters
  "SLC16A1", "SLC16A7", "SLC16A8", "SLC16A3", "SLC16A4",
  #Sodium-glucose cotransporters
  "SLC5A1", "SLC5A2")
genes_transp<-checkGeneSymbols(genes_transp, map=new.hgnc.table)[[3]]
#genes_transp<-genes_transp[order(genes_transp)]
genes_transp

pName<-"Citrate cycle (TCA cycle)"
genes_TCA<-Collect_GenesfromPathway(pName)
genes_TCA

pName<-"Pentose phosphate pathway"
genes_PPP<-Collect_GenesfromPathway(pName)
genes_PPP

pName <- c("Fatty acid degradation",
           "Fatty acid biosynthesis",
           "Fatty acid elongation")
genes_FA<-Collect_GenesfromPathway(pName)
genes_FA<-c(genes_FA,"FABP1","FABP2","FABP3","FABP4","FABP5","FABP6",
            "FABP7","PMP2","FABP9","FABP12")
#"GPAT","AGPAT2","PNPLA3","DGAT2","PLIN2","LPIN1", "NTTP"
genes_FA<-genes_FA[order(genes_FA)]
genes_FA<-unique(genes_FA)
genes_FA

#BiocManager::install("reactome.db")
library(reactome.db)
xx <- as.list(reactomePATHID2EXTID)
names(xx)
geneids<-xx[["R-HSA-8964539"]]
res<-mapIds(org.Hs.eg.db, geneids, 'SYMBOL', 'ENTREZID')
genes_Glutamine <-unname(res)
genes_Glutamine

#####################
FindUniprotGenename <- function(selectedgenes){
  library(org.Hs.eg.db)
  keytypes(org.Hs.eg.db)
  df<-AnnotationDbi::select(org.Hs.eg.db,keys=selectedgenes, 
                            keytype = "SYMBOL", 
                            c("UNIPROT", "GENENAME"))
  df<-na.omit(df)
  df<-df[df$UNIPROT %in% rownames(exprs(eset)),]
  df$ischr21 <-ifelse(df$SYMBOL %in% ch21genes, "chr21", "")
  return(df)
}
##################################$


unlink(paste0(paste0(path,"Glyc\\"),
              list.files(paste0(path,"Glyc\\"))))

unlink(paste0(paste0(path,"\\TCA"),
              list.files(paste0(path,"\\TCA"))))

unlink(paste0(paste0(path,"\\FA"),
              list.files(paste0(path,"\\FA"))))

unlink(paste0(paste0(path,"\\PPP"),
              list.files(paste0(path,"\\PPP"))))

unlink(paste0(paste0(path,"\\Glut"),
              list.files(paste0(path,"\\Glut"))))


#Glycolysis
df_ana<-FindUniprotGenename(genes_glyc)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Glyc\\"), 
         height = 0.75/20, width=1/20, scale=4*20)  
}

#TCA
df_ana<-FindUniprotGenename(genes_TCA)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"TCA\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
  }

#Mito3
df_ana<-FindUniprotGenename(genes_FA)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"FA\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
  }

#Glutamine
df_ana<-FindUniprotGenename(genes_Glutamine)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Glut\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
  }

#PPP
df_ana<-FindUniprotGenename(genes_PPP)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"PPP\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
  }






# ####################################
library(readxl)
library(tokenizers)

complex1genes<-("NDUFA1 NDUFA2 NDUFA3 NDUFA4 NDUFA4L NDUFA4L2 NDUFA5
                NDUFA6 NDUFA7 NDUFA8 NDUFA9 NDUFA10 NDUFA11 NDUFA12
                NDUFA13 NDUFAB1 NDUFAF1 NDUFAF2 NDUFAF3 NDUFAF4 NDUFB1
                NDUFB2 NDUFB3 NDUFB4 NDUFB5 NDUFB6 NDUFB7 NDUFB8 NDUFB9
                NDUFB10 NDUFB11 NDUFC1 NDUFC2 NDUFS1 NDUFS2 NDUFS3 NDUFS4
                NDUFS5 NDUFS6 NDUFS7 NDUFS8 NDUFV1 NDUFV2 NDUFV3 MT-ND1 
                MT-ND2 MT-ND3 MT-ND4 MT-ND4L MT-ND5 MT-ND6")

complex2genes<-("SDHA
                SDHB
                SDHC
                SDHD")

complex3genes <-("MT-CYB CYC1 CYCS UQCRFS1 UQCRB UQCRH UQCRC2 UQCRC1 
                 UQCR UQCR10 TTC19 UQCRQ UQCRHL UQCR11")

complex4genes <-("MT-CO1 MT-CO2  MT-CO3 COX4I1 COX4I2
COX5A COX5B COX6A1 COX6A2 COX6B1 COX6B2 
COX6C COX7A2 COX7A1 COX7B COX7C COX8A 
COA1 COA3 COA4 COA5 COA6 COA7 COX11
COX14 COX15 COX16 COX17 COX18 COX19 COX20 COX10
COX7A2L COX7B2 COX8C
")

complex5genes <- ("ATP5A1  ATP5B  ATP5C1  ATP5C2  ATP5D  ATP5E 
    ATP5F1  ATP5G1  ATP5G2  ATP5G3  ATP5H   ATP5I   ATP5J 
    ATP5J2  ATP5L  ATP5L2  ATP5O  ATP5S
                  ATP12A ATP4A ATP4B              
    ATP6V1A ATP6V1B1 ATP6V1B2           
    ATP6V0C ATP6V1C1 ATP6V1E1           
    ATP6V0B ATP6V1G2 ATP6V0A1      
    ATP6AP1 
    ATP6V0E1 ATP6V0D1  ATP6V1F ATP6V1G1             
    ATP6V0A2  ATP6V0A4 ATP6V1D ATP6V1H            
    ATP6V1E2 ATP6V1G3           
    ATP6V0E2 ATP6V0D2 ATP6V1C2 
    ATP5F1A ATP5F1B ATP5F1C ATP5F1D ATP5F1E ATP5PB             
    ATP5MC1 ATP5MC2 ATP5MC3          
    ATP5ME ATP5PF ATP5PO ATP5MF ATP5PD ATP5MG
    MT-ATP6 MT-ATP8 "               )


complex1gene_list <-unlist(tokenize_ptb(complex1genes))
complex1<-checkGeneSymbols(complex1gene_list,map=new.hgnc.table)[[3]]
complex1

complex2gene_list <-unlist(tokenize_ptb(complex2genes))
complex2<-checkGeneSymbols(complex2gene_list,map=new.hgnc.table)[[3]]
complex2

complex3gene_list <-unlist(tokenize_ptb(complex3genes))
complex3<-checkGeneSymbols(complex3gene_list,map=new.hgnc.table)[[3]]
complex3

complex4gene_list <-unlist(tokenize_ptb(complex4genes))
complex4<-checkGeneSymbols(complex4gene_list,map=new.hgnc.table)[[3]]
complex4

complex5gene_list <-unlist(tokenize_ptb(complex5genes))
complex5<-checkGeneSymbols(complex5gene_list,map=new.hgnc.table)[[3]]
complex5

unlink(paste0(paste0(path,"Mito1\\"),
              list.files(paste0(path,"Mito1\\"))))
unlink(paste0(paste0(path,"Mito2\\"),
              list.files(paste0(path,"Mito2\\"))))
unlink(paste0(paste0(path,"Mito3\\"),
              list.files(paste0(path,"Mito3\\"))))
unlink(paste0(paste0(path,"Mito4\\"),
              list.files(paste0(path,"Mito4\\"))))
unlink(paste0(paste0(path,"Mito5\\"),
              list.files(paste0(path,"Mito5\\"))))
#Mito1
df_ana<-FindUniprotGenename(complex1)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Mito1\\"), 
         height = 0.75/20, width=1/20, scale=4*20)  
}


#Mito2
df_ana<-FindUniprotGenename(complex2)
AnnotationDbi::select(org.Hs.eg.db,keys=complex2, 
                      keytype = "SYMBOL", 
                      c("UNIPROT", "GENENAME"))
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Mito2\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
  }

#Mito3
df_ana<-FindUniprotGenename(complex3)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Mito3\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
  }

#Mito4
df_ana<-FindUniprotGenename(complex4)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Mito4\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
  }

#Mito5
df_ana<-FindUniprotGenename(complex5)
df_ana

for(i in 1:nrow(df_ana)){
  make_boxplot(eset,df_ana$UNIPROT[i],
               subtitle=df_ana$GENENAME[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$UNIPROT[i],".pdf"), 
         path=paste0(path,"Mito5\\"), 
         height = 0.75/20, width=1/20, scale=4*20)
}

exprs(eset)[rownames(exprs(eset))=="P36639"]

