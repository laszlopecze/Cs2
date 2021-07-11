
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
            ~ "Ctrl_UT",
            Treatment.Description != "Untreated" 
            & grepl("Healthy", Sample.Description, fixed = TRUE)
            ~ "Ctrl_T",
            Treatment.Description == "Untreated" 
            & grepl("Trisomy", Sample.Description, fixed = TRUE)
            ~ "DS_UT",
            Treatment.Description != "Untreated" 
            & grepl("Trisomy", Sample.Description, fixed = TRUE)
            ~ "DS_T"))
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
path="C:\\Users\\laslo\\Documents\\Protein\\"

eset <- readRDS(paste0(path,"Protein_eset.rds"))



names(pData(eset))
pData(eset)$Group <-factor(pData(eset)$Group)
pData(eset)$Gender <-factor(pData(eset)$Gender)

summary(pData(eset))



#Contrast5...............................
#Positive interaction- AOAA modifies slightly the metabolite level on Ctrl cell lines,
#but it has a stronger effect on DS
#Negative interaction- AOAA modifies the metabolite level on Ctrl and DS cells,
#but with an opposite tendency



df_interaction_sign<-read_excel(paste0(path,"005_interaction_genelist_sign.xlsx"),
                  sheet = 1)

df_interaction_sign<- as.data.frame(df_interaction_sign)
df_interaction_sign

df_interaction<-read_excel(paste0(path,"005_interaction_genelist.xlsx"),
                                sheet = 1)

df_interaction<- as.data.frame(df_interaction)
df_interaction
names(df_interaction)
df_interaction_sel<-df_interaction %>% 
  filter(P.Value<0.05 & FC<0) %>% 
  select(-Reactome.Pathways) %>% 
  arrange(P.Value)

df_interaction_sel
install.packages("writexl")
library(writexl)
writexl::write_xlsx(df_interaction_sel, paste0(path,"PosInteraction.xlsx"))

df_interaction_p<-df_interaction %>% 
  filter(adj.P.Val<0.8) %>% 
  select(-Reactome.Pathways) %>% 
  arrange(adj.P.Val)
df_interaction_p

df_interaction_sel

sort(p.adjust(df_interaction$P.Value, method = "bonferroni"))


selected_metabolite="P05386 "

################################################################
make_boxplot<-function(x, 
                       selected_metabolite="sucrose",
                       group, pair) {
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
    aes(x=Group, y=Values) +
    geom_boxplot(outlier.colour=NA,outlier.shape = NA, coef=0, 
                 fill="grey", colour="grey20", alpha=0.5) + 
    {if (!missing(pair)) {
                     geom_line(aes(x  = as.numeric(Group)+r, y = Values,
                  group = pData(x)[[pair]]), linetype = "dashed")}} +
    geom_point(aes(x  = as.numeric(Group)+r, fill=Group),
               size=4, shape=21, colour="grey20") +
    labs(title=selected_metabolite, y="Log2 normalised intensity")+
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

levels(pData(eset)$Group)
levels(pData(eset)$Group)<- c("Ctrl AOAA","Ctrl UT", "DS AOAA", "DS UT")
pData(eset)$Group <- factor(pData(eset)$Group, levels = c( "Ctrl UT"  , "Ctrl AOAA","DS UT" ,    "DS AOAA" ))  

rownames(exprs(eset))

names(pData(eset))
make_boxplot(eset,selected_metabolite = "P05386",
             group = "Group", pair = "Cell.Line")
make_boxplot(eset,"Q9NX55",
             group = "Group", pair = "Cell.Line")
make_boxplot(eset,df_interaction_sel$Id[6],
             group = "Group", pair = "Cell.Line")

table(pData(eset)$Cell.Line)

eset_norm_UT<-eset_norm[,eset_norm$Group %in% c("DS UT", "Ctrl UT")]
pData(eset_norm_UT)$Group <- factor(pData(eset_norm_UT)$Group)

make_boxplot(eset_norm_UT,selected_metabolite = "carnosine",
             group = "Group")

make_boxplot(eset_norm,"carnosine",
             group = "Group", pair = "CELL LINE")


make_boxplot(eset_norm,"succinate",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"fumarate",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"malate",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"citrate",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"arginine",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"tiglyl carnitine (C5)",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"tyrosine",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"cysteine",
             group = "Group", pair = "CELL LINE")
make_boxplot(eset_norm,"taurine",
             group = "Group", pair = "CELL LINE")

make_boxplot("glucose")
make_boxplot("pyruvate")
make_boxplot("putrescine")
