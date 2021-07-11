
library(limma)
library(Biobase)
library(readxl)
library(plyr)
library(impute)
library(imputeLCMD)
library(dplyr)

path="C:\\Users\\Acer 3\\Documents\\CsabaII\\Protein\\"


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
            ~ "DS AOAA")) %>%
  mutate(Type=case_when(grepl("Healthy", Sample.Description, fixed = TRUE)
                   ~ "Healthy",
                  grepl("Trisomy", Sample.Description, fixed = TRUE)
                   ~ "Trisomy"))
pData_t


xData<-read_excel(paste0(path,"Batch_corrected_log_intensities_data.xlsx"),
           sheet = 1, range = "A1:AG5539")

names(xData)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
df_ana<-AnnotationDbi::select(org.Hs.eg.db,
                              keys=xData[[1]],
                              keytype = "UNIPROT",
                              c("SYMBOL", "GENENAME"))

df_ana
dim(xData)
dim(df_ana)
head(df_ana)

#add chr location: is chr21 or not?
library(HGNChelper)
new.hgnc.table<-getCurrentHumanMap()

#Genes
ch21genes<-readRDS(paste0(path,"ch21genes.Rds"))
ch21genes<-checkGeneSymbols(ch21genes, map=new.hgnc.table)[[3]]
ch21genes<-ch21genes[!is.na(ch21genes)]
ch21genes

df_ana$ischr21 <-ifelse(df_ana$SYMBOL %in% ch21genes,
                        "chr21", "")
df_ana_xData <-merge(df_ana, xData, by.x="UNIPROT", by.y="Id")

names(df_ana_xData)
dim(df_ana_xData)
df_ana_xData<-na.omit(df_ana_xData)
dim(df_ana_xData)
df_ana_xData[duplicated(df_ana_xData[[3]]),]
sum(duplicated(df_ana_xData[[3]]))
df_ana_xData<-df_ana_xData[!duplicated(df_ana_xData[[3]]),]

xData_m<-as.matrix(df_ana_xData[5:ncol(df_ana_xData)])
rownames(xData_m)<-df_ana_xData[[3]]

names(df_ana_xData)
fData<-df_ana_xData[c(1,2,4)]
rownames(fData)<-df_ana_xData[[3]]

head(rownames(xData_m))

colnames(xData_m)<-make.names(colnames(xData_m))

colnames(xData_m)
rownames(pData_t)
rownames(pData_t)<-colnames(xData_m)

head(pData_t)
head(fData)

eset <- ExpressionSet(assayData = xData_m,
                      phenoData = AnnotatedDataFrame(pData_t),
                      featureData =AnnotatedDataFrame(fData) )
                                                     
                      
head(exprs(eset))
boxplot(exprs(eset))

saveRDS(eset, file = paste0(path,"Protein_eset.rds"))

###############################################################################

##############################################################################
library(dplyr)
library(readxl)
library(Biobase)
library(org.Hs.eg.db)
library(gplots)
library(HGNChelper)
library(Biobase)
library(writexl)
library(limma)


path="C:\\Users\\Acer 3\\Documents\\CsabaII\\Protein\\"

eset <- readRDS(paste0(path,"Protein_eset.rds"))

names(pData(eset))
pData(eset)$Group <-factor(pData(eset)$Group)
pData(eset)$Gender <-factor(pData(eset)$Gender)

summary(pData(eset))

levels(pData(eset)$Group)
pData(eset)$Group <- factor(pData(eset)$Group, 
              levels = c( "Ctrl UT"  , "Ctrl AOAA","DS UT" ,    "DS AOAA" ))  

rownames(exprs(eset))

names(pData(eset))


# View the distribution of the raw data
plotDensities(eset, legend = FALSE)

#Stuxy design
design <- with(pData(eset),model.matrix(~0 + Group + Age + Gender))

design
colnames(design)[1:4] 
colnames(design)[1:4] <- c("Ctrl_UT", "Ctrl_AOAA", "DS_UT","DS_AOAA")
colSums(design)

# Create a contrasts matrix
cm <- makeContrasts(type_UT = DS_UT - Ctrl_UT,
                    type_T = DS_AOAA - Ctrl_AOAA,
                    type_Ctrl = Ctrl_AOAA - Ctrl_UT,
                    type_DS = DS_AOAA - DS_UT,
                    interaction=(DS_AOAA - DS_UT)-
                      (Ctrl_AOAA - Ctrl_UT),
                    levels = design)

# View the contrasts matrix
cm



# Fit the model
fit <- lmFit(eset, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
fit2 <- eBayes(fit2)
fit2
# Summarize results
results <- decideTests(fit2,adjust.method="none")
df_statsum<-summary(results)
Total<-apply(df_statsum,2,sum)
df_statsum2<-rbind(df_statsum,Total)
colnames(df_statsum2)<-c("DS vs. Ctrl UnTreated","DS vs. Ctrl Treated",
                         "AOAA vs. UnTreated Ctrl","AOAA vs. UnTreated DS",
                         "intercept")

stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
head(stats)

# Plot a histogram of the p-values
hist(stats[, "P.Value"])



#User defined function#############################################################
plot_volcano = function(x,
                        y,
                        main = '',
                        xlab = 'log2 FC',
                        ylab = '-log10 adj p',
                        x_cutoffs = log2(c(.8, 1.2)),
                        las = 1,
                        y_cutoff = -log10(.01),
                        plot_lines = T,
                        up_col = '#9E0142',
                        down_col = '#3288BD',
                        to_color = NULL,    # indexes of points to color - will be automatically colored up/down
                        pt_labels = NULL,   # all labels - will automatically only write the sig ones
                        xlim =c(-2, 2), ylim = NULL){
  
  
  
  
  require(ggplot2) 
  require(ggrepel)
  
  sig_up = which(x > x_cutoffs[2] & y > y_cutoff) 
  sig_down = which(x < x_cutoffs[1] & y > y_cutoff)
  
  if (is.null(to_color)){
    to_color = c(sig_up, sig_down)
  }
  
  pt_colors = rep('grey', length(x))
  
  pt_colors[to_color] = ifelse(x[to_color] > 0, up_col,
                               ifelse(x[to_color] < 0, down_col, 'grey60'))
  
  if (!is.null(pt_labels)){
    to_label = c(sig_up, sig_down)
  }
  
  labels=rep("", length(x))
  
  labels[to_label]<- pt_labels[to_label]
  
  df <- data.frame(Log2FC=x,
                   minLog10p=y,
                   pt_labels=labels,
                   pt_colors=pt_colors)
  
  
  ggp<-  ggplot(df)+
    aes(x=Log2FC,y=minLog10p, label = pt_labels) +
    geom_point(shape = 21, colour = "black", fill = pt_colors, size = 3, stroke = 1) +
    geom_text_repel(color=pt_colors, max.overlaps = 30)+
    labs(title = main,
         x=xlab,
         y=ylab)+
    coord_cartesian(xlim = xlim, ylim = ylim)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if(plot_lines==T)      {
    ggp<-ggp + 
      geom_hline(yintercept = y_cutoff, color="grey60", linetype="dashed")+
      geom_vline(xintercept = x_cutoffs[2], color=up_col, linetype="dashed")+
      geom_vline(xintercept = x_cutoffs[1], color=down_col, linetype="dashed")
  }
  
  ggp<-ggp+
    theme(
      axis.line.x.bottom = element_line(color = 'black'),
      axis.line.y.left   = element_line(color = 'black'),
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      panel.border       = element_blank())
  
  return(ggp) 
}
##################################################################

#Contrast1.......................................................................
DS_to_Ctrl.UT<- topTable(fit2, coef=1, adjust.method = "fdr", number= nrow(fit2), confint = T)
names(DS_to_Ctrl.UT)


saveRDS(DS_to_Ctrl.UT, file = paste0(path,"DStoCTR_UT_Prot.rds"))


plotv1 <-plot_volcano(x = DS_to_Ctrl.UT$logFC, 
             y = -log10(DS_to_Ctrl.UT$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.05),
             x_cutoffs = c(-.5, .5), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'DS compared to Ctrl (Untreated) \nmodel fit with age and sex', 
             pt_labels = rownames(DS_to_Ctrl.UT))

print(plotv1)
ggsave(
  filename="DStoCtrl(Untreated).png",
  plot = plotv1,
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

#Contrast2.......................................................................

DS_to_Ctrl.T<- topTable(fit2, coef=2, adjust.method = "fdr", number= nrow(fit2), confint = T)
names(DS_to_Ctrl.T)

saveRDS(DS_to_Ctrl.T, file = paste0(path,"DStoCTR_T_Prot.rds"))


plotv2<-plot_volcano(x = DS_to_Ctrl.T$logFC, 
             y = -log10(DS_to_Ctrl.T$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.05),
             x_cutoffs = c(-.5, .5), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'DS compared to Ctrl (Treated) \nmodel fit with age and sex', 
             pt_labels = rownames(DS_to_Ctrl.T))

print(plotv2)
ggsave(
  filename="DStoCtrl(Treated).png",
  plot = plotv2,
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)
#Contrast3...............................
AOAA_to_UT.CTR<- topTable(fit2, coef=3, adjust.method = "fdr", number= nrow(fit2), confint = T)
names(AOAA_to_UT.CTR)

saveRDS(AOAA_to_UT.CTR, file = paste0(path,"AOAAtoUT_CTR_Prot.rds"))

plotv3<-plot_volcano(x = AOAA_to_UT.CTR$logFC, 
             y = -log10(AOAA_to_UT.CTR$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.1),
             x_cutoffs = c(-.2, .2), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'AOAA compared to Untreated (CTR) \nmodel fit with age and sex', 
             pt_labels = rownames(AOAA_to_UT.CTR))
print(plotv3)
ggsave(
  filename="AOAAtoUT(Ctrl).png",
  plot = plotv3,
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

#Contrast4...............................

AOAA_to_UT.DS<- topTable(fit2, coef=4, adjust.method = "fdr", number= nrow(fit2), confint = T)
names(AOAA_to_UT.DS)

saveRDS(AOAA_to_UT.DS, file = paste0(path,"AOAAtoUT_DS_Prot.rds"))

plotv4<-plot_volcano(x = AOAA_to_UT.DS$logFC, 
             y = -log10(AOAA_to_UT.DS$P.Value), 
             plot_lines = F,
             y_cutoff = -log10(.1),
             x_cutoffs = c(-.5, .5), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'AOAA compared to Untreated (DS) \nmodel fit with age and sex', 
             pt_labels = rownames(AOAA_to_UT.DS))
print(plotv4)

ggsave(
  filename="AOAAtoUT(DS).png",
  plot = plotv4,
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

#Contrast5...............................
#Positive interaction- AOAA modifies slightly the metabolite level on Ctrl cell lines,
#but it has a stronger effect on DS
#Negative interaction- AOAA modifies the metabolite level on Ctrl and DS cells,
#but with an opposite tendency

interaction<- topTable(fit2, coef=5, adjust.method = "fdr", number= nrow(fit2), confint = T)
names(interaction)

plotv5<-plot_volcano(x = interaction$logFC, 
             y = -log10(interaction$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.1),
             x_cutoffs = c(-.2, .2), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'Interaction \nmodel fit with age and sex', 
             pt_labels = rownames(interaction))

print(plotv5)
ggsave(
  filename="Interaction.png",
  plot = plotv5,
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

#add pathway function###################################
add_pathways <- function (df){
startdf<-fData(eset)[rownames(fData(eset)) %in%  rownames(df),]                    
head(startdf)
head(fData(eset))
startdf<-data.frame(GENENAME=rownames(startdf), startdf)

keytypes(org.Hs.eg.db)
dfA<-AnnotationDbi::select(org.Hs.eg.db,
                             keys=rownames(startdf),
                              keytype = "GENENAME",
                              c("GENENAME","PATH"))
head(dfA)
library(KEGG.db)
xx <- as.list(KEGGPATHNAME2ID)

dfA$PATH[2]
length(dfA$PATH)
pathwaylist<-NULL
for(ind in 1:length(dfA$PATH)){
  pathid<-dfA$PATH[ind]
  if (!is.na(pathid)){
    sellistelement=xx[xx==pathid]
    pathwaylist[ind]=names(sellistelement)
  }else{
    pathwaylist[ind]=NA
  }
}

dfA$PATHWAY=pathwaylist

head(dfA)
names(dfA)

dfB<-dfA %>% 
  group_by(GENENAME) %>% 
  mutate(PATHWAYS = paste0(PATHWAY, collapse = ", ")) %>% 
  distinct(PATHWAYS)

startdf$PATHWAYS=dfB$PATHWAYS
names(startdf)
startdf <- dplyr::select(startdf, "PATHWAYS")
print(sum(!rownames(startdf) %in% rownames(df)))

de <- merge(df, startdf, by=0, all=TRUE)
names(de)
sorteddf<-de[order(de$P.Value),]
names(sorteddf)[names(sorteddf)=="Row.names"]<-"Protein"
return (sorteddf)
}
########################################################
View(head(add_pathways(DS_to_Ctrl.UT)))


write_xlsx(list(`DS vs. Ctrl UnTreated`=
                  add_pathways(DS_to_Ctrl.UT) ,
                `DS vs. Ctrl Treated`=
                  add_pathways(DS_to_Ctrl.T) ,
           `AOAA vs. UnTreated Ctrl`=
             add_pathways(AOAA_to_UT.CTR),
           `AOAA vs. UnTreated DS`=
             add_pathways(AOAA_to_UT.DS),
           interaction=
             add_pathways(interaction),
           `Summary Stistics`=
             data.frame(Regulation=rownames(df_statsum2),df_statsum2)
           ),
           paste0(path,"Proteome_tables.xlsx"))

#Modified function from phenobTest
############

pca_plot <- function(x, group, group2, pair, names, ellipse=FALSE, main='', components= c(1, 2),legend=TRUE) {
  #libs
  require(ggplot2)
  require(ggrepel)
  #  require(gridExtra)
  #error control
  stopifnot(length(components)==2)
  stopifnot(is(x, 'ExpressionSet'))
  if (!missing(group)) {
    stopifnot(group %in% colnames(pData(x)))
  }
  if (!missing(group2)) {
    stopifnot(group2 %in% colnames(pData(x)))
  }
  if (!missing(pair)) {
    stopifnot(pair %in% colnames(pData(x)))
  }
  if (!missing(names)) {
    stopifnot(names %in% colnames(pData(x)))
  }
  stopifnot(is(components, 'numeric') | is(components, 'integer'))
  stopifnot(all(!is.na(components)))
  
  #calculate principal components
  pcdat <- prcomp(t(exprs(x)))
  pc <- data.frame(pcdat$x[, components])
  
  pc.lab <- round(pcdat$sdev^2/sum(pcdat$sdev^2)*100,1)[components]
  pc.lab <- paste('PC', components, ' (', pc.lab, '%)', sep='')
  
  colnames(pc) <- c('pc1', 'pc2')
  
  #get plot limits
  lim <- range(pc)
  d <- dist(lim) *.05
  lim <- c(lim[1]-d, lim[2]+d)  
  #
  #parameters for plot
  add <- c()
  if (!missing(group)) {
    colour <- x[[group]]
    if (legend) legend.position <- 'top' else legend.position <- 'none'
  } else {
    colour <- 'black'
    legend.position <- 'none'
  }
  if (!missing(group2)) {
    shape <- x[[group2]]
  } else {
    shape <- ''
    add <- c(add, 'guides(shape=FALSE)')
  }
  if (!missing(pair)) {
    line <- x[[pair]]
    add <- c(add, "geom_line(aes(group=line), colour='grey')")
  } else {
    line <- ''
  }
  if (!missing(names)) {
    label <- x[[names]]
    # add <- c(add, 'geom_text_repel(aes(label=label), size=3, color=as.numeric(as.factor(label)))')
    add <- c(add, 'geom_text_repel(aes(label=label), size=3)')
  } else {
    label <- ''
  }
  dat <- data.frame(pc, colour, shape, line
                    , label
  )
  #
  #ellipse with 95% CI of the mean of each group
  if (ellipse) {
    df.ell <- data.frame()
    for (i in 1:length(levels(colour))) {
      sel <- colour %in% levels(colour)[i]
      df.ell <- rbind(df.ell, cbind(as.data.frame(ellipse(x=cov(data.frame(pc[sel, ])/sqrt(nrow(pc))),
                                                          centre=colMeans(pc[sel, ]))),
                                    colour=levels(colour)[i], shape=''))
    }
    add <- c(add, 'geom_path(data=df.ell, aes(x=pc1, y=pc2, group=colour))')
  }
  #
  #plot
  tmp <- qplot(x=pc1, y=pc2, data=dat, xlab=pc.lab[1], ylab=pc.lab[2], colour=colour, shape=shape, main=main) + 
    geom_point(shape = 21, colour = "black", aes(fill = colour), size = 3, stroke = 1) + 
    coord_cartesian(xlim=lim, ylim=lim) + 
    theme(legend.position=legend.position) +
    ggtitle(main)
  if (length(add)>0) for (i in 1:length(add)) tmp <- tmp + eval(parse(text=add[i]))
  tmp <-tmp + 
    theme_bw() +
    theme(axis.line=element_line(colour="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_blank())
  tmp <-tmp +
    theme(
      axis.line.x.bottom = element_line(color = 'black'),
      axis.line.y.left   = element_line(color = 'black'),
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      panel.border       = element_blank())+
    theme(legend.title=element_blank())+
    theme(legend.position="bottom")
  #theme(legend.justification=c(1,0), legend.position=c(1,0))
  
  return(tmp)
}

pData(eset)
pca_plot(x=eset,
         group='Group',
         names="Cell.Line",
         components= c(3, 2),main="PCA")

ggsave(
  filename="PCA_Group.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 4,
  height = 4,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

pca_plot(x=eset,
         group='Gender',
         names="Cell.Line",
         main="PCA")

ggsave(
  filename="PCA_Gender.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

pca_plot(x=eset, group='Age', names="Cell.Line", main = "PCA")

ggsave(
  filename="PCA_Age.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

names(pData(eset))
pca_plot(x=eset, group="Protein.Sample.Concentration..µg.µl."  
          , names="Cell.Line", main="PCA")

ggsave(
  filename="PCA_Bradford.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)



pca_plot(x=eset, group='Type', names="Cell.Line", components= c(3, 2),main="PCA")

ggsave(
  filename="PCA_Karyotype.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)
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
  vv<-sub('(.{1,45})(\\s|$)', '\\1\n', selected_metabolite)
  selected_metabolite<-gsub('\n$', '', vv)
  print(selected_metabolite)
  
  ggp<- ggplot(df)+
    aes(x=Group, y=Values, fill=Group) +
    geom_boxplot(outlier.colour=NA,outlier.shape = NA, coef=0, 
                  colour="grey20", alpha=0.5) + 
    {if (!missing(pair)) {
                     geom_line(aes(x  = as.numeric(Group)+r, y = Values,
                  group = pData(x)[[pair]]), linetype = "dashed")}} +
    geom_point(aes(x  = as.numeric(Group)+r, fill=Group),
               size=4, shape=21, colour="grey20") +
    labs(title=paste(subtitle,titleadd),
         subtitle= selected_metabolite,
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
################################################################
make_barplot<-function(x, 
                       selected_metabolite="",
                       subtitle="",
                       titleadd="",
                       group,
                       pair) {
  
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
  r <- runif(nrow(df), -0.15, 0.15)
  
  pt1 = t.test(Values ~ Group,
               data=df[df$Group %in% c("Ctrl AOAA", "Ctrl UT"), ],paired=T)
  print(pt1$p.value)
  pt2 = t.test(Values ~ Group,
               data=df[df$Group %in% c("Ctrl UT", "DS UT"), ],paired=FALSE)
  print(pt2$p.value)
  pt3 = t.test(Values ~ Group,
               data=df[df$Group %in% c("DS UT", "DS AOAA"), ],paired=T)
  print(pt3$p.value)
  
  # names(df)
  # df$fourcolor<-mapvalues(df$Group,
  #                             from= c("Ctrl AOAA", "Ctrl UT",  "DS AOAA",   "DS UT" ),
  #                             to=c("blue","green", "red","orange"))
  # print(df)
  vv<-sub('(.{1,25})(\\s|$)', '\\1\n', selected_metabolite)
  selected_metabolite<-gsub('\n$', '', vv)
  print(selected_metabolite)
  
  ggp<- ggplot(df)+
    aes(x=Group, y=Values, fill=Group) +
    stat_summary(fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 fun.min = function(x) mean(x) , 
                 geom = 'errorbar', width=0.3, colour="grey30", size=1)+
    stat_summary(fun= 'mean', geom = 'bar', alpha=0.3,size=1,
                 width=0.8, colour="grey30") +
    
    coord_cartesian(ylim=c(min(Values)-0.01*min(Values),
                           max(Values)+0.02*max(Values)))+
    
    {if (!missing(pair)) {
      geom_line(aes(x  = as.numeric(Group)+r, y = Values,
                    group = pData(x)[[pair]]),
                linetype = "dashed")}} +
    {if (pt1$p.value<0.05) {
      geom_segment(aes(x = 1,
                       y = max(Values)+0.006*max(Values)),
                   xend = 2,
                   yend =  max(Values)+0.006*max(Values))}} +
    {if (pt1$p.value<0.05 & pt1$p.value>=0.01) {
      geom_text(x=1.5, y=max(Values)+0.008*max(Values),
                label="*")}} +
    {if (pt1$p.value<0.01 & pt1$p.value>=0.001) {
      geom_text(x=1.5, y=max(Values)+0.008*max(Values),
                label="* *")}} +
    {if (pt1$p.value<0.001) {
      geom_text(x=1.5, y=max(Values)+0.008*max(Values),
                label="* * *")}} +
    
    {if (pt2$p.value<0.05) {
      geom_segment(aes(x = 1,
                       y = max(Values)+0.016*max(Values)),
                   xend = 3,
                   yend =  max(Values)+0.016*max(Values))}} +
    {if (pt2$p.value<0.05 & pt2$p.value>=0.01) {
      geom_text(x=2, y=max(Values)+0.018*max(Values),
                label="*")}} +
    {if (pt2$p.value<0.01 & pt2$p.value>=0.001) {
      geom_text(x=2, y=max(Values)+0.018*max(Values),
                label="* *")}} +
    {if (pt2$p.value<0.001) {
      geom_text(x=2, y=max(Values)+0.018*max(Values),
                label="* * *")}} +
    
    {if (pt3$p.value<0.05) {
      geom_segment(aes(x = 3,
                       y = max(Values)+0.006*max(Values)),
                   xend = 4,
                   yend =  max(Values)+0.006*max(Values))}} +
    {if (pt3$p.value<0.05 & pt3$p.value>=0.01) {
      geom_text(x=3.5, y=max(Values)+0.008*max(Values),
                label="*")}} +
    {if (pt3$p.value<0.01 & pt1$p.value>=0.001) {
      geom_text(x=3.5, y=max(Values)+0.008*max(Values),
                label="* *")}} +
    {if (pt3$p.value<0.001) {
      geom_text(x=3.5, y=max(Values)+0.008*max(Values),
                label="* * *")}} +
    
    geom_point(aes(x  = as.numeric(Group)+r, fill=Group),
               size=4, shape=21, colour="grey20", stroke=1.1) +
    labs(title=paste(selected_metabolite,titleadd),
       #  subtitle= selected_metabolite,
         y="Log2 normalised intensity")+
    scale_fill_manual(values=c('lightblue','blue', 'orange','red'))+
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
    theme(axis.title.x=element_blank())+
    theme(axis.title.y=element_text(size=10))+
    theme(plot.title = element_text(size=12))
  
  
  ggp <-ggp +
    theme(
      axis.line.x.bottom = element_line(color = 'black'),
      axis.line.y.left   = element_line(color = 'black'),
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      panel.border       = element_blank())+
    theme(legend.title=element_blank())
  return(ggp)
  
}
#########################################################
x11(width = 3.1,
    height = 1.9)
make_barplot(eset,"ADP dependent glucokinase",
             subtitle="Q9BRR6", 
             titleadd = "ischr21"
             , group = "Group", pair = "Cell.Line")

selected_metabolite="ADP dependent glucokinase"


unlink(paste0(paste0(path,"Glyc\\"),
              list.files(paste0(path,"Glyc\\"))))
#Glycolysis
df_ana<-FindUniprotGenename(genes_glyc)
head(df_ana)

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Glyc\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)  
}

library(HGNChelper)
new.hgnc.table<-getCurrentHumanMap()
#Genes
ch21genes<-readRDS(paste0(path,"ch21genes.Rds"))
ch21genes<-checkGeneSymbols(ch21genes, map=new.hgnc.table)[[3]]
ch21genes<-ch21genes[!is.na(ch21genes)]
ch21genes






#Energymetabolites
library(KEGG.db)

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
  df<-df[df$GENENAME %in% rownames(exprs(eset)),]
  df$ischr21 <-ifelse(df$SYMBOL %in% ch21genes, "chr21", "")
  return(df)
}
##################################$

#empty dir
unlink(paste0(paste0(path,"chr21\\"),
              list.files(paste0(path,"chr21\\"))))

unlink(paste0(paste0(path,"Glyc\\"),
              list.files(paste0(path,"Glyc\\"))))

unlink(paste0(paste0(path,"TCA\\"),
              list.files(paste0(path,"TCA\\"))))

unlink(paste0(paste0(path,"FA\\"),
              list.files(paste0(path,"FA\\"))))

unlink(paste0(paste0(path,"PPP\\"),
              list.files(paste0(path,"PPP\\"))))

unlink(paste0(paste0(path,"Glut\\"),
              list.files(paste0(path,"Glut\\"))))

unlink(paste0(paste0(path,"transp\\"),
              list.files(paste0(path,"transp\\"))))
unlink(paste0(paste0(path,"rle\\"),
              list.files(paste0(path,"rle\\"))))

dir.create(paste0(path,"chr21\\"))
dir.create(paste0(path,"transp\\"))
dir.create(paste0(path,"rle\\"))


#Genes on chr21
genes_transp
 
 df_ana<-FindUniprotGenename(ch21genes)
 df_ana
 
 
 for(i in 1:nrow(df_ana)){
   make_barplot(eset,df_ana$GENENAME[i],
                subtitle=df_ana$UNIPROT[i], 
                titleadd = df_ana$ischr21[i],
                group = "Group", pair = "Cell.Line")
   ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
          path=paste0(path,"chr21\\"), 
          scale = 1,
          width = 3.1,
          height = 1.9,
          units = "in",
          dpi = 300,
          limitsize = TRUE)  
 }

#Glycolysis
df_ana<-FindUniprotGenename(genes_glyc)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Glyc\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)  
}

#TCA
df_ana<-FindUniprotGenename(genes_TCA)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"TCA\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  }

#Genes transp

df_ana<-FindUniprotGenename(genes_transp)
df_ana


for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"transp\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)  
}

#FA
df_ana<-FindUniprotGenename(genes_FA)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"FA\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  }

#Glutamine
df_ana<-FindUniprotGenename(genes_Glutamine)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Glut\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  }

#PPP
df_ana<-FindUniprotGenename(genes_PPP)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"PPP\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
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
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Mito1\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)  
}





#Mito2
df_ana<-FindUniprotGenename(complex2)
AnnotationDbi::select(org.Hs.eg.db,keys=complex2, 
                      keytype = "SYMBOL", 
                      c("UNIPROT", "GENENAME"))
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Mito2\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  }

#Mito3
df_ana<-FindUniprotGenename(complex3)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Mito3\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  }

#Mito4
df_ana<-FindUniprotGenename(complex4)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Mito4\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  }

#Mito5
df_ana<-FindUniprotGenename(complex5)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"Mito5\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}

#Rate-limiting



#Rate limiting enzymes
#Glycolysis
#Glucose transporters
rle<-c(
  "SLC2A1" ,  "SLC2A2" ,  "SLC2A3" ,  "SLC2A4", "SLC2A14",
  #hexokinase
  "HK1", "HK2",
  #phosphofructokinase
  "PFKL", "PFKP", "PFKM",
  #Lactate transporters
  "SLC16A1", "SLC16A7", "SLC16A8", "SLC16A3", "SLC16A4",
#Krebcycle
#isocitrate dehydrogenase
"IDH1","IDH2","IDH3A","IDH3B","IDH3G",
#??-ketoglutarate dehydrogenase
"OGDH", "DHTKD1", "OGDHL",
#malate dehydrogenase
"MDH1","MDH2",
#rate limiting enzyme for glycogen synthesis
#glycogen synthase 1
"GYS1", "GYS2",
#rate limiting enzyme for glycogenolysis
"PYGB", "PYGL", "PYGM",
#rate limiting enzyme for the , PPP
"G6PD",
#rate limiting enzyme for de novo pyrimidine synthesis
#carbomoly phosphate synthetase I$$
"CPS1",
#rate limiting enzyme for de novo purine synthesis
"PPAT",
#rate limiting enzyme for fatty acid synthesis
"ACACA", "ACACB",
#fatty acid synthetase
"FASN",
#rate limiting enzyme for ketogenesis, cholesterol synthesis
"HMGCR",
#rate limiting enzyme for catecholamine synthesis
#Tyrosine Hydroxylase
"TH",
#prostaglandibe synthesis
#cyclooxigenases
"PTGS1", "PTGS2",
#alcohol dehydrogenase (ADH) catalyzes
#the rate-limiting step for ethanol metabolism,)
"ADH1A", "ADH1B","ADH1C","ADH4", "ADH5", "ADH6","ADH7",
#11-cis-Retinol dehydrogenase catalyzes the oxidation of cis-retinols, 
#a rate-limiting step in the biosynthesis of 9-cis-retinoic acid.
"RDH5",
#Steroid hormone biosynthesis
#3beta-hydroxy steroid dehydrogenase/isomerase;
"HSD3B1", "HSD3B2",
#Cytochrome c oxidase (COX),
#the rate-limiting enzyme of mitochondrial respiration,
#??-glutamylcysteine synthetase,
#the rate limiting enzyme of chemoprotective glutathione synthesis,
"GCLC"
)



unlink(paste0(paste0(path,"rle\\"),
              list.files(paste0(path,"rle\\"))))

#RLE
df_ana<-FindUniprotGenename(rle)
df_ana

for(i in 1:nrow(df_ana)){
  make_barplot(eset,df_ana$GENENAME[i],
               subtitle=df_ana$UNIPROT[i], 
               titleadd = df_ana$ischr21[i],
               group = "Group", pair = "Cell.Line")
  ggsave(paste0(df_ana$SYMBOL[i],"_",df_ana$UNIPROT[i],".jpg"), 
         path=paste0(path,"rle\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}

## Bimap interface:
x <- org.Hs.egENZYME
head(x)
# Get the entrez gene identifiers that are mapped to an EC number 
mapped_genes <- mappedkeys(x)
head(mapped_genes)
# Convert to a list
xx <- as.list(x[mapped_genes])
xx[xx=="2.3.1.5"]


res2<-unname(res)
res2

pathRTE="C:\\Users\\Acer 3\\Documents\\CsabaII\\Ratelimitingenzimes\\"

RTElist<-as.data.frame(read_excel(paste0(pathRTE,"RTElist.xls")))

mylist <-NULL
for (j in 1:nrow(RTElist) ){
  
  ECname<-gsub("\n","",RTElist[j,1])
  print(ECname)
  
  bb<- as.list(org.Hs.egENZYME2EG)
  geneentrezid<-bb[names(bb)==ECname]
  geneentrezid
  
  res<-AnnotationDbi::select(org.Hs.eg.db, unname(unlist(geneentrezid)),
                             c('SYMBOL',"GENENAME","UNIPROT"), 'ENTREZID')
  
  keytypes(org.Hs.eg.db)
  res
  print(res)
  for (i in 1:nrow(res)){
    res_entrezid=res$ENTREZID[i]
    cc<-as.list(org.Hs.egPATH2EG)
    library(data.table)
    pathcodes<-names(cc)[cc %like% res_entrezid] 
    pathcodes
    xx <- as.list(KEGGPATHNAME2ID)
    pathnames<-names(xx)[which(unname(unlist(xx)) %in% pathcodes)]
    pathnames
    
    res$pathnames[i]<-paste(pathnames, collapse = ', ')
  }
  
  mylist[[j]]<-cbind(ECname,res)
  
}


fulldf<-do.call(rbind.data.frame, mylist)
fulldf

library(writexl)
write_xlsx(fulldf, paste0(pathRTE, "RTEfull.xlsx"))

unname(unlist(cc))


xx <- as.list(KEGGPATHNAME2ID)
xx[1]
xx[which(unname(unlist(xx))=="08630")]

Geneid<-mget(pId, org.Hs.egPATH2EG)

cc
