if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("Biobase")
install.packages("readxl")
BiocManager::install("pcaMethods")
BiocManager::install("impute")

library(limma)
library(Biobase)
library(readxl)
library(plyr)
library(impute)
library(imputeLCMD)

path="C:\\Users\\Acer 3\\Documents\\CsabaII\\MetabolomNEW\\"


#phenotipic data
#pData<-read_excel(paste0(path,"UFRI-01-19VW.xlsx"),
#           sheet = 2, range = "M1:AS22")


pData<-read_excel(paste0(path,"UFRI-01-19VW CO CDT FINAL.xlsx"),
                  sheet = 2, range = "M1:As27")
pData_core <-as.data.frame(t(as.matrix(pData [2:ncol(pData)])))

pData[,1]
colnames(pData_core) <-pData[,1]
colnames(pData_core) <-c("PARENT SAMPLE ID"   ,      "SAMPLE NAME" ,            
                          "AGE"  ,                    "BOX NUMBER"  ,            
                         "CELL LINE"  ,              "CLIENT MATRIX",           
                         "CLIENT SAMPLE ID" ,        "CLIENT SAMPLE NUMBER",    
                         "COMMENTS",                 "DIRECTIONS FOR DUPLICA 1",
                         "DUPLICATES PRESENT" ,      "ESTIMATED PROTEIN SAMP 1",
                          "GENDER"  ,                 "GROUP NAME" ,             
                         "GROUP NUMBER"  ,           "Bradford",
                          "RACE ETHNICITY"   ,        "SAMPLE AMOUNT"  ,         
                         "SAMPLE AMOUNT UNITS" ,     "SAMPLE BOX LOCATION",     
                          "SAMPLE DESCRIPTION"   ,    "TIME POINT"  ,            
                          "TREATMENT"  ,              "UPDATED CLIENT SAMPLE 1", 
                          "VOL EXTRACTED" , "Group"    )

summary(pData_core)

unique(pData_core$AGE)



newAGES=mapvalues(pData_core$AGE, from= c("Newborn", "5 YR",    "5YR",
                                                    "1 YR",    "12 YR",   "3 MO" ,   "2 MO" ,  
                                                    "9 YR" ,   "3 DA"),
             to=c(0,5,5,1,12,0,0,9,0))
as.numeric(as.character(newAGES))
ind=which(colnames(pData_core)=="AGE")
ind
colnames(pData_core)
pData_core[,ind]<-as.numeric(as.character(newAGES))
head(pData_core)
colnames(pData_core)

pData_core$Type<-mapvalues(pData_core$Group,
                           from= c("Ctrl_AOAA", "Ctrl_UT",  "DS_AOAA",   "DS_UT" ),
                           to=c("Healthy", "Healthy",  "Trisomy",   "Trisomy"))

pData_core$Group<-mapvalues(pData_core$Group,
      from= c("Ctrl_AOAA", "Ctrl_UT",  "DS_AOAA",   "DS_UT" ),
      to=c("Ctrl AOAA", "Ctrl UT",  "DS AOAA",   "DS UT"))

#assay data
xData<-as.matrix(read_excel(paste0(path,"UFRI-01-19VW CO CDT FINAL.xlsx"),
                  sheet = 2, range = "N27:AS705" ))
tail(xData)
head(xData)
rownames(pData_core)
colnames(xData) <-rownames(pData_core)
colnames(xData)

fData<-as.data.frame(read_excel(paste0(path,"UFRI-01-19VW CO CDT FINAL.xlsx"),
                            sheet = 2, range = "B27:M705" ))
head(fData)
fData[,1]
rownames(xData)<-fData[,1]
rownames(fData)<-fData[,1]
fData<- fData[,-1] 
head(fData)

esetN <- ExpressionSet(assayData = xData,
            phenoData = AnnotatedDataFrame(as.data.frame(pData_core)),
            featureData = AnnotatedDataFrame(fData))

saveRDS(esetN, file = paste0(path,"Metabolom_esetN.rds"))

###############################################################################

##############################################################################
path="C:\\Users\\Acer 3\\Documents\\CsabaII\\MetabolomNEW\\"


library(limma)
library(Biobase)
library(readxl)
library(plyr)
library(impute)

esetN <- readRDS(paste0(path,"Metabolom_esetN.rds"))

colnames(exprs(esetN))
names(pData(esetN))

#Filter
unique(esetN$Group)
aa<-rowSums(is.na(exprs(esetN[,esetN$Group=="DS UT"])))<5
bb<-rowSums(is.na(exprs(esetN[,esetN$Group=="DS AOAA"])))<5
cc<-rowSums(is.na(exprs(esetN[,esetN$Group=="Ctrl UT"])))<5
dd<-rowSums(is.na(exprs(esetN[,esetN$Group=="Ctrl AOAA"])))<5

length(aa)
sum(aa)


esetN_filtered<-esetN[aa & bb & cc & dd ,]
dim(esetN_filtered)

# eval<-exprs(esetN_filtered)[rownames(exprs(esetN_filtered))=="homocysteine",]
# max(eval)
# plot(eval)
# eval[2]
# eval
# boxplot(eval~pData(esetN)$Group, main="Without normalisation")
# 
# names(eval)
# df_eval<-data.frame(names=names(eval),
#           values=eval,
#           group=pData(esetN)$Group)


pData(esetN_filtered)

# Create new ExpressionSet to store processed data
esetN_norm <- esetN_filtered



#Imputation 
# #install.packages("imputeLCMD") QRILC works only on log2 transformed data
# library(imputeLCMD)
# # perform missing data imputation
# sum(is.na(exprs(esetN_norm)))
# 
# obj.QRILC <- impute.QRILC(log2(exprs(esetN_norm)),tune.sigma = 1)
# #No negative vakue is accepted
# sum(obj.QRILC[[1]]<0)
# View(head(log2(exprs(esetN_norm))))
# View(head(obj.QRILC[[1]]))
# # obj.QRILC[[1]][obj.QRILC[[1]]<0]<-0
# 
# #change back to log2
# exprs(esetN_norm) <- 2^obj.QRILC[[1]]
# 
# sum(is.na(exprs(esetN_norm)))
# 
# exprs(esetN_norm)<0
# sum(exprs(esetN_norm)<0)

#Imputation method B
#impute knn
library(impute)

sum(is.na(exprs(esetN_norm)))

#previously make log2 transformation
m<-log2(exprs(esetN_norm))
sum(impute::impute.knn(m)[[1]]<0)
head(impute::impute.knn(m)[[1]])

plot(m[4,])
plot(impute::impute.knn(m)[[1]][4,])
exprs(esetN_norm) <- impute::impute.knn(m)[[1]]

# #impute smallest values
# 
# dim(exprs(esetN_norm))
# 
# sum(is.na(exprs(esetN_norm)))
# 
# for (i in 1:nrow(exprs(esetN_norm))){
#   exprs(esetN_norm)[i,][is.na(exprs(esetN_norm)[i,])] <- min(exprs(esetN_norm)[i,], na.rm=T)/2
# }
# 
# sum(is.na(exprs(esetN_norm)))

#impute randomForest
# library(randomForest)
# sum(is.na(exprs(esetN_norm)))
# 
# tm<-t(exprs(esetN_norm))
# pData(esetN_norm)$Group
# 
# exprs(esetN_norm)<-t(rfImpute(tm,pData(esetN_norm)$Group))[-1,]
# exprs(esetN_norm)
# View the distribution of the raw data
plotDensities(esetN_norm, legend = FALSE)

# Log tranform
# exprs(esetN_norm) <- log2(exprs(esetN_norm))
# plotDensities(esetN_norm, legend = FALSE)

#Check homocysteine
eval<-exprs(esetN_norm)[rownames(exprs(esetN_norm))=="homocysteine",]
max(eval)
plot(eval)
boxplot(eval~pData(esetN_norm)$Group)


boxplot(exprs(esetN_norm))


# Quantile normalize
exprs(esetN_norm) <- normalizeBetweenArrays(exprs(esetN_norm))
plotDensities(esetN_norm, legend = FALSE)
boxplot(exprs(esetN_norm))

sum(is.na(exprs(esetN_norm)))
# 
# eval<-exprs(esetN_norm)[rownames(exprs(esetN_norm))=="homocysteine",]
# max(eval)
# plot(eval)
# boxplot(exprs(esetN_norm))


eval
pData(esetN)$Group


names(pData(esetN_norm))
pData(esetN_norm)$Group <-factor(pData(esetN_norm)$Group)
pData(esetN_norm)$GENDER <-factor(pData(esetN_norm)$GENDER)
pData(esetN_norm)$AGE <-as.numeric(as.character(pData(esetN_norm)$AGE ))
pData(esetN_norm)$Bradford <-as.numeric(as.character(pData(esetN_norm)$Bradford))


pData(esetN_norm)$AGE

levels(pData(esetN_norm)$Group)
pData(esetN_norm)$Group <- factor(pData(esetN_norm)$Group, levels = c( "Ctrl UT"  , "Ctrl AOAA","DS UT" ,    "DS AOAA" ))  


#Study design
design <- with(pData(esetN_norm),model.matrix(~0 + Group + AGE + GENDER))

design
colnames(design)
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
fit <- lmFit(esetN_norm, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
fit2 <- eBayes(fit2)
fit2
# Summarize results
results <- decideTests(fit2,adjust.method="none")
summary(results)
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


saveRDS(DS_to_Ctrl.UT, file = paste0(path,"DStoCTR_New.rds"))



plot_volcano(x = DS_to_Ctrl.UT$logFC, 
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

ggsave(
  filename="DStoCtrl(Untreated).png",
  plot = last_plot(),
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

plot_volcano(x = DS_to_Ctrl.T$logFC, 
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

ggsave(
  filename="DStoCtrl(Treated).png",
  plot = last_plot(),
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

plot_volcano(x = AOAA_to_UT.CTR$logFC, 
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

ggsave(
  filename="AOAAtoUT(Ctrl).png",
  plot = last_plot(),
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

plot_volcano(x = AOAA_to_UT.DS$logFC, 
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


ggsave(
  filename="AOAAtoUT(DS).png",
  plot = last_plot(),
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

plot_volcano(x = interaction$logFC, 
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

ggsave(
  filename="Interaction.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)


write_xlsx(list(`DS vs. Ctrl UnTreated`=
                  data.frame(Metabolite=rownames(DS_to_Ctrl.UT),DS_to_Ctrl.UT),
                `DS vs. Ctrl Treated`=
                  data.frame(Metabolite=rownames(DS_to_Ctrl.T),DS_to_Ctrl.T),
                `AOAA vs. UnTreated Ctrl`=
                  data.frame(Metabolite=rownames(AOAA_to_UT.CTR),AOAA_to_UT.CTR),
                `AOAA vs. UnTreated DS`=
                  data.frame(Metabolite=rownames(AOAA_to_UT.DS),AOAA_to_UT.DS),
                interaction=
                  data.frame(Metabolite=rownames(interaction),interaction),
                `Summary Stistics`= 
                  data.frame(Regulation=rownames(df_statsum2),df_statsum2)),
                paste0(path,"Metabolome_tables.xlsx"))





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

pca_plot(x=esetN_norm,
         group='Group',
         names="CELL LINE",
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

pca_plot(x=esetN_norm,
         group='GENDER',
         names="CELL LINE",
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

pca_plot(x=esetN_norm, group='AGE', names="CELL LINE", main = "PCA")

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

pca_plot(x=esetN_norm, group="Bradford" , names="CELL LINE", main="PCA")

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



pca_plot(x=esetN_norm, group='Type', names="CELL LINE", components= c(3, 2),main="PCA")

ggsave(
  filename="PCA_Karyotype.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 4,
  height = 3,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)
names(pData(esetN))
# 
# esetN_norm_UT<-esetN_norm[,esetN_norm$Group %in% c("DS_UT", "Ctrl_UT")]
# pca_plot(x=esetN_norm_UT, group='Group', names="CELL LINE")
# pca_plot(x=esetN_norm_UT, group='AGE', names="CELL LINE")
# 
# esetN_norm_DS<-esetN_norm[,esetN_norm$Group %in% c("DS_UT", "DS_AOAA")]
# pca_plot(x=esetN_norm_DS, group='Group', names="CELL LINE")
# 

#heatmap




selected_metabolite="sucrose"

################################################################
make_boxplot<-function(x, 
                       selected_metabolite="",
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
 print(df) 
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
    labs(title=selected_metabolite,
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



################################################################
make_barplot<-function(x, 
                       selected_metabolite="",
                       group, pair) {
  
  # x=esetN_norm
  # selected_metabolite="homocysteine"
  # group = "Group"
  # pair = "CELL LINE"  
  
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
              data=df[df$Group %in% c("Ctrl UT", "DS UT"), ], paired=FALSE)
 print(pt2$p.value)
 pt3 = t.test(Values ~ Group,
              data=df[df$Group %in% c("DS UT", "DS AOAA"), ],paired=T)
 print(pt3$p.value)

  # names(df)
  # df$fourcolor<-mapvalues(df$Group,
  #                             from= c("Ctrl AOAA", "Ctrl UT",  "DS AOAA",   "DS UT" ),
  #                             to=c("blue","green", "red","orange"))
  # print(df)
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
    {if (pt3$p.value<0.01 & pt3$p.value>=0.001) {
      geom_text(x=3.5, y=max(Values)+0.008*max(Values),
                label="* *")}} +
    {if (pt3$p.value<0.001) {
      geom_text(x=3.5, y=max(Values)+0.008*max(Values),
                label="* * *")}} +
    
     geom_point(aes(x  = as.numeric(Group)+r, fill=Group),
               size=4, shape=21, colour="grey20", stroke=1.1) +
    labs(title=selected_metabolite,
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
    theme(axis.title.x=element_blank())
  
  
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
x11(height = 2, width=3)
make_barplot(esetN_norm,"homocysteine",
             group = "Group", pair = "CELL LINE")
pData(esetN_norm)

should_be_increased <- c("urate", "lactate", "succinate",
                         "cysteine", "cystine","kynurenine","phosphate",
                         "pyruvate","AMP","cystathionine",
                         "tiglyl carnitine (C5)","2-hydroxyglutarate")
dir.create(paste0(path,"should_be_increased\\"))
#delete content
unlink(paste0(paste0(path,"should_be_increased\\"),
              list.files(paste0(path,"should_be_increased\\"))))
for( i in 1:length(should_be_increased)){
  make_barplot(esetN_norm,selected_metabolite = should_be_increased[i],
               group = "Group", pair = "CELL LINE")
  ggsave(paste0(make.names(should_be_increased[i]),".jpg"), 
         path=paste0(path,"should_be_increased\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE) 
}

rownames(exprs(esetN_norm))

#should be decreased


should_be_decreased <- c("serine", "tyrosine", "arginine",
                         "homocysteine", "ATP","choline",
                         "carnitine","NAD+",
                         "glutathione, reduced (GSH)",
                         "glutathione, oxidized (GSSG)")
dir.create(paste0(path,"should_be_decreased\\"))
#delete content
unlink(paste0(paste0(path,"should_be_decreased\\"),
              list.files(paste0(path,"should_be_decreased\\"))))
for( i in 1:length(should_be_decreased)){
  make_barplot(esetN_norm,selected_metabolite = should_be_decreased[i],
               group = "Group", pair = "CELL LINE")
  ggsave(paste0(make.names(should_be_decreased[i]),".jpg"), 
         path=paste0(path,"should_be_decreased\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE) 
}

effect_AOAA <- c("citrate", "malate", "phosphate",
                         "acetyl-CoA", "succinate",
                 "dimethylglycine","homocysteine",
                 "hypotaurine","N-acetyltaurine")
effect_AOAA
dir.create(paste0(path,"effect_AOAA\\"))
#delete content
unlink(paste0(paste0(path,"effect_AOAA\\"),
              list.files(paste0(path,"effect_AOAA\\"))))
for( i in 1:length(effect_AOAA)){
  make_barplot(esetN_norm,selected_metabolite = effect_AOAA[i],
               group = "Group", pair = "CELL LINE")
  ggsave(paste0(make.names(effect_AOAA[i]),".jpg"), 
         path=paste0(path,"effect_AOAA\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE) 
}





# esetN_norm_UT<-esetN_norm[,esetN_norm$Group %in% c("DS UT", "Ctrl UT")]
# pData(esetN_norm_UT)$Group <- factor(pData(esetN_norm_UT)$Group)
# 
# make_barplot(esetN_norm_UT,selected_metabolite = "carnosine",
#              group = "Group")

names(fData(esetN_norm))
superpathways<-unique(fData(esetN_norm)$`SUPER PATHWAY`)
superpathways
subpathways<-unique(fData(esetN_norm)$`SUB PATHWAY`)
subpathways


make_barplot(esetN_norm,"carnosine",
             group = "Group", pair = "CELL LINE")


make_barplot(esetN_norm,"succinate",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"fumarate",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"malate",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"citrate",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"arginine",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"tiglyl carnitine (C5)",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"tyrosine",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"cysteine",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"taurine",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"2-hydroxypalmitate",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"glutathione, reduced (GSH)",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"glutathione, oxidized (GSSG)",
             group = "Group", pair = "CELL LINE")

#In a resting cell, the molar GSH:GSSG ratio exceeds 100:1,
#while in various models of oxidative stress,
#this ratio has been demonstrated to decrease
#to values of 10:1 and even 1:1 (13)

x=esetN_norm
m_gsh <-exprs(x)[rownames(exprs(x)) %in% c("glutathione, reduced (GSH)",
                                                          "glutathione, oxidized (GSSG)" ),]
m_gshpergssg <-t(as.matrix((2^m_gsh[2,])/(2^m_gsh[1,])))
View(m_gshpergssg)
dim(m_gshpergssg)
rownames(m_gshpergssg)<-"GSHperGSSH"
eset_gshpergssg <- ExpressionSet(assayData = m_gshpergssg)
pData(eset_gshpergssg) <- pData(esetN_norm)                              

make_barplot(eset_gshpergssg,"GSHperGSSH",
             group = "Group", pair = "CELL LINE")                                 

make_barplot(esetN,"homocysteine",
             group = "Group")

eval<-exprs(esetN_norm)[rownames(exprs(esetN_norm))=="homocysteine",]
max(eval)
eval

eval<-exprs(esetN)[rownames(exprs(esetN))=="homocysteine",]
max(eval)
eval[2]
plot(eval)

eval<-exprs(esetN_filtered)[rownames(exprs(esetN_filtered))=="homocysteine",]
max(eval)
plot(eval)
length(eval)
eval[2]



eval<-exprs(esetN_norm)[rownames(exprs(esetN_norm))=="homocysteine",]
max(eval)
plot(eval)
eval[2]
eval<-exprs(esetN_norm)[rownames(exprs(esetN_norm))=="homocysteine",]
max(eval)
plot(eval)
make_barplot(esetN_norm,"glucose",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"glucose",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"pyruvate",
             group = "Group", pair = "CELL LINE")
make_barplot(esetN_norm,"putrescine",
             group = "Group", pair = "CELL LINE")

unique((fData(esetN_norm))[[2]])

dir.create(paste0(path,"Glyc\\"))
glyc_molecules<-rownames(fData(esetN_norm))[(fData(esetN_norm))[[2]]=="Glycolysis, Gluconeogenesis, and Pyruvate Metabolism"]
glyc_molecules
unlink(paste0(paste0(path,"Glyc\\"),
              list.files(paste0(path,"Glyc\\"))))

for(i in 1:length(glyc_molecules)) {
  make_barplot(esetN_norm,glyc_molecules[i],
               group = "Group", pair = "CELL LINE")
  ggsave(paste0(make.names(glyc_molecules[i]),".jpg"), 
         path=paste0(path,"Glyc\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE) 
  }

dir.create(paste0(path,"TCA\\"))
tca_molecules<-rownames(fData(esetN_norm))[(fData(esetN_norm))[[2]]=="TCA Cycle"]

for(i in 1:length(tca_molecules)) {
  make_barplot(esetN_norm,tca_molecules[i],
               group = "Group", pair = "CELL LINE")
  ggsave(paste0(make.names(tca_molecules[i]),".jpg"), 
         path=paste0(path,"Glyc\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
}

"Methionine, Cysteine, SAM and Taurine Metabolism"
met_molecules<-rownames(fData(esetN_norm))[(fData(esetN_norm))[[2]]=="Methionine, Cysteine, SAM and Taurine Metabolism"]

for(i in 1:22) {
  make_barplot(esetN_norm,met_molecules[i],
               group = "Group", pair = "CELL LINE")
  ggsave(paste0(make.names(met_molecules[i]),".jpg"), 
         path=paste0(path,"Glyc\\"), 
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}

energy_molecules<-rownames(fData(esetN_norm))[(fData(esetN_norm))[[1]]=="Energy"]
energy_molecules

lipid_molecules <-rownames(fData(esetN_norm))[(fData(esetN_norm))[[1]]=="Lipid"]
lipid_molecules

dir.create(paste0(path,"Lipid\\"))
unlink(paste0(paste0(path,"Lipid\\"),
              list.files(paste0(path,"Lipid\\"))))

for(i in 1:length(lipid_molecules)) {
  make_barplot(esetN_norm,lipid_molecules[i],
               group = "Group", pair = "CELL LINE")
  ggsave(paste0(make.names(lipid_molecules[i]),".jpg"), 
         path=paste0(path,"Lipid\\"),
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}

#########################################################

make_grupedbarplot<-function(x, 
                             compounds="",
                             group) {


  require(ggplot2)
  require (tidyr)
  stopifnot(is(x, 'ExpressionSet'))
  if (!missing(compounds)) {
    stopifnot(compounds %in% rownames(exprs(x)))
  }
  if (!missing(group)) {
    stopifnot(group %in% colnames(pData(x)))
  }
  
  
  
  df <-data.frame(exprs(x)[rownames(exprs(x)) %in% compounds,]
                  # ,Group=pData(x)[[group]]
  )
  colnames(df) <-pData(x)[[group]]
  colnames(df) %in% c("Ctrl UT", "DS UT") 
  df$Compound= factor(rownames(df))
  df2<-pivot_longer(df,!Compound, names_to = "Groups", values_to = "values")
  df2$Groups <- factor(df2$Groups,
                       levels = c("Ctrl UT", "Ctrl AOAA" , "DS UT","DS AOAA"))
  df2<-df2[df2$Groups %in% c("Ctrl UT", "DS UT"),]
  
  ggp<- ggplot(df2)+
    aes(x=Compound, y=values, fill=Groups, group=Groups) +
    stat_summary(fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 fun.min = function(x) mean(x) , 
                 geom = 'errorbar',position = position_dodge(0.8),width=0.3, colour="grey30", size=1)+
    stat_summary(fun= 'mean', geom = 'bar',position = "dodge", alpha=0.3,size=1,
                 width=0.8, colour="grey30")+
    coord_cartesian(ylim=c(min(df2$values)-0.01*min(df2$values),
                           max(df2$values)+0.02*max(df2$values)))+
    scale_fill_manual(values=c('lightblue', 'orange'))+
    theme_bw() +
    theme(axis.line=element_line(colour="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_blank())+
    theme(axis.text.x = element_text( color="black", 
                                      size=10),
          axis.text.y = element_text( color="black", 
                                      size=10))+
    #theme(axis.text.x = element_text(angle = 45, vjust=0.5))+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    theme(legend.position = "bottom")+
    theme(axis.title.x=element_blank())+
    labs(y="Log2 normalised intensity")
  
  ggp
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

make_grupedbarplot(esetN_norm,
                   c("kynurenine","cystathionine","arginine",
                     "serine","tyrosine","carnitine","NAD+"),
             group = "Group")
pData(esetN_norm)

dir.create(paste0(path,"Groupbar"))
ggsave("Differences.jpg", 
       path=paste0(path,"Groupbar\\"),
       scale = 1,
       width = 5,
       height = 3,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
