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

path="C:\\Users\\Acer 3\\Documents\\CsabaII\\MetabolomOLD\\"


#phenotipic data
pData<-read_excel(paste0(path,"UFRI-01-19VW.xlsx"),
           sheet = 2, range = "M1:AS22")
pData_core <-t(as.matrix(pData [2:ncol(pData)]))

colnames(pData_core) <-as.data.frame(pData)[,1]
colnames(pData_core) <-c("SAMPLE NAME" ,            "CONTROL"   ,             
                          "AGE"       ,              "BOX NUMBER"  ,           
                          "CELL LINE"   ,            "CLIENT SAMPLE ID" ,      
                          "CLIENT SAMPLE NUMBER" ,   "COMMENTS"  ,             
                          "GENDER"  ,                "GROUP NUMBER" ,          
                          "ORIGINAL MANIFEST LABEL", "BRADFORD MG ML" ,        
                          "RACE ETHNICITY"  ,        "SAMPLE AMOUNT" ,         
                          "SAMPLE AMOUNT UNITS" ,    "SAMPLE BOX LOCATION"  ,  
                          "SAMPLE DESCRIPTION" ,     "TIME POINT" ,            
                         "TREATMENT"   ,            "VOLUME EXTRACTED" ,      
                          "Group HMDB"        )

summary(pData_core)

unique(as.data.frame(pData_core)$AGE)


newAGES=mapvalues(as.data.frame(pData_core)$AGE, from= c("Newborn", "5 YR",    "5YR",
                                                    "1 YR",    "12 YR",   "3 MO" ,   "2 MO" ,  
                                                    "9 YR" ,   "3 DA"),
             to=c(0,5,5,1,12,0,0,9,0))
as.numeric(as.character(newAGES))
ind=which(colnames(pData_core)=="AGE")
ind
colnames(pData_core)
pData_core[,ind]<-as.numeric(as.character(newAGES))
head(pData_core)
#assay data
xData<-as.matrix(read_excel(paste0(path,"UFRI-01-19VW.xlsx"),
                  sheet = 2, range = "N22:AS537" ))
colnames(xData) <-rownames(pData_core)
colnames(xData)

fData<-as.data.frame(read_excel(paste0(path,"UFRI-01-19VW.xlsx"),
                            sheet = 2, range = "B22:M537" ))
head(fData)
fData[,1]
rownames(xData)<-fData[,1]
rownames(fData)<-fData[,1]
fData<- fData[,-1] 
head(fData)

esetO <- ExpressionSet(assayData = xData,
            phenoData = AnnotatedDataFrame(as.data.frame(pData_core)),
            featureData = AnnotatedDataFrame(fData))

saveRDS(esetO, file = paste0(path,"Metabolom_esetO.rds"))

###############################################################################

##############################################################################
path="C:\\Users\\Acer 3\\Documents\\CsabaII\\MetabolomOLD\\"

library(limma)
library(Biobase)
library(readxl)
library(plyr)
library(impute)

esetO <- readRDS(paste0(path,"Metabolom_esetO.rds"))


names(pData(esetO))[names(pData(esetO))=="Group HMDB"]<-"Group"

#Filter
unique(esetO$Group)
aa<-rowSums(is.na(exprs(esetO[,esetO$Group=="DS_UT"])))<5
bb<-rowSums(is.na(exprs(esetO[,esetO$Group=="DS_AOAA"])))<5
cc<-rowSums(is.na(exprs(esetO[,esetO$Group=="Ctrl_UT"])))<5
dd<-rowSums(is.na(exprs(esetO[,esetO$Group=="Ctrl_AOAA"])))<5


esetO_filtered<-esetO[aa & bb & cc & dd ,]
esetO_filtered 


# Create new ExpressionSet to store processed data
esetO_norm <- esetO_filtered



# #install.packages("imputeLCMD")
# library(imputeLCMD)
# # perform missing data imputation
# sum(is.na(exprs(esetO_norm)))
# 
# obj.QRILC <- impute.QRILC(exprs(esetO_norm),tune.sigma = 1)
# sum(obj.QRILC[[1]]<0)
# 
#  obj.QRILC[[1]][obj.QRILC[[1]]<0]<-0
# exprs(esetO_norm) <- obj.QRILC[[1]]
# 
# exprs(esetO_norm)<0
# sum(exprs(esetO_norm)<0)
#impute knn
# library(impute)
# 
# sum(is.na(exprs(esetO_norm)))
# 
# m<-(exprs(esetO_norm))
# 
# 
# exprs(esetO_norm) <- impute::impute.knn(m)[[1]]
# 
# sum(is.na(exprs(esetO_norm)))


dim(exprs(esetN_norm))

sum(is.na(exprs(esetN_norm)))

for (i in 1:nrow(exprs(esetN_norm))){
  exprs(esetN_norm)[i,][is.na(exprs(esetN_norm)[i,])] <- min(exprs(esetN_norm)[i,], na.rm=T)
}

sum(is.na(exprs(esetN_norm)))

#impute randomForest
# library(randomForest)
# sum(is.na(exprs(esetO_norm)))
# 
# tm<-t(exprs(esetO_norm))
# pData(esetO_norm)$Group
# 
# exprs(esetO_norm)<-t(rfImpute(tm,pData(esetO_norm)$Group))[-1,]
# exprs(esetO_norm)
# View the distribution of the raw data
plotDensities(esetO_norm, legend = FALSE)

# Log tranform
exprs(esetO_norm) <- log2(exprs(esetO_norm))
plotDensities(esetO_norm, legend = FALSE)

# Quantile normalize
exprs(esetO_norm) <- normalizeBetweenArrays(exprs(esetO_norm))
plotDensities(esetO_norm, legend = FALSE)


names(pData(esetO_norm))
pData(esetO_norm)$Group <-factor(pData(esetO_norm)$Group)
pData(esetO_norm)$GENDER <-factor(pData(esetO_norm)$GENDER)
pData(esetO_norm)$AGE <-as.numeric(as.character(pData(esetO_norm)$AGE ))
pData(esetO_norm)$`BRADFORD MG ML` <-as.numeric(as.character(pData(esetO_norm)$`BRADFORD MG ML`))


pData(esetO_norm)$AGE

#Stuxy design
design <- with(pData(esetO_norm),model.matrix(~0 + Group + AGE + GENDER))

colnames(design)
colnames(design)[1:4] <- c("Ctrl_AOAA", "Ctrl_UT", "DS_AOAA","DS_UT")
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
fitO <- lmFit(esetO_norm, design)

# Fit the contrasts
fitO2 <- contrasts.fit(fitO, contrasts = cm)

# Calculate the t-statistics for the contrasts
fitO2 <- eBayes(fitO2)
fitO2
# Summarize results
results <- decideTests(fitO2)
summary(results)

stats <- topTable(fitO2, number = nrow(fitO2), sort.by = "none")
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
    geom_text_repel(color=pt_colors)+
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
  
  return(ggp) 
}
##################################################################

#Contrast1.......................................................................
DS_to_Ctrl.UT.O<- topTable(fitO2, coef=1, adjust.method = "fdr", number= nrow(fitO2), confint = T)
names(DS_to_Ctrl.UT.O)

saveRDS(DS_to_Ctrl.UT.O, file = paste0(path,"DStoCTR_Old.rds"))


plot_volcano(x = DS_to_Ctrl.UT.O$logFC, 
             y = -log10(DS_to_Ctrl.UT.O$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.05),
             x_cutoffs = c(-.5, .5), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'DS compared to Ctrl (Untreated) \nmodel fit with age and sex', 
             pt_labels = rownames(DS_to_Ctrl.UT.O))

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

DS_to_Ctrl.T.O<- topTable(fitO2, coef=2, adjust.method = "fdr", number= nrow(fitO2), confint = T)
names(DS_to_Ctrl.T.O)

plot_volcano(x = DS_to_Ctrl.T.O$logFC, 
             y = -log10(DS_to_Ctrl.T.O$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.05),
             x_cutoffs = c(-.5, .5), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'DS compared to Ctrl (Treated) \nmodel fit with age and sex', 
             pt_labels = rownames(DS_to_Ctrl.T.O))

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
AOAA_to_UT.CTR.O<- topTable(fitO2, coef=3, adjust.method = "fdr", number= nrow(fitO2), confint = T)
names(AOAA_to_UT.CTR.O)

plot_volcano(x = AOAA_to_UT.CTR.O$logFC, 
             y = -log10(AOAA_to_UT.CTR.O$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.1),
             x_cutoffs = c(-.2, .2), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'AOAA compared to Untreated (CTR) \nmodel fit with age and sex', 
             pt_labels = rownames(AOAA_to_UT.CTR.O))

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

AOAA_to_UT.DS.O<- topTable(fitO2, coef=4, adjust.method = "fdr", number= nrow(fitO2), confint = T)
names(AOAA_to_UT.DS.O)

plot_volcano(x = AOAA_to_UT.DS.O$logFC, 
             y = -log10(AOAA_to_UT.DS.O$P.Value), 
             plot_lines = F,
             y_cutoff = -log10(.1),
             x_cutoffs = c(-.5, .5), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'AOAA compared to Untreated (DS) \nmodel fit with age and sex', 
             pt_labels = rownames(AOAA_to_UT.DS.O))


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

interaction.O<- topTable(fitO2, coef=5, adjust.method = "fdr", number= nrow(fitO2), confint = T)
names(interaction.O)

plot_volcano(x = interaction.O$logFC, 
             y = -log10(interaction.0$P.Value), 
             plot_lines = T,
             y_cutoff = -log10(.1),
             x_cutoffs = c(-.2, .2), xlim = c(-2,2),
             xlab = 'log2 FC',
             ylab = '-log10 p',
             up_col = 'red',
             down_col = 'blue',
             main = 'Interaction \nmodel fit with age and sex', 
             pt_labels = rownames(interaction.O))

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


#Modified function from phenobTest
############

pca_plot <- function(x, group, group2, pair
          , names, ellipse=FALSE, main=''
          , components= c(1, 2),legend=TRUE) {
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
  tmp
}

pca_plot(x=esetO_norm, group='Group'
         , names="CELL LINE", components= c(3, 2)
         ,main="PCA")

ggsave(
  filename="PCA_Group.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

pca_plot(x=esetO_norm, group='GENDER'
         , names="CELL LINE", main="PCA")

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

pca_plot(x=esetO_norm, group='AGE', 
         names="CELL LINE", main = "PCA")

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

pca_plot(x=esetO_norm, group="BRADFORD MG ML" ,
         names="CELL LINE", main="PCA")

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

names(pData(esetO))

esetO_norm_UT<-esetO_norm[,esetO_norm$Group %in% c("DS_UT", "Ctrl_UT")]
pca_plot(x=esetO_norm_UT, group='Group', names="CELL LINE")
pca_plot(x=esetO_norm_UT, group='AGE', names="CELL LINE")

esetO_norm_DS<-esetO_norm[,esetO_norm$Group %in% c("DS_UT", "DS_AOAA")]
pca_plot(x=esetO_norm_DS, group='Group', names="CELL LINE")


selected_metabolite="sucrose"

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

levels(pData(esetO_norm)$Group)
levels(pData(esetO_norm)$Group)<- c("Ctrl AOAA","Ctrl UT"
                                    , "DS AOAA", "DS UT")
pData(esetO_norm)$Group <- factor(pData(esetO_norm)$Group, levels = c( "Ctrl UT"  , "Ctrl AOAA","DS UT" ,    "DS AOAA" ))  

rownames(exprs(esetO_norm))
make_boxplot(esetO_norm,selected_metabolite = "lactate",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"benzoate",
             group = "Group", pair = "CELL LINE")

esetO_norm_UT<-esetO_norm[,esetO_norm$Group %in% c("DS UT", "Ctrl UT")]
pData(esetO_norm_UT)$Group <- factor(pData(esetO_norm_UT)$Group)

make_boxplot(esetO_norm_UT,selected_metabolite = "carnosine",
             group = "Group")

make_boxplot(esetO_norm,"carnosine",
             group = "Group", pair = "CELL LINE")


make_boxplot(esetO_norm,"succinate",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"fumarate",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"malate",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"citrate",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"arginine",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"tiglyl carnitine (C5)",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"tyrosine",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"cysteine",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm,"taurine",
             group = "Group", pair = "CELL LINE")

make_boxplot(esetO_norm,"glucose"
             ,group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm, "pyruvate",
             group = "Group", pair = "CELL LINE")
make_boxplot(esetO_norm, "putrescine",
             group = "Group", pair = "CELL LINE")
