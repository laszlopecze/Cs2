
library(limma)
library(Biobase)
library(readxl)
library(plyr)
library(impute)
library(imputeLCMD)
library(dplyr)
library(tidyr)

path="C:\\Users\\Acer 3\\Documents\\CsabaII\\Fluxomics\\"


#phenotipic data

pData<-read_excel(paste0(path,"Flux_results.xls"),
                  sheet = 1, range = "E2:AJ6")
pData<-as.matrix(pData)
colnames(pData)<-make.names(colnames(pData))
pData <-t(pData)
pData<-as.data.frame(pData)
pData
names(pData)<-c("Age","Sex","Treatment","Cell_line")

unique(pData$Age)
pData$Age=mapvalues(pData$Age, from= c("Newborn", "5 YR",    
                                                    "1 YR",    "12 YR",   "3 MO" ,   "2 MO" ,  
                                                    "9 YR" ,   "3 DA"),
                  to=c(0,5,1,12,0,0,9,0))
unique(pData$Cell_line)

DS_cell_lines=c("AG07096",
                "AG05397",
                "3-FCYPR10000285#",
                "Detroit 532",
                "GM02571",
                "GM004616",
                "3-FCYPR10000369*",
                "Detroit 539")

pData$Type <-ifelse(pData$Cell_line %in% DS_cell_lines,"Trisomy","Healthy")
pData
pData_t<-pData %>% 
  mutate(Group=
  case_when(Type=="Healthy" 
            & grepl("Untreated", Treatment, fixed = TRUE)
            ~ "Ctrl UT",
            Type == "Healthy" 
            & grepl("Treated", Treatment, fixed = TRUE)
            ~ "Ctrl AOAA",
            Type == "Trisomy" 
            & grepl("Untreated", Treatment, fixed = TRUE)
            ~ "DS UT",
            Type == "Trisomy" 
            & grepl("Treated", Treatment, fixed = TRUE)
            ~ "DS AOAA")) %>% 
  mutate(AOAA_Treatment=
           case_when(grepl("Untreated", Treatment, fixed = TRUE)
                     ~ "Untreated",
                     grepl("Treated", Treatment, fixed = TRUE)
                     ~ "Treated",
                     grepl("Untreated", Treatment, fixed = TRUE)
                     ~ "Untreated",
                     grepl("Treated", Treatment, fixed = TRUE)
                     ~ "Treated"))
pData_t$Group
pData_t

fxData<-read_excel(paste0(path,"Flux_results.xls"),
           sheet = 1, range = "A8:AJ467")
fxData<-as.data.frame(fxData)
names(fxData)<-make.names(names(fxData))
names(fxData)[names(fxData)=="radioisotope.compound.in.analysis"] <-"Compound"
names(fxData)[names(fxData)=="PATHWAY"] <-"Pathway"
names(fxData)[names(fxData)=="Metabolite.Flux.ID.number"] <-"Flux.ID"

names(fxData)

#remove duplicates
fxData<-fxData[!duplicated(fxData[[3]]),]

xData <-as.matrix(fxData[5:36])
fData <-as.data.frame(fxData[1:4])
fData$Compound<-gsub("\\+", " m\\+",fData$Compound)
fData$Compound <-gsub("  ", " ",fData$Compound)
fData$Compound <-gsub("_", " ",fData$Compound)
fData$Compound
rownames(xData)<-fData$Compound
rownames(fData)<-fData$Compound

fData$Pathway[fData$Pathway %in% as.character(1:50)]<-NA

fData<-fill(fData,c(Pathway,Flux.ID))


colnames(xData)
rownames(pData_t)
identical(rownames(pData_t),colnames(xData))

head(pData_t)

eset <- ExpressionSet(assayData = xData,
                      phenoData = AnnotatedDataFrame(pData_t),
                      featureData=AnnotatedDataFrame(fData))
                                                     
                      
head(exprs(eset))
boxplot(exprs(eset))

saveRDS(eset, file = paste0(path,"Flux_eset.rds"))

###############################################################################

##############################################################################
path="C:\\Users\\Acer 3\\Documents\\CsabaII\\Fluxomics\\"

#install.packages("sjPlot")
library(sjPlot)
library(limma)
library(Biobase)
library(readxl)
library(plyr)
library(impute)

esetN <- readRDS(paste0(path,"Flux_eset.rds"))

colnames(exprs(esetN))
names(pData(esetN))
rownames(exprs(esetN))
exprs(esetN)[exprs(esetN)== 0]<-NA

#Filter
unique(esetN$Group)
aa<-rowSums(is.na(exprs(esetN[,esetN$Group=="DS UT"])))<5
aa
bb<-rowSums(is.na(exprs(esetN[,esetN$Group=="DS AOAA"])))<5
cc<-rowSums(is.na(exprs(esetN[,esetN$Group=="Ctrl UT"])))<5
dd<-rowSums(is.na(exprs(esetN[,esetN$Group=="Ctrl AOAA"])))<5

length(aa)
sum(aa)

esetN_filtered<-esetN[aa & bb & cc & dd ,]
esetN_filtered

#View(exprs(esetN_filtered))


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



# Create new ExpressionSet to store processed data
esetN_norm <- esetN_filtered
dim(esetN_norm)

#View(exprs(esetN_norm))


# #install.packages("imputeLCMD")
# library(imputeLCMD)
# # perform missing data imputation
# sum(is.na(exprs(esetN_norm)))
# 
# obj.QRILC <- impute.QRILC(exprs(esetN_norm),tune.sigma = 1)
# sum(obj.QRILC[[1]]<0)
# 
#  obj.QRILC[[1]][obj.QRILC[[1]]<0]<-0
# exprs(esetN_norm) <- obj.QRILC[[1]]
# 
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
# 
# 
# dim(exprs(esetN_norm))
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

# # Log tranform
# exprs(esetN_norm) <- log2(exprs(esetN_norm))
# plotDensities(esetN_norm, legend = FALSE)

boxplot(exprs(esetN_norm))
# Quantile normalize
exprs(esetN_norm) <- normalizeBetweenArrays(exprs(esetN_norm))
plotDensities(esetN_norm, legend = FALSE)
boxplot(exprs(esetN_norm))





names(pData(esetN_norm))
pData(esetN_norm)$Group <-factor(pData(esetN_norm)$Group)
pData(esetN_norm)$Sex <-factor(pData(esetN_norm)$Sex)
pData(esetN_norm)$Age <-as.numeric(as.character(pData(esetN_norm)$Age ))


pData(esetN_norm)$Age

levels(pData(esetN_norm)$Group)
pData(esetN_norm)$Group <- factor(pData(esetN_norm)$Group, levels = c( "Ctrl UT"  , "Ctrl AOAA","DS UT" ,    "DS AOAA" ))  


#Study design
design <- with(pData(esetN_norm),model.matrix(~0 + Group + Age + Sex))

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
df_statsum<-summary(results)
Total<-apply(df_statsum,2,sum)
df_statsum2<-rbind(df_statsum,Total)
colnames(df_statsum2)<-c("DS vs. Ctrl UnTreated","DS vs. Ctrl Treated",
                      "AOAA vs. UnTreated Ctrl","AOAA vs. UnTreated DS",
                      "intercept")
df_statsum2 
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

results <- decideTests(fit)
summary(results)


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


head(DS_to_Ctrl.UT)
library(writexl)
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
           paste0(path,"Flux_tables.xlsx"))





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
  return(tmp)
}

pca_plot(x=esetN_norm, group='Group', names="Cell_line", components= c(3, 2),main="PCA")

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

pca_plot(x=esetN_norm, group='Sex', names="Cell_line", main="PCA")

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

pca_plot(x=esetN_norm, group='Age', names="Cell_line", main = "PCA")

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

pca_plot(x=esetN_norm, group="Type" , names="Cell_line", main="PCA")

ggsave(
  filename="PCA_Type.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 4,
  height = 3,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

pca_plot(x=esetN_norm, group="AOAA_Treatment" , names="Cell_line", main="PCA")


ggsave(
  filename="PCA_Treatment.png",
  plot = last_plot(),
  path = path,
  scale = 1,
  width = 8,
  height = 6,
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
  
  ggp
  
}



################################################################
make_barplot<-function(x, 
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
  r <- runif(nrow(df), -0.15, 0.15)
  
  pt1 = t.test(Values ~ Group,
               data=df[df$Group %in% c("Ctrl AOAA", "Ctrl UT"), ],paired=T)
  print(pt1$p.value)
  pt2 = t.test(Values ~ Group,
               data=df[df$Group %in% c("Ctrl UT", "DS UT"), ],paired=F)
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

make_grupedbarplot<-function(x, 
                       isotopologues="",
                       group) {

  require(ggplot2)
  require (tidyr)
  stopifnot(is(x, 'ExpressionSet'))
  if (!missing(isotopologues)) {
    stopifnot(isotopologues %in% rownames(exprs(x)))
  }
  if (!missing(group)) {
    stopifnot(group %in% colnames(pData(x)))
  }

  
 
  df <-data.frame(exprs(x)[rownames(exprs(x)) %in% isotopologues,]
                 # ,Group=pData(x)[[group]]
                  )
  colnames(df) <-pData(x)[[group]]
  df$isotop=rownames(df)
  titlestring=gsub(" m\\+.*","",df$isotop[1])
  df$isotop<-gsub(".*m\\+","m\\+",df$isotop)
  df2<-pivot_longer(df,!isotop, names_to = "Groups", values_to = "values")
  df2$Groups <- factor(df2$Groups,
                       levels = c("Ctrl UT", "Ctrl AOAA" , "DS UT","DS AOAA"))
  
  
  ggp<- ggplot(df2)+
    aes(x=isotop, y=values, fill=Groups, group=Groups) +
    stat_summary(fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 fun.min = function(x) mean(x) , 
                 geom = 'errorbar',position = position_dodge(0.8),width=0.3, colour="grey30", size=1)+
    stat_summary(fun= 'mean', geom = 'bar',position = "dodge", alpha=0.3,size=1,
                 width=0.8, colour="grey30")+
    coord_cartesian(ylim=c(min(df2$values)-0.01*min(df2$values),
                           max(df2$values)+0.02*max(df2$values)))+
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
    theme(legend.position = "bottom")+
    theme(axis.title.x=element_blank())+
    labs(title=titlestring,
         y="Log2 normalised intensity")
  
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


compounds<-rownames(exprs(esetN_norm))
compounds

make_grupedbarplot(
  esetN_norm,
  isotopologues=c("Glucose-6-phosphate m+0",           
                  "Glucose-6-phosphate m+1",           
                  "Glucose-6-phosphate m+2",           
                  "Glucose-6-phosphate m+3",           
                  "Glucose-6-phosphate m+4",           
                   "Glucose-6-phosphate m+5",           
                  "Glucose-6-phosphate m+6" ),
  group="Group")
dir.create(paste0(path,"Groupbar"))
unlink(paste0(paste0(path,"Groupbar\\"),
              list.files(paste0(path,"Groupbar\\"))))

ggsave("Glucose-6-phosphate.jpg", 
       path=paste0(path,"Groupbar\\"),
       scale = 1,
       width = 5,
       height = 3,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

 
make_grupedbarplot(
  esetN_norm,
  isotopologues=c("Glutathione m+0" , "glutathione m+1" ,               
                   "glutathione m+2" ,    "glutathione m+3"                
                  , "glutathione m+4", "glutathione m+5"),
  group="Group")

make_grupedbarplot(
  esetN_norm,
  isotopologues=c("DHAP m+0",                          
                  "DHAP m+1",                          
                   "DHAP m+3"  ),
  group="Group")

make_grupedbarplot(
  esetN_norm,
  isotopologues=c("3-Phosphoglycerate m+0" ,           
                  "3PG m+1" ,                          
                   "3PG m+2" ,                          
                   "3PG m+3"   ),
  group="Group")

make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "2PG m+0" ,                          
                  "2PG m+1" ,                          
                  "2PG m+2" ,                          
                  "2PG m+3"   ),
  group="Group")

make_grupedbarplot(
  esetN_norm,
  isotopologues=c("Phosphoenolpyruvate m+0" ,"PEP m+2",                          
                   "PEP m+3"),
  group="Group")
                        
make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "Pyruvate m+0" ,                     
                  "Pyruvate m+1" ,                     
                  "Pyruvate m+2" ,                     
                 "Pyruvate m+3"  ),
  group="Group")

ggsave("Pyruvate.jpg", 
       path=paste0(path,"Groupbar\\"),
       scale = 1,
       width = 5,
       height = 3,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
 
make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "Lactate m+0" ,                      
                   "Lactate m+1" ,                      
                   "Lactate m+2" ,                      
                   "Lactate m+3"  ),
  group="Group") 

compounds<-rownames(exprs(esetN_norm))
compounds

make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "citrate/isocitrate m+0" ,           
                   "citrate/isocitrate m+1",            
                   "citrate/isocitrate m+2" ,           
                   "citrate/isocitrate m+3" ,           
                   "citrate/isocitrate m+4",            
                   "citrate/isocitrate m+5" ,           
                   "citrate/isocitrate m+6"   ),
  group="Group")


ggsave("Citrate.jpg", 
       path=paste0(path,"Groupbar\\"),
       scale = 1,
       width = 5,
       height = 3,
       units = "in",
       dpi = 300,
       limitsize = TRUE)



make_grupedbarplot(
  esetN_norm,
  isotopologues=c(  "Cis-aconitate m+0" ,                
                   "cis-aconitate m+1" ,                
                    "cis-aconitate m+2" ,                
                   "cis-aconitate m+3" ,                
                    "cis-aconitate m+4",                 
                   "cis-aconitate m+5" ,                
                   "cis-aconitate m+6"   ),
  group="Group")   

make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "aKG (Alpha-Ketoglutaric acid) m+0" ,
                    "aKG m+2",                           
                     "aKG m+3"  ,                         
                    "aKG m+4"  ,                         
                     "aKG m+5"   ),
  group="Group")  

make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "Succinic acid m+0" ,                
                   "Succinic acid m+1" ,                
                   "Succinic acid m+2" ,                
                   "Succinic acid m+3"     ),
  group="Group")               
                          

make_grupedbarplot(
  esetN_norm,
  isotopologues=c("Malate m+0","Malate m+1","Malate m+2",
                  "Malate m+3", "Malate m+4" ),
  group="Group")  

compounds<-rownames(exprs(esetN_norm))
compounds

make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "5Ribose 5P/Ribulose-5P/Xu-5P m+0",  
                  "5Ribose 5P/Ribulose-5P/Xu-5P m+1" , 
                  "5Ribose 5P/Ribulose-5P/Xu-5P m+2",  
                  "5Ribose 5P/Ribulose-5P/Xu-5P m+3",  
                  "5Ribose 5P/Ribulose-5P/Xu-5P m+4",  
                  "5Ribose 5P/Ribulose-5P/Xu-5P m+5"),
  group="Group") 

make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "Erythrose 4-phosphate m+0" ,        
                   "erythrose 4-phosphate m+1" ,        
                   "erythrose 4-phosphate m+2" ,        
                   "erythrose 4-phosphate m+3" ),
  group="Group") 

make_grupedbarplot(
  esetN_norm,
  isotopologues=c( "Sedoheptulose-7 phosphate"  ,       
                   "Sedoheptulose-7 phosphate m+2"  ,   
                   "Sedoheptulose-7 phosphate m+4" ,    
                   "Sedoheptulose-7 phosphate m+5",     
                   "Sedoheptulose-7 phosphate m+6",     
                   "Sedoheptulose-7 phosphate m+7"   ),
  group="Group") 
 
       
   

colnames(pData(esetN_norm))
make_barplot(esetN_norm,"Cystathionine",
             group = "Group", pair = "Cell_line")

make_barplot(esetN_norm,"ADP +4" ,
             group = "Group", pair = "Cell_line")

dir.create(paste0(path,"Figures\\"))

compounds<-rownames(exprs(esetN_norm))
compounds
#delete content
unlink(paste0(paste0(path,"Figures\\"),
              list.files(paste0(path,"Figures\\"))))
for( i in 1:nrow(exprs(esetN_norm))){
  make_barplot(esetN_norm,selected_metabolite = compounds[i],
               group = "Group", pair = "Cell_line")
  ggsave(paste0(make.names(compounds[i]),".jpg"), 
         path=paste0(path,"Figures\\"),
         scale = 1,
         width = 3.1,
         height = 1.9,
         units = "in",
         dpi = 300,
         limitsize = TRUE) 
}

