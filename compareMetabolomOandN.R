

pathO="C:\\Users\\Acer 3\\Documents\\CsabaII\\MetabolomOLD\\"
DStoCTR_Old <- readRDS(paste0(pathO,"DStoCTR_Old.rds"))


pathN="C:\\Users\\Acer 3\\Documents\\CsabaII\\MetabolomNEW\\"


DStoCTR_New <- readRDS(paste0(pathN,"DStoCTR_New.rds"))

head(DStoCTR_New)
# Load library
install.packages("VennDiagram")
library(VennDiagram)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")

# Generate 3 sets of 200 words
setO <- rownames(DStoCTR_Old)
setN <- rownames(DStoCTR_New) 
venn.diagram(
  x = list(setO, setN),
  category.names = c("Metab Old" , "Metab New "),
  filename =paste0(pathN, '#14_venn_diagramm.png'),
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("lightblue","yellow"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.6,
 cat.fontface = "bold",
cat.default.pos = "outer",
  cat.pos = c(-27, 27),
 cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans")

head(DStoCTR_Old)
upregOLD <- DStoCTR_Old[DStoCTR_Old$P.Value<0.1 & DStoCTR_Old$logFC>0,]
upregOLD
upregNEW <- DStoCTR_New[DStoCTR_New$P.Value<0.1 & DStoCTR_New$logFC>0,]
upregNEW

# Generate 3 sets of 200 words
setO <- rownames(upregOLD)
setN <- rownames(upregNEW) 
venn.diagram(
  x = list(setO, setN),
  category.names = c("Metab Old" , "Metab New "),
  filename =paste0(pathN, '#upregulated.png'),
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("lightblue","yellow"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-12, 12),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans")

#common elements
intersect(setO,setN)

head(DStoCTR_Old)
downregOLD <- DStoCTR_Old[DStoCTR_Old$P.Value<0.1 & DStoCTR_Old$logFC<0,]
downregOLD
downregNEW <- DStoCTR_New[DStoCTR_New$P.Value<0.1 & DStoCTR_New$logFC<0,]
downregNEW


# Generate 3 sets of 200 words
setO <- rownames(downregOLD)
setN <- rownames(downregNEW) 
venn.diagram(
  x = list(setO, setN),
  category.names = c("Metab Old" , "Metab New "),
  filename =paste0(pathN, '#downregulated.png'),
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("lightblue","yellow"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-12, 12),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans")

#common elements
intersect(setO,setN)

#is this random
SA<-sample(rownames(DStoCTR_New), 50, replace = FALSE, prob = NULL)
SB<-sample(rownames(DStoCTR_Old), 20, replace = FALSE, prob = NULL)

#common elements
intersect(SA,SB)


