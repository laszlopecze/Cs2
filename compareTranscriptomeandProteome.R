library(VennDiagram)
Sys.setenv(LANG = "en")


pathT="C:\\Users\\Acer 3\\Documents\\CsabaII\\mRNA\\"
DStoCTR_UT_mRNA <- readRDS(paste0(pathT,"DStoCTR_UT_mRNA.rds"))
DStoCTR_T_mRNA <- readRDS(paste0(pathT,"DStoCTR_T_mRNA.rds"))
AOAAtoUT_CTR_mRNA <- readRDS(paste0(pathT,"AOAAtoUT_CTR_mRNA.rds"))
AOAAtoUT_DS_mRNA <- readRDS(paste0(pathT,"AOAAtoUT_DS_mRNA.rds"))


pathP="C:\\Users\\Acer 3\\Documents\\CsabaII\\Protein\\"


DStoCTR_UT_Prot <- readRDS(paste0(pathP,"DStoCTR_UT_Prot.rds"))
DStoCTR_T_Prot <- readRDS(paste0(pathP,"DStoCTR_T_Prot.rds"))
AOAAtoUT_CTR_Prot <- readRDS(paste0(pathP,"AOAAtoUT_CTR_Prot.rds"))
AOAAtoUT_DS_Prot <- readRDS(paste0(pathP,"AOAAtoUT_DS_Prot.rds"))

#path for results
pathC="C:\\Users\\Acer 3\\Documents\\CsabaII\\Comparison_mRNA_Prot_0p05\\"
dir.create(pathC)

#User defined function#########################################################

makeVenn<-function(df1, df2, 
                   regulation=c("up","down","total"),
                   signl=0.05)
{
  
if (regulation=="up"){
  df_mRNA <- df1[df1$P.Value<signl & df1$logFC>0,]
  df_Prot <- df2[df2$P.Value<signl & df2$logFC>0,]
  dfname<-deparse(substitute(df1))
  dfname2<-gsub("mRNA", "upreg",dfname)
}else if(regulation=="down"){
  df_mRNA <- df1[df1$P.Value<signl & df1$logFC<0,]
  df_Prot <- df2[df2$P.Value<signl & df2$logFC<0,]
  dfname<-deparse(substitute(df1))
  dfname2<-gsub("mRNA", "downreg",dfname)
}else{
  df_mRNA <- df1
  df_Prot <- df2
  dfname<-deparse(substitute(df1))
  dfname2<-gsub("mRNA", "all",dfname)
}
  setT <- rownames(df_mRNA)
  setP <- rownames(df_Prot) 
  venn.diagram(
    x = list(setT, setP),
    category.names = c("Transcriptome" , "Proteome"),
    filename =paste0(pathC, '#',dfname2,'.png'),
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
  ce<-intersect(setT,setP)
  tt<-DStoCTR_UT_mRNA[rownames(DStoCTR_UT_mRNA) %in% ce,3:4]
  pp<-DStoCTR_UT_Prot[rownames(DStoCTR_UT_Prot) %in% ce,1:2]
  res<-merge(tt,pp,by="row.names",all.x=TRUE)[-1]
  return(res)

}

du1<-makeVenn(DStoCTR_UT_mRNA, DStoCTR_UT_Prot, regulation="up")
dd1<-makeVenn(DStoCTR_UT_mRNA, DStoCTR_UT_Prot, regulation="down")
makeVenn(DStoCTR_UT_mRNA, DStoCTR_UT_Prot, regulation="total")


du2<-makeVenn(DStoCTR_T_mRNA, DStoCTR_T_Prot, regulation="up")
dd2<-makeVenn(DStoCTR_T_mRNA, DStoCTR_T_Prot, regulation="down")
makeVenn(DStoCTR_T_mRNA, DStoCTR_T_Prot, regulation="total")

du3<-makeVenn(AOAAtoUT_CTR_mRNA, AOAAtoUT_CTR_Prot, regulation="up")
dd3<-makeVenn(AOAAtoUT_CTR_mRNA, AOAAtoUT_CTR_Prot, regulation="down")
makeVenn(AOAAtoUT_CTR_mRNA, AOAAtoUT_CTR_Prot, regulation="total")

du4<-makeVenn(AOAAtoUT_DS_mRNA, AOAAtoUT_DS_Prot, regulation="up")
dd4<-makeVenn(AOAAtoUT_DS_mRNA, AOAAtoUT_DS_Prot, regulation="down")
makeVenn(AOAAtoUT_DS_mRNA, AOAAtoUT_DS_Prot, regulation="total")


library(writexl)
write_xlsx(list(`DS vs. Ctrl UnTreated`=du1,
                `DS vs. Ctrl Treated`=du2,
                `AOAA vs. UnTreated Ctrl`= du3,
                `AOAA vs. UnTreated DS`=du4),
           paste0(pathC,"BothUpregulated_tables.xlsx"))

write_xlsx(list(`DS vs. Ctrl UnTreated`=dd1,
                `DS vs. Ctrl Treated`=dd2,
                `AOAA vs. UnTreated Ctrl`= dd3,
                `AOAA vs. UnTreated DS`=dd4),
           paste0(pathC,"BothDownregulated_tables.xlsx"))

