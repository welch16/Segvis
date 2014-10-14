## ----<style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------------------
  BiocStyle::latex()

## ----load,eval=FALSE--------------------------------------------------------------------
#    library(profile)

## ----real_load,include=FALSE,echo=FALSE,eval=TRUE---------------------------------------
  library(profile)
  rdatadir = "../inst/auxRData"

## ----create_prof,include=TRUE,echo=TRUE,eval=TRUE---------------------------------------
  name = "H3k27ac"
  file ="../inst/extdata/encode_K562_H3k27ac_first3chr.sort.bam"
  maxBw = 501
  fl = 200
  ourProfile = Profile(regionName = name,file = file,fileFormat= "bam",maxBandwidth = maxBw,fragLen = fl,remChr = "")
  ourProfile

## ----prof_mat,include=TRUE,echo=TRUE,eval=FALSE,tidy=TRUE-------------------------------
#    ourProfileMatrix1 = ProfileMatrix(ourProfile,1,mc)
#    ourProfileMatrix2 = ProfileMatrix(ourProfile,51,mc)
#    ourProfileMatrix3 = ProfileMatrix(ourProfile,301,mc)
#    ourProfileList = ProfileMatrixList(ourProfileMatrix1,ourProfileMatrix2,ourProfileMatrix3)
#    names(ourProfileList) = c("1","51","301")

## ----plot1,include=TRUE,echo=TRUE,eval=FALSE,tidy=TRUE----------------------------------
#    p1 = plot.profiles(ourProfileList)
#    p2 = plot.profiles(ourProfileList,coord=-windowExt:windowExt)

## ----plot1_load,include=FALSE,echo =FALSE,eval=TRUE-------------------------------------
  load(file = file.path(rdatadir,"bw_figs.RData"))

## ----plot1_both,include=TRUE,echo=FALSE,eval=TRUE,out.width='8 cm',fig.show = 'hold',fig.align='center'----
print(p1)
print(p2)

## ----plot2_load,include = FALSE,echo = FALSE,eval=TRUE----------------------------------
  load(file=file.path(rdatadir,"trim_figs.RData"))

## ----buildmatrix, include=TRUE,echo=TRUE,eval=FALSE,tidy = TRUE-------------------------
#    ourMatrices = lapply(ourProfiles,ProfileMatrix,251,mc)
#    names(ourMatrices) = names
#    ourList_notScaled = ProfileMatrixList(ourMatrices)

## ----normalize,include=TRUE,echo=TRUE,eval=FALSE,tidy=TRUE------------------------------
#    ourMatrices = lapply(ourMatrices,normalize.matrix)
#    ourList = ProfileMatrixList(ourMatrices)

## ----include=TRUE,echo=TRUE,eval=FALSE,tidy=TRUE----------------------------------------
#    q1 = plot.profiles(ourList_notScaled,coord= -windowExt:windowExt)
#    q2 = plot.profiles(ourList,coord= -windowExt:windowExt)

## ----plot2_both,include=TRUE,echo=FALSE,eval=TRUE,out.width='8 cm',fig.show = 'hold',fig.align='center'----
    print(q1)
    print(q2)

## ----include=TRUE,echo=TRUE,eval=FALSE--------------------------------------------------
#  q3 = plot.profiles(ourList,coord= -windowExt:windowExt,trim = 1)
#  q4 = plot.profiles(ourList,coord= -windowExt:windowExt,trim = .25)

## ----plot3_both,include=TRUE,echo=FALSE,eval=TRUE,out.width='8 cm',fig.show = 'hold',fig.align='center'----
  print(q3)
  print(q4)

## ----plot_format,include=TRUE,echo=TRUE,eval=TRUE,out.width='8 cm',fig.show ='hold',fig.align='center',tidy=TRUE----
  q2 +theme(legend.position = 'bottom')+scale_color_brewer(palette = 'Dark2')
  q2 +theme_bw()+theme(legend.position = 'bottom')+scale_color_brewer(palette = 'Accent')+
    geom_vline(xintercept=0,linetype = 'dashed')+scale_linetype_manual(values = rep('solid',3))+
    ylab('Nomalized signal')+xlab('Distance from summit')

## ----include=TRUE,echo=TRUE,eval=FALSE--------------------------------------------------
#    ourList = lapply(ourList,function(x,dnase_regions){
#      regions(x)$dnase = countOverlaps(regions(x),dnase_regions)
#      return(x)},dnase_regions)

## ----include=TRUE,echo=TRUE,eval = FALSE,tidy=TRUE--------------------------------------
#    p = plot.profiles(ourMatrices,coord = -windowExt:windowExt)
#    p1 = plot.profiles(ourMatrices,coord = -windowExt:windowExt,condition = dnase > 0 & pol2b > 0)
#    p2 = plot.profiles(ourMatrices,coord = -windowExt:windowExt,condition = dnase > 0 & pol2b == 0)
#    p3 = plot.profiles(ourMatrices,coord = -windowExt:windowExt,condition = dnase == 0 & pol2b > 0)

## ----include=FALSE,echo=FALSE,eval=TRUE-------------------------------------------------
  load(file = file.path(rdatadir,"subset_figs.RData"))

## ----plot_subset1,include=TRUE,echo=FALSE,eval=TRUE,out.width='8 cm',fig.show ='hold',fig.align='center',tidy=TRUE----
p +theme_bw()+theme(legend.position = 'bottom')+scale_color_brewer(palette = 'Dark2')+
  geom_vline(xintercept=0,linetype = 'dashed')+scale_linetype_manual(values = rep('solid',4))+
  ylab('Nomalized signal')+xlab('Distance from summit')+ggtitle('Without subsetting')
p1 +theme_bw()+theme(legend.position = 'bottom')+scale_color_brewer(palette = 'Dark2')+
  geom_vline(xintercept=0,linetype = 'dashed')+scale_linetype_manual(values = rep('solid',4))+
  ylab('Nomalized signal')+xlab('Distance from summit')+ggtitle('Overlaps with Dnase sensitive regions and Pol2b peaks')
p2 +theme_bw()+theme(legend.position = 'bottom')+scale_color_brewer(palette = 'Dark2')+
  geom_vline(xintercept=0,linetype = 'dashed')+scale_linetype_manual(values = rep('solid',4))+
  ylab('Nomalized signal')+xlab('Distance from summit')+ggtitle('Only overlaps with Dnase sensitive regions')
p3 +theme_bw()+theme(legend.position = 'bottom')+scale_color_brewer(palette = 'Dark2')+
  geom_vline(xintercept=0,linetype = 'dashed')+scale_linetype_manual(values = rep('solid',4))+
  ylab('Nomalized signal')+xlab('Distance from summit')+ggtitle('Only overlaps with Pol2b peaks')

## ----sessionInfo,include=TRUE,echo =TRUE,eval=TRUE,results="asis"-----------------------
  toLatex(sessionInfo())

