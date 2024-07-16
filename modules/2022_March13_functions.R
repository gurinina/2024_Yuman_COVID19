mycollapse2 = function(sresp,limit = 20){
lens = sapply(sresp,length)
wlens = which(lens > limit)
wresp= lapply(sresp[wlens], function(x) x = x[1:limit])
lsresp = lapply(sresp,function(x) paste(x,collapse = "|"))
lwresp = lapply(wresp,function(x) paste(x,collapse = "|"))
lsresp[wlens]=lwresp
lsresp
}
  ##################################################
stringWINDOW = function(x, width = 25){

  strng = paste(strwrap(x,width = width), collapse="\n")

  strng

}

barplotHEIGHT = function(df){

  scoreRange = range(df$score)

  dataRange <- scoreRange[2] - scoreRange[1]

  barLens <- round((df[, "score"] - scoreRange[1])/dataRange *100)
  w <- 150
  barWidth = 15
  h <- barWidth * length(unlist(barLens)) + 50
  h
}


genebarHEIGHT = function(leadInfo) {
  scoreRange = range(leadInfo$score)

  dataRange <- scoreRange[2] - scoreRange[1]

  barLens <- round((leadInfo[, 2] - scoreRange[1]) / dataRange *
                     100)
  w <- 150
  barWidth = 15
  h <- barWidth * length(barLens) + 50
  h
}

###############

compSCORE <- function(mat,coln, sig = 1){

  df = data.frame(score = mat[,coln],stringsAsFactors = F)
  library(dplyr)
  df$gene = rownames(mat)
  rownames(df) = df$gene
  df$index=0
  wdf = which(df$score >= sig)
  df$index[wdf]=1
  df = df[,c('index','score','gene')]
  df = df %>% arrange(desc(score))
  df
}


compRANK <- function(mat,coln, rank = 30){

  df = data.frame(score = mat[,coln],stringsAsFactors = F)
  library(dplyr)
  df$gene = rownames(mat)
  rownames(df) = df$gene
  df$index=0
  wdf = which(df$score <= rank)
  df$index[wdf]=1
  df = df[,c('index','score','gene')]
  df = df %>% arrange(score)
  df
}



############ GENE ANNOTATION
# generates annotation for the HOP fitness plots for a single screen from a matrix
# mat - data matrix; each colname is a screen. values are defect scores, rownames are gene names
# fdat - gene annotation data.frame file with genes, ORFS, descriptors, etc.
# cmp = name of screen of interest; column name
# sgdlink - provides link to SGD for gene names if desired
# enrichtype = FD for FITNESS DEFECT; FA for FITNESS ADVANTAGE (resistant)
# out put is shown on the HOP fitness tap underneath the plots; used as input for the ggfit plots below
####
####

myhyperlinkblank = function(x,href,y=x){
link = paste0("<a href=",href)
link2 = paste0(link,x)
t=" target='_blank'"
link3 = paste(link2,noquote(t))
link4 = paste0(link3,">",y)
link5= paste0(link4,"</a>")
link5
}
####

# ############
geneAnno <- function(mat,fdat = NULL,cmp,xvar = T,
                    enrichtype = "fd",arrange = F){

  mystrsplit = function(str,splt,index){
  xsplt = sapply(strsplit(str,splt), function(x) x = x[index])
  xsplt
  }

  if(is.null(fdat)) {
  fdat<-  readRDS("24Jan8/orth.RDS")
          }

  mat = mat [order(rownames(mat)),,drop = FALSE]
  w1 = which(colnames(mat) %in% cmp)
  req(length(w1) > 0)
   validate(
      need(length(w1) == 1 ,message =
          "please enter a valid compound"))

  dmat = data.frame(gene = rownames(mat),FD = mat[,w1],stringsAsFactors = F,check.names = F)
  dmat = dmat %>% arrange(desc(FD))

  m = match(dmat$gene,fdat$best_human_gene)

  dmat[,c("GENE","descriptor")] = NA
  if(xvar) {

    dmat[,c("gene","descriptor","GENE")] =  fdat[m,c("best_human_gene","descriptor","genecards_link")]
       dmat$xvar = as.numeric(1:nrow(dmat))
   dmat = dmat[,c("xvar","FD","gene","GENE","descriptor")]
      } else if(xvar==F){
        dmat[,c("gene","descriptor","GENE")] = fdat[m,c("best_human_gene","descriptor","genecards_link")]
        dmat = dmat[,c("FD","gene","GENE","descriptor")]}

  dmat

}
##########
mystrsplit = function(str,splt,index){
  xsplt = sapply(strsplit(str,splt), function(x) x = x[index])
  xsplt
}

####
visSetup = function(enrichInfo, edgeMat, fontsize = 22, fontface = "Arial") {

  library(igraph)
  library(visNetwork)
  n = enrichInfo
  e = edgeMat

  w = which(names(n) == "id")
  coln = (1:ncol(n))[-w]
  n = n[, c(w, coln)]

  if (is.null(e) & !is.null(n)) {
    gr = make_empty_graph(nrow(enrichInfo))
    v = gr
    V(v)$color.background = n$cluster
    v = set_vertex_attr(v, "label", value = n$formattedLabel)

  }


  if (!is.null(e) & !is.null(n)) {
    w = which(names(e) == "label")
    let = graph_from_data_frame(e[, -w], vertices = n, directed = F)
    v = set_vertex_attr(let, "label", value = n$formattedLabel)
    V(v)$color.background = n$cluster
  }

  vis = toVisNetworkData(v)
  vis$nodes = data.frame(vis$nodes, stringsAsFactors = F)
  if (!is.null(vis$edges))
    vis$edges = data.frame(vis$edges, stringsAsFactors = F)
  m = match(vis$nodes$id, n$id)
  vis$nodes$label = n$formattedLabel[m]
  vis$nodes$FDR = n$FDR[m]
  vis$nodes$FDR = signif(vis$nodes$FDR, 5)


  w = which(duplicated(vis$nodes$FDR))
  if (length(w) > 0) {
    vis$nodes =  vis$nodes %>% group_by(FDR) %>% mutate(
      jitter = if (n() > 1)
        abs(jitter(FDR))
      else
        (FDR)
    )

    w = which(names(vis$nodes) == "FDR")
    vis$nodes = vis$nodes[, -w]
    w = which(names(vis$nodes) == "jitter")
    names(vis$nodes)[w] = "FDR"
    vis$nodes$FDR = signif(vis$nodes$FDR, 6)
  }

  w = which(duplicated(vis$nodes$FDR))



  vis$nodes$term = n$term[m]
  vis$nodes$interSect = n$interSect[m]
  vis$nodes$nQuery= n$nQuery[m]
  vis$nodes$nGenes = n$nGenes[m]
  vis$nodes$geneSetFraction = n$geneSetFraction[m]
  vis$nodes$querySetFraction = n$querySetFraction[m]
  vis$nodes$filename =  n$filename[m]

  vis$nodes$formattedLabel = n$formattedLabel[m]
  vis$nodes$overlapGenes = n$overlapGenes[m]
  vis$nodes$label = vis$nodes$formattedLabel
  vis$nodes$color.border = "black"
  vis$nodes = vis$nodes %>% arrange(label)
  if (nrow(vis$edges) > 0)
    vis$edges$color = "black"

  vis$nodes$borderWidthSelected	= 4
  w = which(names(vis$nodes) %in% c("fomattedLabel", "color", "cluster"))
  if (length(w) > 0)  vis$nodes = vis$nodes[, -w]
  vis$nodes$color.highlight.border = "#000066"
  vis$nodes$color.highlight.background = "#c0b3ff"
  vis$nodes$color.hover.background = "#000066"
  vis$nodes$color.hover.border = "#c0b3ff"
  vis$nodes$font.face = fontface
  vis$nodes$shape = "dot"
  vis$nodes$font.size = fontsize
  vis$nodes$font.bold = F
  vis$nodes$borderWidth = 2
  #vis$nodes$vadjust = "mono"
  vis$nodes$borderWidthSelected = 4
  vis$nodes$labelHighlightBold = T
  w = which.min(vis$nodes$FDR)
  if (length(w) > 0) {
    vis$nodes$font.size[w] = fontsize + 4
    vis$nodes$borderWidth[w] = 4
    vis$nodes$font.bold[w] = T
  }
  vis$nodes = vis$nodes %>% arrange(FDR)
  vis

}


setnicepar = function(mar = c(3, 3, 2, 1), mgp = c(2, 0.4, 0),
  tck = -.01, cex.axis = 0.9,
  las = 1, mfrow = c(1, 1), ...) {
  par(mar = mar, mgp = mgp, tck = tck,
    cex.axis = cex.axis, las = las,
    mfrow = mfrow, ...)
}



######### PLOTTING FUNCTION
############
# generate leading edge barplots for all clusters
# scoreMat - dataframe: ID, score, optional 3rd column indicating which ones to mark (TRUE/FALSE)
# geneSets - a list of all gene sets, where each list element is a vector of gene IDs
# gsInfo - dataframe of gene set info: id=gene set id, cluster=cluster id, es=enrichment score, fdr=FDR value
# scoreName - score label to use in leading edge plots
# plotCol - the colour of the bars, in hexidecimal format without the '#' character
# RETURNS a dataframe of barplot node info, each row corresponds to a different node;
# contains "id", "image" (google chart URL), "w" (plot width), "h" (plot height), "cluster" columns
#genLeadingEdgePlot.all <- function(scoreMat, geneSets, gsInfo, scoreName, plotCol="BEBEBE") {
####

geneBARPLOT = function(overlapGenes, desc = FALSE){
library(dplyr)
strSPLIT <- function(str,splt,index){
xsplt = sapply(strsplit(str,splt), function(x) x = x[index])
xsplt
}
splt = strsplit(overlapGenes, "\\|")
gene = lapply(splt,strSPLIT,"\\(",1)
score = lapply(splt,strSPLIT,"\\(",2)
score2 = lapply(score,strSPLIT,"\\)",1)
score3 = lapply(score2,as.numeric)
leadInfo = data.frame(gene = unlist(gene),score = unlist(score3),stringsAsFactors = F)
if(nrow(leadInfo) > 10) leadInfo = leadInfo[1:10,]

if(desc){ leadInfo = leadInfo %>% arrange(desc(score))}else {

  leadInfo = leadInfo %>% arrange(score)

}

leadInfo

}


genebarHEIGHT = function(leadInfo) {
  scoreRange = range(leadInfo$score)

  dataRange <- scoreRange[2] - scoreRange[1]

  barLens <- round((leadInfo[, 2] - scoreRange[1]) / dataRange *
                     100)
  w <- 150
  barWidth = 15
  h <- barWidth * length(barLens) + 50
  h
}



############
mymeltdf = function(mat, row, df = phop){
  library(reshape2)

  library(dplyr)
  mx = melt(mat[row,], value.name = "fitness_defect")
  mx$gene = row
  mx$screen = rownames(mx)

  mx = mx[,c("gene","screen","fitness_defect")]
  m = match(mx$screen,df$name)
  table(is.na(m))
  #mx$id = df$id[m]
  mx$PCID = df$PCID[m]
  mx$cmp = df$cmp[m]

  mx$site = df$site[m]

  mx

}

##############

p10 <- function(matI,x,sig = 1,ylab1 = "fitness defect",xlab1 = "gene",las1 = 2,font1 = 3.3,cex1 = 1,pch1 = 17,
                cex.main1 = 1.5,cex.axis1 = 1.5,cex.lab1=1.5,cex.text1 = 1,...)  {

  mycolors =
    c(
      "darkorange",
      "dodgerblue",
      "limegreen"
      )

  col = ifelse (matI[, x] > sig, 1,      ifelse (matI[, x] < -sig, 3, 2))
  w <- which(matI[, x]  > sig | matI[, x]  < -sig, arr.ind = T)
  palette(mycolors)
  posw = ifelse (w > nrow(matI) - 0.1 * nrow(matI) ,2, ifelse (w <  nrow(matI) + 0.1 * nrow(matI), 4 ,  2))
  plot(
    matI[, x],
    col = col ,cex.lab = cex.lab1,
    main = colnames(matI)[x],
    ylab = ylab1, xaxt = "n",
    xlab = xlab1,
    las = las1,font = 2, font.lab = 2, font.main = 2,
    cex = cex1,
    cex.axis = cex.axis1,
    cex.main = cex.main1,
    pch = pch1,
    ...
  )
  if (length(w != 0))
    text(
      w,
      matI[w, x],
      names(w),
      pos = posw,
      cex = cex.text1,
      font = font1,
      ...
    )
  abline(
    h = sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
  abline(
    h = -sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
}

###########
mycolors = c(
  "darkorange1",
  "dodgerblue",
  "darkgreen",
  "navy",
  "mediumpurple"  ,
  "royalblue3",
  "darkolivegreen4",
  "firebrick",
  "cyan4",
  "hotpink3",
  "plum4",
  "blue",
  "magenta4",
  "skyblue3",
  "green4",
  "red3",
  "steelblue3",
  "tomato",
  "purple4",
  "goldenrod3",
  "steelblue",
  "darkred",
  "lightpink3",
  "darkorchid",
  "lightblue3",
  "dimgrey",
  "chocolate1",
  "seagreen3",
  "darkkhaki",
  "darksalmon")
