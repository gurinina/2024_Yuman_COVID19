
# Helper function to generate img src tag
generate_img_tag <- function(path, border = "3px solid black",
                             width = "200px", height = "200px") {

  img_tag <- sprintf(
    '<img src="%s" style="border: %s; width: %s; height: %s;">',
    path, border, width, height
  )
  return(img_tag)
}

slider = function(x=input$corrFIT, y= input$scoreFIT, mfit){

  score =  as.numeric(y)

  w = which(mfit[,2] == score)
  corr = mfit[w,1:3]

  cos = as.numeric(x)

  w = which(mfit[,3] == cos)

  cort = mfit[w,1:3]

  return(list(corr=cort,score=corr))
}

 slid = function(x=input$corrFIT, y= input$scoreFIT, mfit){

   score =  as.numeric(y)

   w = which(mfit[,2] == score)
   corr = mfit[w,1:3]

   cos = as.numeric(x)


   w = which(mfit[,3] == cos)


   cort = mfit[w,1:3]

  return(list(corr=cort,score=corr))
}

mymeltdf = function(mat, row, df = phop){

  library(reshape2)

  library(dplyr)
  mx = melt(mat[row,], value.name = "fitness_defect")
  mx$gene = row
  mx$screen = rownames(mx)

  mx = mx[,c("gene","screen","fitness_defect")]
  m = match(mx$screen,df$name)

  mx$PCID = df$PCID[m]
  mx$cmp = df$cmp[m]

  mx$site = df$site[m]

  mx

}

stringWINDOW = function(x, width = 25){

  strng = paste(strwrap(x,width = width), collapse="\n")

  strng

}

############ GRAPHING FUNCTION
meltDF = function(mat, row, df = phiphop) {
  mx = reshape2::melt(mat[row, ], value.name = "fitness_defect")
  mx = data.frame(
    screen = rownames(mx),
    gene = row,
    fitness_defect = mx$fitness_defect,
    stringsAsFactors = F
  )

  mx$gene = paste0("<a href=https://www.yeastgenome.org/locus/",mx$gene,">",mx$gene,"</a>")

  mx = mx[, c("gene", "screen", "fitness_defect")]
  m = match(mx$screen, df$name)
  table(is.na(m))



  mx$target = df$target_html[m]

  mx$signature = df$signature[m]
  mx$shape = 19


  g = grep(row,df$target)
  if(length(g) > 0) mx$shape[g] = 17

  mx$PCID = df$PCID[m]
  mx$FDA = factor(df$FDA[m])
  mx$mechanism = df$mechanism[m]
  mx$target = df$target[m]
  n =   c("screen",
    "gene",
    "fitness_defect",
    "PCID",
    "mechanism",
    "target",
    "signature"   ,
     "FDA",
    "shape"
  )
  g = grep("fitness",colnames(mx))

  mx = mx[,n]
  mx
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


# hop = geneAnno(mat = x,colnames(x)[w1],fdat = fdat, sgdlink = T)
#
# pclick <- nearPoints(df = hop,  xvar = "gene", yvar = "FD",
# coordinfo = input$plot_click,threshold = 10)
# outhop = geneAnno(mat = xne,fdat = fdat,sgdlink = F,cmp = cmp)
# hotabp = geneAnno(mat = xne, fdat = fdat, arrange = T,
# sgdlink = T, cmp = input$cmpSERV)

myhyperlinkblank = function(x,href,y = x){
link = paste0("<a href=",href)
link2 = paste0(link,x)
t=" target='_blank'"
link3 = paste(link2,noquote(t))
link4 = paste0(link3,">",y)
link5= paste0(link4,"</a>")
link5
}

hgeneAnno <- function(mat, fdat = NULL, cmp, arrange = TRUE){

if(is.null(fdat)) {
    fdat <- readRDS("data/orthologs.RDS")
}

mat = mat [order(rownames(mat)),,drop = FALSE]

w1 = which(colnames(mat) %in% cmp)

req(length(w1) > 0)

validate(need(length(w1) == 1,
  message = "please enter a valid compound"))

dmat = data.frame(gene = rownames(mat),FD = mat[,w1],
  stringsAsFactors = F,check.names = F)

if(arrange) {dmat = dmat %>% arrange(desc(FD))}

m = match(dmat$gene,fdat$human_gene)

dmat[,c("GeneCards","descriptor")] = NA

dmat[,c("gene","GeneCards","descriptor")] =  fdat[m,c("human_gene","GeneCards","descriptor")]

dmat = dmat[,c("FD","gene","GeneCards","descriptor")]

dmat$FD = round(dmat$FD,3)

dmat

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

###############

compSCORE <- function(mat,coln, sig){

  df = data.frame(score = mat[,coln],stringsAsFactors = F)
  library(dplyr)
  df$gene = rownames(mat)
  rownames(df) = df$gene
  df$index=0
  wdf = which(df$score > sig)
  df$index[wdf]=1
  df = df[,c('index','score','gene')]
  df = df %>% arrange(desc(score))
  df
}

  ###### bulk of code


 #####################################################
 #####################################################
 #####################################################

############
# generates parameters for drawing GO enrichment networks
# enrichInfo - dataframe with enrichment stats for gene sets (one per row), with the following columns:
#             filename, term, geneSetFraction, FDR, overlapGenes, maxOverlapGeneScore
#           - see documentation for the output of hyperG() for descriptions of these columns
#           - rows with the same value, x, in the filename column specify enrichment results for
#             the set of query genes in queryGeneSets with name=x
# edgeMat - dataframe of edge set info: id=gene set id, cluster=cluster id, es=enrichment score, fdr=FDR value
#
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
