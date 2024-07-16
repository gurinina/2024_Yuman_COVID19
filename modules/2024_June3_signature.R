
fdat      <-     readRDS("24Jan8/2021_July19_fdat.RDS")

noness = fdat %>% filter(essential == "noness")
ess = fdat %>% filter(essential == "ess")

xsig      <-     readRDS("24Jan8/x21.RDS")

dfsig   <-  readRDS("24Jan8/dfsig.RDS")

choices = sort(dfsig$signature)

selected = "DNA damage"

sigNetworkModuleUI <- function(id, label = "sigNetworkModule") {

  ns <- NS(id)

tagList(

fluidPage(

fluidRow(

  column(width  = 9,

box(title = "Select response signature:",

  selectizeInput(ns("respMOD"), width = "500px",
label = NULL, selected = selected,
choices = choices, multiple = F),

status = "primary", solidHeader = T, width = "100%", height = 150)),

  column(width  = 3,

box(title = "Reset response menu:",

  br(),

  fluidRow(column(align = "center", width = 12,

 actionButton(ns("resetSig"), label = "reset respMenu"))),

status = "primary", solidHeader = T, width = "100%", height = 150))),

fluidRow(

  column(width = 5,

box(title = "Fitness threshold:",

  sliderInput(ns("scorethreshMOD"),label  = NULL, min = 0, max = 5.0,
value = 1.0, step = 0.1,ticks = TRUE),

status = "primary", solidHeader = T,width="100%", height = 170)),

  column(width = 4,

box(title = "FDR threshold:",

  sliderInput(ns("fdrMOD"), label = NULL, min = 0, max = 0.5,
value = 0.2, step = 0.1, ticks = TRUE),

status = "primary", solidHeader = T,width="100%", height = 170)),

  column(width = 3,

    box(title= "Background genome:", align = "center",

  prettyRadioButtons(ns("poolset"), label = NULL,
    choices = c("HIP" = "essential",
    "HOP" = "nonessential","
     HIPHOP" = "all"),

    status = "primary", shape = "square", animation = "tada",
    bigger = T,selected = "all", inline = T),

    solidHeader = T,status = "primary",width = "100%",height = 170))),

fluidRow(

  box(title = "GO enrichment analysis of the response signature:",

fluidPage(
fluidRow(downloadButton(ns("downloadSignatures"),
  "download response signatures")),

fluidRow(br(),
  visNetworkOutput(ns("network_Proxy"), width = "100%",height = 855))),

  width = 9, status = "primary", solidHeader = TRUE, height = 1000),

  box(title = "Enrichment details:",width = 3,

br(),

uiOutput(ns("Gotable")),

  solidHeader = T,status = "primary",background = "navy", height = 270),

  box(title = "Leading edge genes:",width = 3,

uiOutput(ns("LeadingEdge")),

  solidHeader = T,status = "primary", height = 345),

box(title = "Response signature genes:", width = 3,

   uiOutput(ns("sigOutput")),

   status = "primary", solidHeader = T, height = 345)),

fluidRow(

  box(title = "GO enrichment table:",

DT::dataTableOutput(ns("enrichSig")),

  status = "primary", solidHeader = TRUE,width = 12)))

)

}
###########
sigNetworkModule = function(input, output, session,
xRespDat, xRespInput, inputTab,
message = "No GO enrichment, try relaxing the FDR or scorethreshold"){

  datRECVD  =   reactive({xRespDat()})
  respRECVD = reactive({xRespInput()})
  tabRECVD  = reactive({inputTab()})

###############################################
############## START OBSERVEEVENTS ############
###############################################

observeEvent(respRECVD() ,{

req(!is.null(respRECVD()))

drugMOD = input$respMOD

sigRECVD =  respRECVD()

ns <- session$ns
req(sigRECVD)
req(length(sigRECVD) > 0)

if(sigRECVD != drugMOD){

  updateSelectizeInput(session, 'respMOD', choices = choices, selected = sigRECVD)

  updateSliderInput(session,'scorethreshMOD',
    label = NULL, min = 0, max = 5,
    value = 1, step = 0.5)

  updateSliderInput(session,'fdrMOD',label = NULL,
    min = 0, max = 0.5,
    value = 0.2, step = 0.05)

}

},ignoreInit = F, ignoreNULL = T)
########################

observeEvent(input$respMOD,{

ns <- session$ns

req(input$respMOD)

updateSliderInput(session,'scorethreshMOD',
  label = NULL, min = 0, max = 5.0,
  value = 1.0, step = 0.1)

updateSliderInput(session,'fdrMOD',label = NULL, min = 0,
  max = 0.5, value = 0.2, step = 0.05)

tabINP = tabRECVD()

returnedMOD$tab = tabINP

returnedMOD$sig = input$respMOD

},ignoreInit = F, ignoreNULL = T)

observeEvent(input$respMOD,{

req(input$respMOD)
tabMOD = tabRECVD()
returnedMOD$sig = input$respMOD

}, ignoreInit = F, ignoreNULL = T)

observeEvent(input$resetSig,{

cmpRECVD = colnames(datRECVD())

updateSelectizeInput(session, 'respMOD',
  choices = choices, selected = selected)

updateSliderInput(session,'scorethreshMOD',
  label = NULL, min = 0, max = 5.0,
  value = 1.0, step = 0.1)

updateSliderInput(session,'fdrMOD',label = NULL, min = 0, max = 0.5,
  value = 0.2, step = 0.05)

  },ignoreInit = F, ignoreNULL = T)

#######################################################
################## END OBSERVEEVENTS ##################
#######################################################

#########################################
############# START GO ENRICHMENT #######
#########################################

goEnrich <- reactive({

  ns <- session$ns
  req(input$fdrMOD)
  req(input$scorethreshMOD)

  req(input$respMOD)
#  req(respRECVD())
 # req(datRECVD())
  x <-  xsig

  req(input$poolset)


we = which(rownames(x) %in% ess$sgd_gene)
wne = which(rownames(x) %in% noness$sgd_gene)

  w = which(colnames(x) %in% input$respMOD)
#
#   w = which(colnames(x) %in% respRECVD())
#   validate(need(length(w) != 0, "Please select a compound"))

  req(length(w)!=0)

  xinp = x

  thresh = input$scorethreshMOD

  w2 = which(xinp[,w] >= thresh)

  validate(need(length(w2)!=0, message = "No scores above threshold"))

  req(length(w2)!=0)

  df = compSCORE(mat = xinp,coln = colnames(xinp)[w],sig = thresh)

  print(c("input$respMOD",input$respMOD))

 #  w = which(dfsig$signature %in% respRECVD())

  w = which(dfsig$signature %in% input$respMOD)

  ww = which(df$gene%in%dfsig$gene[w])

  m = match(df$gene[ww],dfsig$gene[w])
  df$score[ww] = 0
  df$score[ww] = dfsig$score[w][m]

  curr_exp = "network"

  FDR = input$fdrMOD

  network = runGOENRICH(fdrThresh = FDR,
    curr_exp =  curr_exp,score = df,
    bp_path = "2023_January30_GOBP.RDS",
    go_path = "2023_January30_GOID_GOBP_SGD.txt")

  enrichInfo = network$enrichInfo

  validate(need(!is.null(network$enrichInfo), message = message))

  return(network)

})

####################################
############# VISNET ###############
####################################

visNet <- reactive({

  req(goEnrich()$enrichInfo)

  enrich = goEnrich()$enrichInfo

  edge = goEnrich()$edgeMat

  vis = visSetup(enrichInfo = enrich,edgeMat = edge, fontsize = 18, fontface = "Courier")

  vis

  })

##################################################
################# START PROXY ####################
##################################################

output$network_Proxy <- renderVisNetwork({

  req(visNet()$nodes)

  vis = visNet()

  n = visNet()$nodes
  req(n)

  w =  nrow(n)
  n <- n %>% arrange(term)

  names = n$id

  if(nrow(vis$edges)==0) {

    visNetwork(vis$nodes, width = "100%") %>%

      visNodes(shadow=list(enabled=T,size=25),borderWidth=1) %>%

      visOptions(

        highlightNearest = list(enabled = T, degree = 5, hover = T),

        nodesIdSelection = list(enabled = TRUE, values = names,
          style = 'width: 500px; height = 31px; font-size: 18px;
          color: #000066;border: 3px solid #4d88ff;'),

        selectedBy = list(variable="FDR",
          style = 'width: 500px; height = 31px; font-size: 18px; c
             olor: #000066;border: 3px solid #4d88ff;')) %>%

      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);;}")

  } else {

    visNetwork(vis$nodes, vis$edges, width = "100%") %>%

      visNodes(shadow=list(enabled=T,size=25)) %>%

      visOptions(

        highlightNearest = list(enabled = T, degree = 5, hover = T),


        nodesIdSelection = list(enabled = TRUE, values = names,
                                style = 'width: 500px; height = 31px; font-size: 18px;
            color: #000066;border: 3px solid #4d88ff;'),

        selectedBy = list(variable="FDR",
                          style = 'width: 500px; height = 31px; font-size: 18px;
            color: #000066;border: 3px solid #4d88ff;')) %>%

      visIgraphLayout(type = "full") %>%

      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);;}")
  }

})

observeEvent(input$network_Proxy_selectedBy, {

  req(input$network_Proxy_selectedBy)

  n = visNet()$nodes

  w = which(n$FDR %in% as.numeric(input$network_Proxy_selectedBy))

  id = n$id[w]


  visNetworkProxy("network_Proxy") %>%

    visSelectNodes(id = id)

}, ignoreNULL = T, ignoreInit = FALSE)

##################################################
################### END PROXY ####################
##################################################

##################################################
################### END PROXY ####################
##################################################

##################################################
################  START  GO ENRICHMENT  ##########
##################################################

output$Gotable = renderUI({

req(input$network_Proxy_selected)

ns <- session$ns

DT::dataTableOutput(ns("detail"))

  })

#######################################################
################## END OBSERVEEVENTS ##################
#######################################################

output$detail <- DT::renderDataTable({

req(input$network_Proxy_selected)

vis = visNet()

n = visNet()$nodes

w = which(vis$nodes$id %in% c(input$network_Proxy_selected))

req(length(w) > 0)

term = vis$nodes$label[w]

nam = c("term","nGenes","geneSetFraction","FDR")

m = match(nam,names(vis$nodes))

n = vis$nodes[w,m]

term = vis$nodes$label[w]

nam = c("term","nGenes","geneSetFraction","FDR")

m = match(nam,names(vis$nodes))

names(vis$nodes)[m] = c("GO term","geneSet size","% of geneSet","FDR")

n = vis$nodes[w,m]

req(nrow(n)!=0)

t = t(n[,2:4])

datatable(t,width=220,selection = "single",caption = htmltools::tags$caption(term,
   style = "caption-side: top; text-align: center; color:black;background:white;font-weight:bold;"),
 options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F,rowCallback = JS(
   "function(row, data) {",
   "for (i = 1; i < data.length; i++) {",
   "if (data[i]>1000 | data[i]<1){",
   "$('td:eq('+i+')', row).html(data[i].toExponential(1));",
   "}",
   "}",
   "}")),
  height = 400,colnames = "") %>%
  formatStyle( target = "row", color = "black",backgroundColor = "white",
  columns = c(1,2),fontWeight = "bold")
  })
##################################

############################

##################################
#############   NETWORK   ###############
#########################################

#########################################
#############   NETWORK   ###############
#########################################
############################
############################
output$structure = renderUI({
 req(input$node_proxy_selected)
 ns <- session$ns

 htmlOutput(ns("chemistry"))

   })
############################
############################
output$chemistry = renderText({

 req(input$node_proxy_selected)

 n = nodes

 w = which(nodes$id %in% c(input$node_proxy_selected))

 req(length(w) > 0)

 n = nodes[w,]

 validate(need(nrow(n)!=0, message = "click signature or structure node for detail"))


 validate(need(n$img!="", message = "select fragment node to view structure"))

 struct = n$img
 h3 = paste0("<img src=",struct,">")
 h3
   })

###########################
##### START signature plot output
###########################
 output$sigOutput = renderUI({

 ns <- session$ns
 plotOutput(ns("barSigPlot"), height = 275)

   })
 output$LeadingEdge = renderUI({

 ns <- session$ns
 plotOutput(ns("barSigPlot1"), height = 275)

   })
#########################################
#############   NETWORK   ###############
#########################################

   output$barSigPlot <- renderPlot({

 ns <- session$ns

 w = which(dfsig$signature %in% input$respMOD)

 validate(
   need(length(w)!= 0, "Please choose a signature")
 )

 d = dfsig[w,]
 d$score = as.numeric(d$score)
 d = d %>% arrange(desc(score))
 nrows = nrow(d)

 if(nrows > 10) d = d[1:10,]

 d = d %>% arrange(score)

 tit = stringWINDOW(d$signature[1], width = 30)

 barplot(d$score,names.arg = d$gene,las=1,horiz=T,col="navy",main = tit, xlab = "median fitness defect score")

   })

  output$barSigPlot1 <- renderPlot({
req(input$network_Proxy_selected)

vis = visNet()

w = which(vis$nodes$id %in% c(input$network_Proxy_selected))

req(length(w) > 0)

n = vis$nodes[w,]

validate(need(nrow(n)!=0, message = "click node for detail"))

maxGenes=10
oGeneStr=n$overlapGenes
oGenes <- unlist(strsplit(oGeneStr, "\\|"))
genes <- strsplit(oGenes, "\\(")
scores <- sapply(genes, function(vec) { vec[length(vec)] })
genes <- sapply(genes, function(vec) { paste(vec[-length(vec)], collapse="(") })
scores <- unlist(strsplit(scores, ")"))
scoreMat <- data.frame(score=as.numeric(scores), stringsAsFactors=F)
rownames(scoreMat) <- genes
scoreMat$gene <- genes
if (nrow(scoreMat) > maxGenes) {
scoreMat <- scoreMat[1:maxGenes, ]
}

scoreMat <- scoreMat %>% arrange(score)
plotCol="#DACCFF"


tit = stringWINDOW(n$term, width = 25)

barplot(scoreMat$score,main = tit, horiz=T,col=plotCol,las=1,names.arg=scoreMat$gene, xlab = "Fitness defect score")


   })
#########################################
#############   NETWORK   ###############
#########################################
  hgtSigPlot = reactive({

req(input$respMOD)

w = which(dfsig$signature %in% input$respMOD)

validate(
  need(length(w)!= 0, "Please choose a signature")
)

d = dfsig[w,]

d = d %>% arrange(score)

nrows = nrow(d)

if(nrows > 10) d = d[1:10,]

height = barplotHEIGHT(d)

height = height*2

})
###########################
##### END signature plot output
###########################

 ########################################################################
 ################  END  GO NETWORK, PLOTS & DATATABLES  #################
 ########################################################################

 #################################################
 #################################################
 ########## START GO ENCHICHMENT TABLE ###########
 #################################################
 #################################################

   enrid = reactive({

req(visNet()$nodes)

enrich = visNet()$nodes

req(length(nrow(enrich)) > 0)

row <- input$enrichSig_rows_selected

out = outenrich()

id = out$id[row]

id


  })


############################
outenrich = reactive({

req(input$respMOD)
req(goEnrich()$enrichInfo)

enrich = goEnrich()$enrichInfo

w = which(names(enrich) %in% c("querySetFraction", "geneSetFraction" ,
   "foldEnrichment" , "P" , "FDR" ))
enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")] =
  format(round(enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")],2),
nsmall = 1,scientific = F)

enrich[,c("P", "FDR")] =
  format(signif(enrich[,c("P", "FDR")],2),nsmall = 1,scientific = T)

enrich
  }
  )
 ############################
  observeEvent(enrid(), {
req(enrid())
id = enrid()

ns <- session$ns
visNetworkProxy(ns("network_Proxy")) %>%
  visSelectNodes(id = id)
  })
############################

############################
output$enrichSig = DT::renderDataTable({

out = outenrich()
w = which(names(out) %in% c("GOID","term","querySetFraction", "geneSetFraction" ,
  "foldEnrichment" , "P" , "FDR","overlapGenes" ))

out = out[,c("GOID","term","querySetFraction", "geneSetFraction" ,
"foldEnrichment" , "FDR","overlapGenes" )]


opts = list(pageLength = 10,
autoWidth = T,scrollX = T,columnDefs = list(
 list(className = 'dt-left',targets = c(0,1,6)),
 list(className = 'dt-center',targets = c(2:5)),
 list(width = c('40px'),targets = c(2:4)),
 list(width = c('300px'),targets = c(6)),
 list(width = c('50px'),targets = c(0,5)),
 list(width = c('400px'),targets = 1)
   ))

df =  DT::datatable(out, options = opts,rownames = F, escape = F, selection = "single") %>%
  formatStyle(c(1:7),fontWeight = 'bold', fontSize = '14px')
 })


 ############################
 ############################ END GO ENCHICHMENT TABLE and downloads
 ############################

output$downloadSignatures <- downloadHandler(
filename = function() {
  paste0("response_signatures:", Sys.Date(), ".txt")
},
content = function(file) {
  write.table(dfsig, file, row.names = F,sep="\t")
}
  )
 #################################################
 ########## END GO ENCHICHMENT TABLE #############
 #################################################

####


 ################################################
 ###### START UPDATES & ObserveEVENTS ###########
 ##################
  returnedMOD = reactiveValues(
sig = NULL,
tab = NULL
  )

 ################################################
 ###### END UPDATES & ObserveEVENTS #############
 ################################################



  return(

  list(
sig = reactive({returnedMOD$sig }),
tab =  reactive({ returnedMOD$tab })
  )

)
}















 ############################
 ############################
 ############################

##################################
##################################
##########   MAP   ##############
#################################
#################################
#################################




