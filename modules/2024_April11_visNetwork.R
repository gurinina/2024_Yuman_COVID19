
source(file.path("modules","2022_March13_GOENRICH.R"))
source(file.path("modules", "2024_June19_hgeneAnnofunctions.R"))


phop        <-      readRDS("24Jan8/phop.RDS")
pelke       <-      readRDS("24Jan8/pelke.RDS")
pmarjhip    <-      readRDS("24Jan8/phip.RDS")

d3          <-      readRDS("24Jan8/d4.RDS")
dn3         <-      readRDS("24Jan8/dn4.RDS")
de3         <-      readRDS("24Jan8/de4.RDS")
delke       <-      readRDS("24Jan8/delke.RDS")
dint        <-      readRDS("24Jan8/dint.RDS")

###### 940 pixels = 12 columns

visNetworkModuleUI <- function(id, label = "visNetwork") {

  ns <- NS(id)

  tagList(

    fluidPage(title = "GO enrichments",

      fluidRow(
                column(width  = 3,

                       box(title = "Selected input datasets:",
                           br(),

      prettyRadioButtons("inpsite", label = NULL,
        choices = c(
          "2008 Erikson PLOS Genetics" = "elke",
          "2021 HIPHOP Marjan" = "marjhip",
          "2020:2024 HOP Marjan" = "marjhop"),

      outline = TRUE, fill = FALSE, thick = TRUE, shape = "square",
      bigger = FALSE, selected = "marjhop", inline = FALSE),

  status = "primary", solidHeader = T, width = "100%", height = 190)),


             column(width  = 5,
               box(title = "Select screen:",
                   br(),
                   selectizeInput(ns('cmpMOD'),
                                  '', choices = NULL, multiple = F
                   ),
                   status = "primary", solidHeader = T,width="100%", height = 200)
                ),



  column(width = 4,

      tabBox(title = tagList(shiny::icon("info-circle"), "GO enrichment legend"),

      tabPanel("GO network","To gain insight into b
        biological processes, Gene Ontology (GO) Analysis determines
        whether defined genesets or GO Terms ab
        out biological processes are over-represented within the
        geneset of interest, in this case the set
        of genes passing the fitness score threshold.
        To increase specificty, GO Terms were restr
        icted to genesets composed of at least 5 but
        not more than 300 genes. The significance of t
        he enrichments was estimated
        using the hypergeometric test."),
                 tabPanel("Leading edge genes","The enrichment
        of a GO gene set is driven by the gene subset
        present in the query set derived from the genes
        passing the fitness score threshold.
        In some cases, a similar set of query set genes
        drives the enrichment of multiple
        GO gene sets. Of these GO sets, we highlight the
        most significantly enriched,
        and the rest are considered redundant. The barplot
        shows the top contributing genes
        driving the enrichment for a given GO term
        (represented as a node in the network).
        The top 10 genes are shown."),
        width = NULL,height = 200))),


              fluidRow(

                column(width = 3,
                       box(title = "Reset menu:", align = "center",
                           br(),
                           actionButton(ns("resetGO"),label = "reset cmpMenu"),

                           solidHeader = TRUE, width = "100%",height = 170, status = "primary")#box
                ),


                column(width = 5,
                       box(title = "Fitness threshold:",#input score threshold = 3 ~ P < 0.001
                           sliderInput(ns("scorethreshMOD"),label  = NULL, min = 0, max = 5.0,
                                       value = 1.0, step = 0.1,ticks = TRUE),
                           status = "primary", solidHeader = T,width="100%", height = 170)
                ),#coln


                column(width = 4,
                       box(title = "FDR threshold:",
                           sliderInput(ns("fdrMOD"),label = "(FDR default = 0.2)", min = 0, max = 0.5,
                                       value = 0.2, step = 0.1, ticks = TRUE),
                           status = "primary", solidHeader = T,width="100%", height = 170)
                )),

              fluidRow(

                box(title =  "GO enrichment network:",
                    h5("select node to view GO term set details"),
                    visNetworkOutput(ns("network_proxy"),width = "100%",height = 652),
                    width = 7, status = "primary", solidHeader = TRUE, height = 800),

                box(title = "GO enrichment details:",width = 5,
                    uiOutput(ns("goTable")),solidHeader = T,
                    status = "primary",background = "navy", height = 310),


                box(title = "Leading edge genes:",width = 5,
                    column(width = 12, align = "center",
                    uiOutput(ns("LeadingEdge"))),solidHeader = T,
                    status = "primary", height = 470)
              ),

            fluidRow(
              box(title = "GoSet Enrichment Table:", DT::dataTableOutput(ns("enrichTable")),
                  status = "primary", solidHeader = TRUE,width = 12)
            ),

            fluidRow(
              box(title ="GO term enrichment table:",
                  column(width  = 12, align = "center",
                         br(),
                         downloadButton(ns('enrichdownload'), 'download GO enrichment table')),
                  status = "primary", solidHeader = T,width = 4)


          )
  ))
}

visNetworkModule = function(input,output,session, xinput,cmpInput,tabsInput, threshInput, siteInput,message = "No GO enrichment, try relaxing the FDR or scorethreshold"){

#####################   REACTIVE VARIABLES FROM SERVER   ############
   cmpRECVD = reactive({cmpInput()})
   datRECVD = reactive({xinput()})
   tabRECVD = reactive({tabsInput()})
threshRECVD = reactive({threshInput()})
  siteRECVD = reactive({siteInput()})

########################################################################
########################## START OBSERVEEVENT ##########################
########################################################################

observeEvent({datRECVD()
							threshRECVD()
							cmpRECVD()},{
  req(datRECVD())
  req(!is.null(tabRECVD()))
  req(cmpRECVD())
  req(siteRECVD())
ns <- session$ns
tabSERV = tabRECVD()
threshSERV = threshRECVD()
threshMOD = input$scorethreshMOD

cmpSERV = colnames(datRECVD())

selected = cmpRECVD()
site = siteRECVD()

tabNOTgoenrich = tabSERV!= "goenrich"

if(tabNOTgoenrich){
updateSelectizeInput(session, 'cmpMOD', choices = cmpSERV, selected = selected)

updateSliderInput(session,'scorethreshMOD',label = NULL,
    min = 0, max = 5.0,value = threshSERV, step = 0.1)

updateSliderInput(session,'fdrMOD',label = NULL,
    min = 0, max = 0.5,value = 0.2, step = 0.01)

updatePrettyRadioButtons(session,"inpsite", label = NULL,
    c("2008 Erikson PLOS Genetics" = "elke",
      "2021 HIPHOP Marjan only" = "marjhip",
      "2020:2023 HOP Marjan only" = "marjhop"),
       selected = siteRECVD())
}

},ignoreInit = F, ignoreNULL = F)

  observeEvent(input$cmpMOD, {
  ns <- session$ns

  tabINP = tabRECVD()

  req(tabINP == "goenrich" | tabINP == "hiphop" |
        tabINP == "genebyscreen"| tabINP == "signature" )

  choices = input$cmpMOD

  updateSliderInput(session,'scorethreshMOD',label = NULL,
    min = 0, max = 5.0,value = 1.0, step = 0.1)

  updateSliderInput(session,'fdrMOD',label = NULL,
    min = 0, max = 0.5,value = 0.2, step = 0.01)

  if(tabINP == "goenrich" ){

    returnMOD$thresh = input$scorethreshMOD
    returnMOD$cmp = input$cmpMOD
    returnMOD$site = input$inpsite
  }
},ignoreInit = F, ignoreNULL = F)

  observeEvent(input$scorethreshMOD,{
  req(threshRECVD())
  ns <- session$ns
  req(tabRECVD() == "goenrich" | tabRECVD() == "hiphop" | tabRECVD() == "genebyscreen")
  threshSERV = threshRECVD()
  tabGO = tabRECVD() == "goenrich"
  req(length(tabGO) >0)

  if(tabGO){

    returnMOD$thresh = input$scorethreshMOD
    returnMOD$cmp = input$cmpMOD
    returnMOD$site = input$inpsite
  }
},ignoreInit = F, ignoreNULL = T)

observeEvent(siteRECVD(), {

  req(datRECVD())
  req(cmpRECVD())
  drug = colnames(datRECVD())
  req(siteRECVD())
  cmp = cmpRECVD()
  site = siteRECVD()

  print(c("cmp",cmp))

  site = req(input$inpsite)

  if(site == "marjhop") {choices = phop$name} else if(site == "marjhip"){
    choices = pmarjhip$name} else if(site == "elke"){choices = pelke$name}

  w = which(cmp %in% phop$name)

  wj = which(cmp %in% pmarjhip$name)

  we = which(cmp %in% pelke$name)

  if(length(w) > 0){

    choices = phop$name} else if(length(wj) > 0) {

    choices = pmarjhip$name} else if(length(we)==0){

    choices = pelke$name}

  if(site == "marjhop") {

    choices = phop$name} else if(site == "marjhip"){

    choices = pmarjhip$name} else if(site == "elke"){

    choices = pelke$name}

  if(site == "marjhop") {

    selected ='24Jan05:methyl-methanesulfonate-mms_800uM'} else if(site == "marjhip"){

    selected ="23May09:doxorubicin-hydrochloride_13uM"} else if(site == "elke"){

    selected = "sertraline hydrochloride"}

    updateSelectizeInput(session,"cmpMOD",

       label = NULL,choices = drug,selected = selected,server = T)

    updatePrettyRadioButtons(session,"site", label = NULL,
      c("2008 Erikson PLOS Genetics" = "elke",
        "2021 HIPHOP Marjan only" = "marjhip",
        "2020:2023 HOP Marjan only" = "marjhop"), selected = siteRECVD())

}, ignoreNULL = TRUE, ignoreInit = FALSE)


observeEvent(input$site, {

  site = input$site

  if(site == "marjhop") {choices = phop$name} else if(site == "marjhip"){

    choices = pmarjhip$name} else if(site == "elke"){choices = pelke$name}

  if(site == "marjhop") {selected ='24Jan05:methyl-methanesulfonate-mms_800uM'} else if(site == "marjhip"){

    selected ="23May09:doxorubicin-hydrochloride_13uM"} else if(site == "elke"){

    selected = "sertraline hydrochloride"}

    print(c("site",site))

   updateSelectizeInput(session,"cmpMOD",
     label = NULL,choices = choices, selected = selected, server = T)

}, ignoreNULL = TRUE, ignoreInit = FALSE)


observeEvent(input$resetGO,{

    site=input$site

    req(site)

  if(site == "marjhop") {
    choices = phop$name} else if(site == "marjhip"){
    choices = pmarjhip$name} else if(site == "elke"){
    choices = pelke$name}

  if(site == "marjhop") {selected ='24Jan05:methyl-methanesulfonate-mms_800uM'} else if(site == "marjhip"){
    selected ="23May09:doxorubicin-hydrochloride_13uM"} else if(site == "elke"){
    selected = "sertraline hydrochloride"}

   updatePrettyRadioButtons(session,"site", label = NULL,
      c("2008 Erikson PLOS Genetics" = "elke",
        "2021 HIPHOP Marjan only" = "marjhip",
        "2020:2023 HOP Marjan only" = "marjhop"), selected = site)

   updateSelectizeInput(session,'cmpMOD',label = NULL,choices = choices,
        selected = selected)

   updateSliderInput(session,'scorethreshMOD',label = NULL,
    min = 0, max = 5.0,value = 1.0, step = 0.1)

   updateSliderInput(session,'fdrMOD',label = NULL,
    min = 0, max = 0.5,value = 0.2, step = 0.01)

  },ignoreInit = F, ignoreNULL = T)

#######################################################
#######################################################
#################### END OBSERVEEVENTS ################
#######################################################
#######################################################

####################################      NETWORK     #########################################
goEnrich <- reactive({

  x = datRECVD()
  w = which(colnames(x) %in% cmpRECVD())

  validate(
    need(length(w) != 0, "Please select a compound")
        )

  thresh = input$scorethreshMOD

  w2 = which(x[,w] >= thresh)

  validate(need(length(w2)!=0, message = "No scores above threshold"))

  req(length(w2)!=0)

  df = compSCORE(mat = x,coln = colnames(x)[w],sig = thresh)

  curr_exp = colnames(x)[w]

  FDR = input$fdrMOD

  network = runGOENRICH(fdrThresh = FDR, curr_exp =  curr_exp,score = df,
    bp_path = "2024_February26_GOBP.RDS",go_path = "2024_February26_GOBP_DF.txt")

  validate(
    need(!is.null(network$enrichInfo), message = message)
  )

  enrichInfo = network$enrichInfo

  req(!(is.null(enrichInfo)))

  edgeMat = network$edgeMat

  return(network)
  })

####################################      NETWORK     #########################################
visNet <- reactive({

  req(goEnrich()$enrichInfo)

  enrich = goEnrich()$enrichInfo

  edge = goEnrich()$edgeMat

  vis = visSetup(enrichInfo = enrich,edgeMat = edge, fontsize = 20, fontface = "Courier")


  vis

  })

####################################    NETWORK   #########################################
output$network_proxy <- renderVisNetwork({
  req(visNet()$nodes)
  vis = visNet()
  ns <- session$ns
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
       style = 'width: 500px; height = 31px; font-size: 18px;
       color: #000066;border: 3px solid #4d88ff;')) %>%

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

####################################      NETWORK     #########################################
observeEvent(input$network_proxy_selectedBy, {

  req(input$network_proxy_selectedBy)

  n = visNet()$nodes

  w = which(n$FDR %in% as.numeric(input$network_proxy_selectedBy))

  id = n$id[w]

  ns <- session$ns
  visNetworkProxy(ns("network_proxy")) %>%
    visSelectNodes(id = id)

               }
  )

output$goTable = renderUI({

    req(input$network_proxy_selected)

    ns <- session$ns

    DT::dataTableOutput(ns("gotermTable"))

  })

output$gotermTable <- DT::renderDataTable({

    req(input$network_proxy_selected)

    vis = visNet()

    n = visNet()$nodes

    w = which(vis$nodes$id %in% c(input$network_proxy_selected))

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

    datatable(t,width=220,caption = htmltools::tags$caption(term,
         style = "caption-side: top; text-align: center;
         color:black;background:white;font-weight:bold;"),

        options=list(paging=F,scrollY=F,dom="t",
          scroller=F,searching=F,ordering=F,rowCallback = JS(
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
########################      GO TERM TABLE     #########################################
output$LeadingEdge = renderUI({

     ns <- session$ns
     plotOutput(ns("barSigPlot"), height = hgt())

   })

output$barSigPlot <- renderPlot({

req(input$network_proxy_selected)

vis <- visNet()

w <- which(vis$nodes$id %in% c(input$network_proxy_selected))

req(length(w) > 0)

n <- vis$nodes[w, ]

validate(need(nrow(n) != 0, message = "click node for detail"))

# automatically 10 rows and arrange by score
scoreMat <- geneBARPLOT(n$overlapGenes)

plotCol <- "#DACCFF"

tit <- stringWINDOW(n$term, width = 25)

barplot(scoreMat$score,
main = tit, horiz = T, col = plotCol,
las = 1, names.arg = scoreMat$gene, xlab = "Fitness defect score"
)
})

hgt = reactive({

  req(input$network_proxy_selected)

  vis = visNet()

  n = visNet()$nodes

  w = which(vis$nodes$id %in% c(input$network_proxy_selected))

  req(length(w) > 0)

  n = vis$nodes[w,]

  validate(need(nrow(n)!=0, message = "click node for detail"))

  o = geneBARPLOT(n$overlapGenes)

  height = genebarHEIGHT(o)

  height = height*2


  height
})

###########################  GO ENRICHMENT TABLE  ###################################
  enrichReact = reactive({

    req(visNet()$nodes)

    enrich = visNet()$nodes

    req(length(nrow(enrich)) > 0)

    row <- input$enrichTable_rows_selected

    out = outEnrich()

    id = out$id[row]

    id


  })

###########################  GO ENRICHMENT TABLE  ###################################
outEnrich = reactive({

    req(input$cmpMOD)
    req(goEnrich()$enrichInfo)

    enrich = goEnrich()$enrichInfo

  #########
  #########
    enrich[,c("querySetFraction","geneSetFraction",

      "foldEnrichment")] = format(round(enrich[,c("querySetFraction","geneSetFraction",

      "foldEnrichment")],2),nsmall = 1,scientific = F)

    enrich[,c("P", "FDR")] = format(signif(enrich[,c("P", "FDR")],2),

       nsmall = 1,scientific = T)

    w = which(names(enrich) %in% c("formattedLabel","cluster",

      "size","filename","maxOverlapGeneScore"))

    enrich = enrich[,-w]

    enrich = enrich[,c("id",  "GOID","term","nGenes","nQuery",
       "nOverlap", "querySetFraction","geneSetFraction",
       "foldEnrichment","P","FDR","overlapGenes")]

    enrich
  }
  )

###########################  GO ENRICHMENT TABLE  ###################################
observeEvent(enrichReact(), {
    req(enrichReact())
    id = enrichReact()

    ns <- session$ns
    visNetworkProxy(ns("network_proxy")) %>%
      visSelectNodes(id = id)
  })

###########################  GO ENRICHMENT TABLE  ###################################

output$enrichTable = DT::renderDataTable({

  out = outEnrich()

  out = out[,c("GOID","term","querySetFraction", "geneSetFraction" ,
  "foldEnrichment" , "FDR","overlapGenes" )]

  opts = list(pageLength= 10,
        autoWidth = T,scrollX = T,

    columnDefs = list(
    list(className = 'dt-nowrap',targets = c(1,5))),
    list(className = 'dt-left',width='100%',targets = c(0:6)))

  df =  DT::datatable(out, options = opts,rownames = F,
    escape = F, selection = "single") %>% formatStyle(c(1:7),fontWeight = 'bold')

 })

###########################  GO ENRICHMENT TABLE  ###################################
output$enrichdownload <- downloadHandler(
  filename = function() {
    paste0("enrich:",input$cmpMOD,"_", Sys.Date(), ".txt")
  },
  content = function(file) {
    write.table(outEnrich(), file, row.names = F,sep="\t")
  }
  )

#####################################################################################@####
##################### VARIABLES TO RETURN TO SERVER ##################################@####
#####################################################################################@####
returnMOD = reactiveValues(cmp = NULL,thresh = NULL, site = NULL)
##########################################################################################
  return(
    list(
    cmp = reactive({ returnMOD$cmp }),
    thresh =  reactive({ returnMOD$thresh }),
    site = reactive({ returnMOD$site })
    )
  )

}

