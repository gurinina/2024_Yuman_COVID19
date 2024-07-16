library(shiny)
library(dplyr)
library(shinyWidgets)
library(DT)
library(shinyjs)
library(shinydashboard)
library(bcrypt)
library(visNetwork)
library(shinycssloaders)
library(igraph)
library(ggplot2)
library(ggrepel)
library(shinyjs)
library(shinyBS)
library(stringr)
library(shinyjqui)
#######

options(rsconnect.max.bundle.size=3145728000)

source(file.path("modules", "2022_March13_GOENRICH.R"))
source(file.path("modules" ,"2024_June19_hgeneAnnofunctions.R"))
source(file.path("modules" ,"2024_March22_functions.R"))
source(file.path("modules" ,"2024_June3_signature.R"))
source(file.path("modules" ,"2024_April3worksyeast_visNetwork.R"))

xsig        <-      readRDS("24Jan8/x21.RDS")
orth        <-      readRDS("24Jan8/orthologs.RDS")
   w        <-      which(orth$hessential=="ess")
noness      <-      orth[-w,]
ess         <-      orth[w,]

d3          <-      readRDS("24Jan8/d4.RDS")
dn3         <-      readRDS("24Jan8/dn4.RDS")
de3         <-      readRDS("24Jan8/de4.RDS")
delke       <-      readRDS("24Jan7/delke.RDS")
dint        <-      readRDS("24Jan8/dint.RDS")

phop        <-      readRDS("24Jan8/phop3.RDS")
pelke       <-      readRDS("data/pelke.RDS")
pmarjhip    <-      readRDS("24Jan8/phip.RDS")

xcoinh      <-      readRDS("24Jan8/hcoinh.RDS")
xpval       <-      readRDS("24Jan8/hcoinh_pval.RDS")
xcofit      <-      readRDS("24Jan8/xcofit.RDS")
xpfit       <-      readRDS("24Jan8/xcofit_pval.RDS")

xcohh       <-      readRDS("24Jan8/hxhiphop_coinh.RDS")
xpvhh       <-      readRDS("24Jan8/hxhiphop_coinh_pval.RDS")
xfithh      <-      readRDS("24Jan8/xhiphop_cofit.RDS")
xpithh      <-      readRDS("24Jan8/xhiphop_cofit_pval.RDS")

xefitpv     <-    readRDS("24Jan8/helke_cofit_pval.RDS")
xefit       <-    readRDS("24Jan8/helke_cofit.RDS")
xelke       <-    readRDS("24Jan8/helke.RDS")
xelkepv     <-    readRDS("24Jan8/helke_pval.RDS")

d3          <-     d3[,phop$name]
de3         <-     de3[,pmarjhip$name]
delke       <-     delke[,pelke$name]

st          <-     rev(seq(0,1,0.05))

sq          <-     rev(seq(0,1,0.001))

cmp         <-     '24Jan05:methyl-methanesulfonate-mms_800uM'


header <- dashboardHeader(title =  span("Chemogenomic Profiling COVID19",
          style = "font-weight:bold; font-size:18px; text-align: right;"),
          titleWidth = 1200)

sidebar <- dashboardSidebar(width = 220,
   tags$style(HTML(".main-sidebar {
   background-color: #000066 !important; color: #fff !important;
   font-weight:bold !important}.treeview-menu>li>a {font-weight:
   color: #fff !important;}")),

sidebarMenu(id = 'tabs',selected = "HOmozygous Profiles",

menuItem("HOmozygous Profiles", tabName="hiphop",icon = icon("bullseye"), selected = T),

menuItem('GO enrichments', tabName = 'goenrich', icon = icon('dna')),

menuItem('Signatures', tabName = 'signature', icon = icon('database')),

menuItem("Gene fitness profiles", tabName="genebyscreen", icon = icon("chart-line"))

))

body <-  dashboardBody(

   tags$style(HTML(".main-sidebar { font-size: 16px;

   background-color: #000066 !important; color: #fff;
   font-weight:bold} .treeview-menu>li>a { font-size: 16px;
   font-weight:bold !important; color: #fff !important;}")),

  tags$style("


        #controls {
          background-color: #dcd0ff;
          opacity: 1;
        }
        #controls:hover{
          opacity: 1;
        }

        #econtrols {
          background-color: #dcd0ff;
          opacity: 1;
        }
        #econtrols:hover{
          opacity: 1;
          }"),

  tags$head(tags$style(HTML('
        .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
             background-color: #4d88ff;
             color: #ffffff;
         }
             '))),

  tags$head(tags$style(HTML('
        .skin-blue .main-header .logo {
             background-color: #000066;
             }
             .skin-blue .main-header .logo:hover {
             background-color: #4d88ff;;
             }
             '))),

  tags$head(tags$style(HTML('
        .skin-blue .main-header .navbar{
             background-color: #000066;
             }
             '))),

  tags$head(tags$style(HTML('
       navbar-custom-menu {
             background-color: #000066;
             }
             .skin-blue .main-header .logo:hover {
             background-color: #4d88ff;;
             }
             '))),

  tags$head(tags$style(HTML('
.skin-blue .main-header .navbar .sidebar-toggle {
    color: #fff
 background-color: #000066;}}'))),


tags$style("

        #controls {
          background-color: #dcd0ff;
          opacity: 1;}

        #controls:hover{
          opacity: 1;}

        #econtrols {
          background-color: #dcd0ff;
          opacity: 1;}

        #econtrols:hover{
          opacity: 1;}}"),


tags$head(tags$style(HTML("#point_einfo tr.selected {background-color:white},
  #point_info tr.selected {background-color:white}"))),

#makes pngs same size as window after rescaling
tags$script("$(document).on('shiny:connected', function(event) {
  var myWidth = $(window).width();
  Shiny.onInputChange('shiny_width',myWidth)  });"),

useShinyjs(),

tags$script("$(document).on('shiny:connected', function(event) {
  var myHeight = $(window).height();
  Shiny.onInputChange('shiny_height',myHeight)});"),


tags$script(func <- JS('function(event, ui){return $(event.target).offset();}')),

tags$script(
  HTML("// Get mouse coordinates
           var mouseX, mouseY;
           $(document).mousemove(function(e) {
           }).mouseover();

           // Function to position draggable, place on current mouse coordinates
           Shiny.addCustomMessageHandler ('placeDraggable',function (message) {
           var element = $('#click_info').parent();});

           // Show or hide draggable
          Shiny.addCustomMessageHandler ('hideDraggable',function (message) {
            if(message.hide == true){
              $('#click_info').parent().hide();} else{
              $('#click_info').parent().show();}});


         // Get mouse coordinates
           #var emouseX, emouseY;
           $(document).mousemove(function(e) {
           }).mouseover();

           // Function to position draggable, place on current mouse coordinates
           Shiny.addCustomMessageHandler ('eplaceDraggable',function (message) {
           var element = $('#eclick_info').parent();
          # element.css({'top': emouseY + 'px', 'left' : emouseX + 'px'})});

           // Show or hide draggable
          Shiny.addCustomMessageHandler ('ehideDraggable',function (message) {
            if(message.hide == true){
              $('#eclick_info').parent().hide();
            } else{
              $('#eclick_info').parent().show();}")),

tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "chemogenomics.css")),

tags$style(HTML(".selectize-input{white-space: nowrap;
    font-size: 14px !important; font-weight:bold!important;")),

tabItems(

  tabItem(tabName = "genebyscreen",

   fluidRow(

    column(width = 3,

      box(title = "Enter gene of interest:",

        textInput(inputId = "gene", label = NULL, value   = "RAD51"),

      status = "primary", solidHeader = T, width = "100%", height = 250)),

    column(width = 9,

      box(title = "Gene descriptor and function:",

        dataTableOutput('geneinfo'),

      status = "primary", solidHeader = T, width = "100%", height = 250))),

  fluidRow(column(width = 12,

      box(title = "Fitness defect score across conditions:",

        h5("click OR drag to select screens of interest;
          corresponding profiles can be found in compound
          menu om hiphop tab; cyan: p-value < 0.001; triangle:
          screens identifying drug target candidates"),

        plotOutput("genebydrug",
          click = clickOpts(id = "clickpts"),
          brush = brushOpts(id = "brushpts"), height = 800),

        status = "primary", solidHeader = T, width = "80%"))),

  fluidRow(

      box(title = "Experimental detail; select row to view
        corresponding HOP chemogenomic profile",

        dataTableOutput("tabpts"),

      status = "primary", solidHeader = T, width = 12)),

  fluidRow(

      box(title = "Cofit genes, select row to view corresponding
        gene fitness profile",

        dataTableOutput('coFitNess'),

     status = "primary", solidHeader = T, width = 12)),

   fluidRow(
      column(width = 6,
             uiOutput("mxFIT")),
      column(width = 6,
             uiOutput("fiToutput"))),


  fluidRow(
      box(title = "GO enrichment network: select a node to view details;
        right click to save image",
        width = 8,
        status = "primary",
        solidHeader = TRUE,
        height = 800,
        uiOutput("GOterms"),
        visNetworkOutput("network_Proxy", width = "85%", height = 652)
      ),

      box(
        title = "GO term set enrichment details:",
        width = 4,
        uiOutput("GoTable"),
        solidHeader = T,
        status = "primary",
        background = "navy",
        height = 200
      ),

      box(
        title = "Top-contributing genes:",
        width = 4,
        br(),
        uiOutput("leadingEdge"),
        solidHeader = T,
        status = "primary",
        height = 580
      )
    ),
    #fluidRow

  fluidRow(

    box(title = "GO enrichment table:",

      DT::dataTableOutput("EnrichTable"),

    status = "primary", solidHeader = T, width = 12)),

  fluidRow(

      column(width = 3, align = "center",

      box(title = "Download CoFitness:",

        br(),

        downloadButton("downloadcofitness", "HOP CoFitness"),

      status = "primary", solidHeader = T, width = "100%", height = 150)),

      column(width = 3,

        box(title = "Download EnRichments:", align = "center",

        br(),

        downloadButton("enrichdownload", "GO EnRichments"),

        status = "primary", solidHeader = T, width = "100%", height = 150)))),

  # ok
tabItem("hiphop", class = "active",

  fluidRow(

    column(width = 3,

      box(title = "Select input datasets:",

        br(),

        prettyRadioButtons("site", label = NULL,

          choices = c("2008 Erikson PLOS Genetics" = "elke",
                      "2021 HIPHOP Marjan" = "marjhip",
                      "2020:2024 HOP Marjan" = "marjhop"),

          outline = T, fill = F, thick = T, status = "primary", shape = "square",

            bigger = F, selected = "marjhop", inline = F),

      status = "primary", solidHeader = T, width = "100%", height = 160)),

column(width = 9,

        box(title = "Select screen:",

          br(),
          selectizeInput('cmpSERV',
          label = NULL, choices = phop$name, multiple = F,

          selected = cmp,
          options = list('plugins' = list('autofill_disable'))),

        status = "primary", solidHeader = T, width = "100%", height = 160))),

    fluidRow(

      column(width = 3,

        box(title = "Reset compound menu:", align = "center",

          br(),

          actionButton("resetCmp", "Reset cmpMenu"),

        status = "primary", solidHeader = T, width = "100%", height = 200)),













        column(width = 3,

          box(title = "Fitness score threshold:",

            fluidRow(column(width = 12, align = "center",

              sliderInput("scorethreshSERV", label = NULL, min = 0, value = 1.0,

                step = 0.1,  max = 5.0))),


#
  fluidRow(column(
              width = 12,
              align = "center",
              actionLink(
                inputId = "resetscore",
                label = "reset",
                style = "height: 100px;width: 150px;font-size:120%;text-align:center",
                size = "lg"
              )
            )),
            status = "primary",
            solidHeader = T,
            width = "100%",
            height = 200
          )
        ),
        #coln
        uiOutput("ylimits")
      ),#fluid row

      uiOutput("hiphoppanel")
    ),#tabItem


tabItem("signature",

  fluidRow(
    sigNetworkModuleUI("sigNetworkModule1")),

  fluidRow(
    box(title = "HOP profiles in response signature:",

     dataTableOutput("screenResp"),

     status = "primary", solidHeader = T, width = 12))),

tabItem("goenrich",

  fluidRow(

    visNetworkModuleUI("visNetwork1")))))#fluidRow tabItems
# ok


ui<-dashboardPage(header, sidebar, body, skin = "blue",
  setBackgroundColor(color = "#e6e6ff", shinydashboard = TRUE))

server <- function(input, output, session) {

############################

observeEvent(input$site, {

  req(input$site)

  site = input$site

  req(site)

  if(site == "marjhop") {
    choices = phop$name} else if(site == "marjhip"){
    choices = pmarjhip$name} else if(site == "elke"){
    choices = pelke$name}

  if(site == "marjhop") {

  selected = "24Jan05:methyl-methanesulfonate-mms_800uM"

} else if(site == "marjhip"){

    selected = "23May09:doxorubicin-hydrochloride_13uM"} else if(site == "elke"){

    selected =  "sertraline hydrochloride" }

    updateSelectizeInput(session,"cmpSERV", label = NULL,

    choices = choices,selected = selected,server = T)

}, ignoreNULL = TRUE, ignoreInit = FALSE)


  ###########################################################################
  ########################### START HIPHOPPANEL #############################
  ###########################################################################

  output$hiphoppanel  = renderUI({

    welke  = which(colnames(delke)%in% input$cmpSERV)
    whop   = which(colnames(d3) %in% input$cmpSERV)
    wboth  = which(colnames(dint) %in% input$cmpSERV)

    options <-  list(shiny = list(abs_position = list(
      dragcreate = func, # send returned value back to shiny when interaction is created.
      drag = func)))# send returned value to shiny when dragging.

    if(length(whop) > 0){

      tabItem("HOP",

              fluidRow(

                box(title = "Signature & compound information:",

                    HTML("<h5><b>Click datatable row to view response signature:</b></h5>"),

                    DT::dataTableOutput("targethip"),

                    status = "primary", solidHeader = TRUE, width = 12)), # fluidRow ok

              fluidRow(

                box(title = "HOP profile: mouse click on points to view gene description",

                    plotOutput("nfd", width = "90%", height = 800,

                               click = clickOpts(id = "plot_click")),

                    status = "primary", solidHeader = TRUE, width = 12),

                jqui_draggable(absolutePanel(id = "controls", width = 600,
                                             draggable = TRUE,uiOutput("click_info")),
                               options = list(cancel = ".shiny-input-container"))),# fluidRow ok

              fluidRow(

                box(title="HOP genes: mouse click on row links to gene fitness profile;

      mouse click on genename links to SGD",

      dataTableOutput("hoptab"),

      status = "primary", solidHeader = TRUE,width = 12), # fluidRow ok

      box(title="download HOP:",

          column(width =6,

                 downloadButton("downloadhop", "HOPdownload")),

          column(width = 6,

                 downloadButton("ExportnonPlot", "HOPfitness plot")),

          status = "primary", solidHeader = TRUE,width = 12)),

      fluidRow(

        box(title="Compound information: CoInhibitory (CoI) value,
      pvalue, PubChem ID (CID) link, Mechanism of Action (MoA) link,
      indication and mechanism",

      h5(strong("Mouse click on row to view CoInhibitory (CoI)
        Homozygous Profile (HOP) and corresponding GO enrichments in next tab")),

      br(),

      dataTableOutput('coInhib'),

      status = "primary", solidHeader = TRUE,width = 12)), # fluidRow ok

      fluidRow(align = "center",

               box(title="download CoInhibition:",

                   downloadButton('downloadcoinhib', 'CoInhibition'),

                   status = "primary", solidHeader = TRUE,width = 3))) # tabItem ok

    } else if(length(wboth) > 0){

      tabItem("HIPHOP",

              fluidRow(

                box(title = "HIP fitness profile:",

                    plotOutput("efd", width = "90%",height = 800,

                               click = clickOpts(id="plot_eclick")),

                    status = "primary", solidHeader = TRUE, width = 6),

                box(title = "HOP fitness profile:
      mouse click on points to view gene description",

      plotOutput("nfd", width = "90%", height = 800,

                 click = clickOpts(id="plot_click")),

      status = "primary", solidHeader = TRUE, width = 6),

      jqui_draggable(absolutePanel(id = "controls", width = 600,
                                   draggable = TRUE,uiOutput("click_info")),
                     options = list(cancel = ".shiny-input-container"))), # fluidRow oK

      fluidRow(

        box(title="HIP genes: mouse click on row links to gene fitness profile",

            dataTableOutput("hiptab", width = "100%"),

            status = "primary", solidHeader = TRUE,width = 6),

        box(title="HOP genes: mouse click on genename links to SGD",

            dataTableOutput("hoptab",width = "100%"),

            status = "primary", solidHeader = TRUE, width = 6)), # fluidRow oK

      fluidRow(

        box(title="download HIP:",

            column(width = 6,

                   downloadButton("downloadhip", "HIPdownload")),

            column(width = 6,

                   downloadButton("ExportessPlot", "HIPfitness plot")),

            status = "primary", solidHeader = TRUE,width = 6),

        box(title="download HOP:",

            column(width =6,

                   downloadButton("downloadhop", "HOPdownload")),

            column(width = 6,

                   downloadButton("ExportnonPlot", "HOPfitness plot")),

            status = "primary", solidHeader = TRUE,width = 12)), # fluidRow oK

      fluidRow(

        box(title="Compound information: CoInhibitory (CoI) value,
      pvalue, PubChem ID (CID) link, Mechanism of Action (MoA) link,
      indication and mechanism",

      h5(strong("Mouse click on row to view CoInhibitory (CoI)
        Homozygous Profile (HOP) and corresponding GO enrichments in next tab")),

      br(),

      dataTableOutput('coInhib'),

      status = "primary", solidHeader = TRUE,width = 12)), # fluidRow ok

      fluidRow(align = "center",

               box(title="download CoInhibition:",

                   downloadButton('downloadcoinhib', 'CoInhibition'),

                   status = "primary", solidHeader = TRUE,width = 3))) # tabItem ok

    } else if(length(welke) > 0){

      tabItem ("HOPELKE",

               fluidRow(

                 box(title = "HOP profile: mouse click on points to view gene description",

                     plotOutput("nfd", width = "90%", height = 800,

                                click = clickOpts(id = "plot_click")),

                     status = "primary", solidHeader = TRUE, width = 12),

                 jqui_draggable(absolutePanel(id = "controls", width = 600,
                                              draggable = TRUE,uiOutput("click_info")),
                                options = list(cancel = ".shiny-input-container"))),# fluidRow ok

               fluidRow(

                 box(title="HOP genes: mouse click on row links to gene fitness profile;

      mouse click on genename links to SGD",

      dataTableOutput("hoptab"),

      status = "primary", solidHeader = TRUE,width = 12), # fluidRow ok

      box(title="download HOP:",

          column(width =6,

                 downloadButton("downloadhop", "HOPdownload")),

          column(width = 6,

                 downloadButton("ExportnonPlot", "HOPfitness plot")),

          status = "primary", solidHeader = TRUE,width = 12)),

      fluidRow(

        box(title = "Compound information: CoInhibitory (CoI) value,
    pvalue, PubChem ID (CID) link, Mechanism of Action (MoA) link,
    indication and mechanism",

    h5(strong("Mouse click on row to view CoInhibitory (CoI)
        Homozygous Profile (HOP) and corresponding GO enrichments in next tab")),

    br(),

    dataTableOutput('coInhib'),

    status = "primary", solidHeader = TRUE,width = 12)), # fluidRow ok

    fluidRow(align = "center",

             box(title="download CoInhibition:",

                 downloadButton('downloadcoinhib', 'CoInhibition'),

                 status = "primary", solidHeader = TRUE,width = 3))) # tabItem ok

    }})

#####################################################################
########################## END HIPHOPPANEL ##########################
#####################################################################

######################################################################
###################### START CALL MODULES ############################
######################################################################

tabSEND = reactiveValues(

  tab = NULL

)

xinput = reactive({

  req(input$site)

  site = input$site

  if(site == "marjhop"){

    xinput = d3[,phop$name]

  } else if(site == "marjhip"){

    xinput = dint

  } else {xinput = delke}

})


pdata = reactive({

  req(input$site)

  site = input$site

  if(site == "marjhop"){

    pdata = phop

  } else if(site == "elke"){

    pdata = pmarjhip

  } else {pdata = pelke}

})

threshSEND = reactiveValues(thresh = NULL)

cmpSEND = reactiveValues(cmp = NULL)

respSEND = reactiveValues(resp = NULL)

observeEvent(input$tabs, {

  tabSEND$tab = input$tabs
  cmpSEND$cmp = input$cmpSERV

},  ignoreInit = F, ignoreNULL = T)

returnedMOD = callModule(visNetworkModule,'visNetwork1',xinput = xinput,
  cmpInput = reactive(cmpSEND$cmp),
  threshInput = reactive(threshSEND$thresh),
  siteInput = reactive(input$site),
  tabsInput = reactive(tabSEND$tab))

cmpRETURN = reactive({
  print(paste("cmpRETURNED",returnedMOD[[1]]()))
  returnedMOD[[1]]()})

threshRETURN = reactive({
  print(paste("threshRETURNED",returnedMOD[[2]]()))
  returnedMOD[[2]]()})

siteRETURN = reactive({
  print(paste("siteRETURNED",returnedMOD[[3]]()))
  returnedMOD[[3]]()})

returnedSIG <- callModule(
      sigNetworkModule,
      'sigNetworkModule1',
      xRespDat = reactive(xsig),
      xRespInput = reactive(respSEND$resp),
      inputTab  = reactive(tabSEND$tab),
      message = "No GO enrichment, try another signature")

sigRETURN = reactive({
    print(paste("sigRETURN", returnedSIG[[1]]()))
    returnedSIG[[1]]()})

tabRETURN = reactive({
    print(paste("tabservRETURN", returnedSIG[[2]]()))
    returnedSIG[[2]]()})

########################################################################
########################## END CALL MODULES ############################
########################################################################

########################################################################
########################## START OBSERVEEVENT ##########################
########################################################################

observeEvent(threshRETURN(), {

  req(threshRETURN())

  tabNOThop = input$tabs != "hiphop"

  threshMOD = threshRETURN()

  req(threshMOD)

  req(input$scorethreshSERV)

  req(threshMOD!="")

  req(!is.null(input$cmpSERV))

  threshDIFF = threshMOD!=input$scorethreshSERV

  if(threshMOD!=input$scorethreshSERV){

    updateSliderInput(session, "scorethreshSERV", label = NULL,

    value = threshMOD, step = 0.1, min = 0, max = 5.0)}

},ignoreInit = F, ignoreNULL = T)

observeEvent(cmpRETURN(),
               {
req(cmpRETURN())

tabNOThop = input$tabs != "hiphop"

threshMOD = threshRETURN()

cmpMOD = cmpRETURN()

req(cmpMOD)

req(cmpMOD!="")

req(!is.null(input$cmpSERV))

cmpDIFF = cmpMOD!=input$cmpSERV

req(tabNOThop == TRUE | tabNOThop == FALSE)

if(cmpMOD %in% phop$name){
  choices = phop$name}else if(cmpMOD %in% pmarjhip$name){
  choices = pmarjhip$name}else if(cmpMOD %in% pelke$name){
  choices = pelke$name}

if(cmpMOD %in% phop$name){
    site = "marjhop"} else if(cmpMOD %in% pmarjhip$name){
    site = "marjhip"} else if(cmpMOD %in% pelke$name){
    site = "elke"}

        updatePrettyRadioButtons(session,"site", label=NULL,
        c("2008 Erikson PLOS Genetics" = "elke",
          "2021 HIPHOP Marjan" = "marjhip",
          "2020:2024 HOP Marjan" = "marjhop"),selected = site)

if(cmpMOD!=input$cmpSERV){

  updateSelectizeInput(session,"cmpSERV",

  label = NULL, choices = choices,

  selected = cmpRETURN(),server = T)
}

},ignoreInit = F, ignoreNULL = T)

observeEvent(input$cmpSERV, {

    req(input$cmpSERV)

    req(!is.null(input$tabs))

    req(input$tabs!= "")

    updateSliderInput(session, "scorethreshSERV", "",
      value = 1.0, step = 0.1, min = 0, max = 5.0)

    cmpSEND$cmp = input$cmpSERV

    tabSEND$tab = input$tabs

    threshSEND$thresh = input$scorethreshSERV

  },ignoreInit = F, ignoreNULL = T)

  observeEvent(input$scorethreshSERV,{

    req(input$cmpSERV)

    cmpSEND$cmp = input$cmpSERV

    tabSEND$tab = input$tabs

    threshSEND$thresh = input$scorethreshSERV


  },ignoreInit = F, ignoreNULL = T)

#######################################################
################## END OBSERVEEVENTS ##################
#######################################################

#######################################################
############# START SIGNATURE EVENTS ##################
#######################################################

output$screenResp = renderDataTable({

    w = which(phop$signature %in% sigRETURN())

    validate(need(length(w) > 0 , message =
   "please enter a valid signature"))


    df = data.frame(

      screen = phop$name[w],

      mechanism = phop$mechanism[w],

      PCID = phop$PCID[w],


      compound = phop$cmp[w],

      stringsAsFactors = F
    )

     opts = list(

     pageLength = 25,

     autoWidth = F,

     scrollX = F,

     columnDefs = list(

     list(className = 'dt-left', targets = c(0, 1, 2, 3)),

    list(width = c('100px'), targets = c(2)),

    list(width = c('500px'), targets = c(1)),

    list(width = c('300px'), targets = 0)))

    df = df %>% dplyr::arrange(compound)

     df =  DT::datatable(

     df,

     options = opts, rownames = F, escape = F,
     selection = "single") %>%

    formatStyle(c(1:4), fontWeight = 'bold', fontSize = '14px')

       })


  observeEvent(input$screenResp_rows_selected, {

      w = which(phop$signature %in% sigRETURN())

         print(paste("sigselectedRETURN", sigRETURN()))

         validate(need(length(w) > 0 , message =
   "please enter a valid signature"))

    d = phop$name[w]

      df = data.frame(

        screen = phop$name[w],

         mechanism = phop$mechanism[w],

      PCID = phop$PCID[w],


      compound = phop$cmp[w],

      stringsAsFactors = F
    )

    df = df %>% dplyr::arrange(compound)

    row = input$screenResp_rows_selected

    w = which(df$screen %in% df$screen[row])

    validate(need(length(w) > 0 , message = "please enter a valid signature"))

    selected = df$screen[row]

    cmpSEND$cmp = selected

   updateSelectizeInput(session, "cmpSERV", label = NULL,
     choices = phop$name, selected = selected,server = T)

   newtab <- switch(input$tabs,

  "signature" = "hiphop",

   "hiphop" = "signature")

  updateTabItems(session, "tabs", newtab)

  })

output$targethip <- DT::renderDataTable({

  w = which(phop$name %in% input$cmpSERV)

  d = phop[w, c(
      "name",
      "cmp",
      "signature",
      "PCID",
      "image",
      "mechanism"
    )]

    PCID = '<a href="https://pubchem.ncbi.nlm.nih.gov" target="_blank">PCID</a>'

    names(d)[5] = "structure"

    names(d)[1] = "screen"

    names(d)[4] = PCID

    names(d)[2] = "compound"

    opts = list(
        dom = 'Bfrtip',
        paging = F,
        target = "cell",
        searching = F,
        info = F,
        autowidth = T,
        scrollX = TRUE,
        ordering = F,
        columnDefs = list(
          list(className = 'dt-left', targets = c(0,1,4)),
          list(className = 'dt-center', targets = c(2,3)),
          list(width = c('200px'), targets = c(0,2)),
          list(width = c('100px'), targets = c(1,3)),
          list(width = c('400px'), targets = c(5))))


    datatable(
      d,
      options = opts,
      escape = F,
      class = 'table-bordered stripe table-condensed',
      rownames = F
    )
  })

  observeEvent(input$targethip_rows_selected, {

    w = which(phop$name %in% input$cmpSERV)

    row = input$targethip_rows_selected

    respSEND$resp = phop$signature[w][row]

    newtab <- switch(input$tabs,

    "hiphop" = "signature",

    "signature" = "hiphop")

    updateTabItems(session, "tabs", newtab)

  })


#######################################################
############### END SIGNATURE EVENTS ##################
#######################################################

############# RESET COMPOUND #############

 observeEvent(input$resetCmp,{

   site = input$site

   req(input$site)

   updateSelectizeInput(session,'cmpSERV',label = NULL,choices = NULL,
         selected = "")

  if(site == "marjhop"){

   updateSelectizeInput(session,"cmpSERV",

      label = NULL, choices = phop$name,

      selected = "24Jan05:methyl-methanesulfonate-mms_800uM",server = T)

  } else if(site == "marjhip"){

  updateSelectizeInput(session,"cmpSERV", label = NULL,

      choices = pmarjhip$name,

      selected = "23May09:doxorubicin-hydrochloride_13uM",server = T)

     } else if(site == "elke"){

  updateSelectizeInput(session,"cmpSERV",

      label = NULL, choices = pelke$name,

      selected = "sertraline hydrochloride",server = T)}

  },ignoreInit = F, ignoreNULL = T)


##################################################
##################################################
################ MOUSE-OVERS #####################
##################################################
##################################################
##################################################

show_click = reactiveVal(NULL)

print_click = reactiveVal(NULL)

output$click_info <- renderUI   (show_click())

output$point_info <- DT::renderDataTable({

  df = print_click()

 print(paste("info",rownames(df)))

  opts = list(autowidth = TRUE, scrollX = T, searching = F, paging = F,
    ordering = F, info = FALSE, columnDefs = list(list(
    width = c('10%'), targets = c(0, 1, 2) ),
    list(width = c('70%'), targets = c(3))))

  df =  DT::datatable(df, options = opts, rownames = F, escape = F,
    selection = "single") %>% formatStyle(c(1:4),fontWeight = 'bold',
    fontSize = '12px',target="row")})

 observeEvent(input$plot_click, {

   welke = which(colnames(delke)%in% input$cmpSERV)

   whop = which(colnames(d3) %in% input$cmpSERV)

   wboth = which(colnames(dn3) %in% input$cmpSERV)

   if(length(whop)>0) {

     x = d3} else if(length(wboth)>0) {

     x = dn3} else {x = delke}

   req(input$cmpSERV)

   w1 = which(colnames(x) == input$cmpSERV)

  req(length(w1) > 0)

  validate(need(length(w1) == 1 ,message = "please enter a valid compound"))

  hop <- hgeneAnno(mat = x, cmp = colnames(x)[w1], fdat = orth, arrange = FALSE)

  pclick <- nearPoints(df = hop,  xvar = "gene", yvar = "FD",

    coordinfo = input$plot_click,threshold = 10)

  hideTooltip <- function( hide ){

    session$sendCustomMessage(type = 'hideDraggable',

    message = list('hide'=hide))}

    g = grep("FD",names(pclick))

    if(length(g) > 0) pclick[,g] = format(round(pclick[,g],2),
      nsmall = 1,scientific = F)

    g = which(names(pclick) %in% c("xvar","gene"))

    if(length(g) > 0) pclick = pclick[,-g]

    if( nrow(pclick) == 0 ) {

      show_click(NULL)

    hideTooltip(TRUE) # Hide tooltip if there's no info to show

    return()

    } else {

    session$sendCustomMessage(type = 'placeDraggable', message = list())

    show_click(tagList({DT::dataTableOutput("point_info", width = "100%")}))

    print_click({pclick})}})

show_eclick = reactiveVal(NULL)

print_eclick = reactiveVal(NULL)

output$eclick_info <- renderUI   (show_eclick())

output$point_einfo <- DT::renderDataTable({

  df = print_eclick()

   opts = list(autowidth = TRUE, scrollX = T, searching = F, paging = F,
    ordering = F, info = FALSE, columnDefs = list(list(
    width = c('10%'), targets = c(0, 1, 2) ),
    list(width = c('70%'), targets = c(3))))

  df =  DT::datatable(df, options = opts, rownames = F, escape = F,

    selection = "single") %>%

 formatStyle(c(1:4),fontWeight = 'bold', fontSize = '12px',target="row")})

 observeEvent(input$plot_eclick, {

  x = de3

  req(input$cmpSERV)

  w1 = which(colnames(x) == input$cmpSERV)

  req(length(w1) > 0)

  validate(need(length(w1) == 1 ,message =  "please enter a valid compound"))

  hideTooltip <- function( hide ){
    session$sendCustomMessage(type = 'ehideDraggable', message = list('hide'=hide))
  }

  hop =  geneAnno(mat = x,cmp=colnames(x)[w1],fdat = orth,xvar=F)

    peclick <- nearPoints(df = hop,  xvar = "gene", yvar = "FD",

      coordinfo = input$plot_eclick,threshold = 5)

 g = grep("FD",names(peclick))

 if(length(g) > 0) peclick[,g] = format(round(peclick[,g],2),

   nsmall = 1,scientific = F)

 g = which(names(peclick) %in% c("xvar","gene"))

   if(length(g) > 0) peclick = peclick[,-g]

   if( nrow(peclick) == 0 ) {show_eclick(NULL)

   hideTooltip(TRUE) # Hide tooltip if there's no info to show

   return() } else {

   session$sendCustomMessage(type = 'eplaceDraggable', message = list())

   show_eclick(tagList(

      { DT::dataTableOutput("point_einfo", width = "100%")}))}})

####################################################################### #######################################################################
#######################  END ALL MOUSEOVVERS ##########################
#######################################################################
#######################################################################

####################################################################### #######################################################################
####################### START GENEBYDRUG ##############################
#######################################################################
#######################################################################

data_for_genebydrug <- reactive({

    req(input$gene)

    gene = toupper(input$gene)

    xdat = xinput()

    df = pdata()

    w = which(rownames(xdat) == gene)

    validate(need(length(w) == 1 ,message = "please enter a valid gene"))

    mx2 = mymeltdf(xdat,row = gene,df = df)

    mx2

})

##################################################
################ GENE BY SCREEN PLOT #############
##################################################

output$genebydrug <- renderPlot({

  x4 = data_for_genebydrug()

  mx2 = x4

  med = median(mx2$fitness_defect)

  tit = toupper(input$gene)

  mx2$fitness_defect = round(mx2$fitness_defect,2)

  mx2$sig = 0

  mx2$shape = 19

  g = which(mx2$fitness_defect >= 15)

  if(length(g) > 0) mx2$sig[g] = 17

  wx = which(mx2$fitness_defect >= 1)

  if(length(wx) > 0) mx2$sig[wx] = 1

  g  = ggplot(mx2,aes(x =  screen,y=fitness_defect,col = factor(sig))) + theme_bw() +
    geom_point(aes(col = factor(sig),shape = factor(shape),size = 15)) +
    scale_shape_manual(values = c(17,19))

  g1 = g + theme(legend.position="none") +
    theme(panel.grid.minor =   element_blank()) +
    theme(panel.grid.major = element_blank()) +
    theme(axis.ticks = element_blank(), axis.text.x =   element_blank())

  g1 = g1 + labs(y = "fitness defect score") + labs(x = "compound") +
    geom_hline(yintercept=median(mx2$fitness_defect),color = "black",
    linetype = "dashed",size =1) + ggtitle(tit)

  g2 = g1 + geom_text_repel(size = 6,
    data = subset(mx2, sig == TRUE),
    aes(x =  screen,y = fitness_defect, label = screen),
    point.padding = 0.25,
    segment.alpha = 0.2, col = "black")

  g2 = g2 + theme(axis.text=element_text(size=18),
    axis.title=element_text(size=24,face="bold"),
    plot.title = element_text(size = 24,face = "bold"))

    g2})

####################################################################### #######################################################################
####################### START BRUSH PTS ###############################
#######################################################################
#######################################################################

  observeEvent(input$brushpts, {

    mx2 = data_for_genebydrug()
    mx2$fitness_defect = round(mx2$fitness_defect,2)

    nam = names(mx2)

    g = grep("cmp",nam)

    pts = brushedPoints(mx2, input$brushpts, xvar = "screen", yvar = "fitness_defect")


    output$tabpts = DT::renderDataTable({

    pts


    }, escape = F,rownames = F, selection = "single",server = F)


    choices = pts$screen


  },ignoreInit = T,ignoreNULL = F)


observeEvent(input$tabpts_rows_selected, {

  row = input$tabpts_rows_selected

  req(data_for_genebydrug())

  req(input$brushpts)

  mx2 = data_for_genebydrug()

  mx2$fitness_defect = round(mx2$fitness_defect,2)

  nam = names(mx2)

  g = grep("cmp",nam)

  pts = brushedPoints(mx2, input$brushpts, xvar = "screen", yvar = "fitness_defect")

  selected = pts$screen[row]

  site = input$site

  w = which(rownames(d3) == toupper(input$gene))

  welke = which(colnames(delke)%in% selected)

  whop = which(colnames(d3) %in% selected)

  wboth = which(colnames(dint) %in% selected)

  print(c("tabpts",selected))


  if(length(whop) > 0) {site = "marjhop"}

  if(length(wboth) > 0) {site = "marjhip"}

  if(length(welke) > 0) {site = "elke"}

  req(site)

  updatePrettyRadioButtons(session,"site", label = NULL,
    c("2008 Erikson PLOS Genetics" = "elke",
      "2021 HIPHOP Marjan" = "marjhip",
      "2020:2024 HOP Marjan" = "marjhop"), selected = site)

  if(selected%in% phop$name) {
    choices = phop$name} else if(selected%in% pmarjhip$name) {
    choices = pmarjhip$name}else if(selected%in% pelke$name) {
    choices = pelke$name}

  updateSelectizeInput(session,"cmpSERV",label = NULL,
      choices = choices, selected = selected,server = T)

  newtab <- switch(input$tabs,
  "genebyscreen" = "hiphop",
  "hiphop" = "genebyscreen")

  updateTabItems(session, "tabs", newtab)

})


##################################################
################ GENE BY SCREEN ##################
#################################################


####################################################
################### START BRUSH ####################
####################################################

####################################################
################### END BRUSH ######################
####################################################

output$geneinfo  <- renderDataTable({

 w = which(orth$best_human_gene %in% rownames(dint))

 fd = orth[w,] %>% distinct(best_human_gene,.keep_all = T)

  w1 = which(fd$best_human_gene %in% toupper(input$gene))

  validate(
    need(length(w1) == 1 ,message =
           "please enter a valid gene"))

  f = fd[w1,c("best_human_gene","genecards_link","descriptor")]

  names(f)[1] = c("GENE")

  f

},escape = F,options=list(paging = F,searching = F,

  info = F,scrollX = T,scrollY=TRUE,scrollCollapse=TRUE,autowidth=TRUE),

  rownames = F)

########################################################################
############################ END GENEBYSCREEN ##########################
########################################################################

##################################################
########### START GO ENRICHMENT BY RANK ##########
##################################################

output$fiToutput = renderUI({


  xfit = outFit$dFit
  xfit$cofit = as.numeric(xfit$cofit)
  xfit$pvalue = as.numeric(xfit$pvalue)

  mfit = as.matrix(xfit[,2:3])

  mfit[,2] = 1:nrow(mfit)
  colnames(mfit)[2]="rank"


  tagList(
    fluidRow(
      box(title = "Set cofitness rank threshold:",
      HTML("<h5><b>Top: 1000 genes ranked by cofitness, from 1 to 1000, left to right;
           # on scale marks # of strains used in GO enrichment.
           Bottom: Corresponding cofitness values 1 to 1000, right to left.
           Click on arrow to animate GO enrichment
           dependency on cofitness rank.</b></h5>"),
    br(),

      br(),
      HTML("<h5><b>Cofitness values:</b></h5>"),

    st =rev(seq(0,1,0.05)),

    sq = rev(seq(0,1,0.001)),

    sliderTextInput(inputId = "corrFIT", label = NULL,

    choices = st,grid = T,animate = T),

      HTML("<h5><b>Corresponding cofitness rank:</b></h5>"),

      sliderInput("scoreFIT",label = "",
                value = 1, min = 1,max = 1000, step = 1,animate = TRUE),

      HTML("<h5><b>Set FDR threshold for GO enrichment by RANK:</b></h5>"),

          sliderInput( "fdrSERV", label = NULL,
            min = 0, max = 0.5,value = 0.1, step = 0.05),


    status = "primary", solidHeader = T, width="100%", height = 650)))
})


rankFIT = debounce(reactive(input$scoreFIT),500)
coFIT = debounce(reactive(input$corrFIT),500)

scoreVal = reactiveValues(
  previous = NULL,
  current  = NULL
)

corrVal = reactiveValues(
  previous = NULL,
  current  = NULL
)



# 0,49 621 0.28 721 rank 0.74 271 0.56 441,[ 0,908 93 ]0.86 93 102 0.899 102 0.867
   observeEvent(input$corrFIT,{

     req(input$scoreFIT)

     req(input$corrFIT)

     print(c("input$scoreFIT",input$scoreFIT))
     print(c("input$corrFIT",input$corrFIT))

     corr = input$corrFIT

     mfit = outMit$mFit

     sss = rev(seq(0, 1, length.out = nrow(mfit)))

     ss = round(sss, 4)

     w = which(mfit[,3] == corr)

     req(length(w) > 0)

     score = mfit[w,"rank"]

     w1 = slid(x =input$corrFIT,y = input$scoreFIT, mfit = mfit)

     score1 = w1$corr[2]

     print(c("obfit","score1",score1,"score",score, "wscore",
             w1$score,"current" ,corrVal$current,
             "corr",w1$corr,"scoreFIT",input$scoreFIT,"corrFITr",
             input$corrFIT,"previp" ,corrVal$previous,"currs" ,
             scoreVal$current,"current" ,scoreVal$previous))

     test = c(as.numeric(input$scoreFIT) != scoreVal$current)

     test2 = as.numeric(input$scoreFIT) != score

     print(c("test",test))

     print(c("test2",test2))

     print(c("scorefit",score,input$scoreFIT,scoreVal$current))

     cond =   is.null(scoreVal$current)

  if( is.null(scoreVal$current)){

     updateSliderInput(session = session,
          inputId = 'scoreFIT',
          value = score) # updateSliderInput

   } else if(!is.null(scoreVal$current) & as.numeric(input$scoreFIT) != score){

     updateSliderInput(session = session,
          inputId = 'scoreFIT',
          value = score)}

     scoreVal$previous = scoreVal$current

     scoreVal$current = score

    }, ignoreNULL = F, ignoreInit = T)
#
#
#
    observeEvent(input$scoreFIT,{

      mfit = outMit$mFit

      req(input$scoreFIT)
      req(input$corrFIT)
#
# print(c("input$scoreFITRANK",input$scoreFIT))
#  print(c("input$corrFITRANK",input$corrFIT))
      w1 = slid(x =input$corrFIT,y= input$scoreFIT, mfit = mfit)
      score1 = w1$score[1]


      rank=input$scoreFIT
      w = which(mfit[,2] == rank)
      req(length(w)>0)
      score = mfit[w,3]

     print(c("obscore","score1",score1,"score",score, "wscore",w1$score,corrVal$current,
             "corr",w1$corr,"scoreFIT",input$scoreFIT,"corrFITr",
             input$corrFIT,"previcorrscore" ,corrVal$previous,"currscorecurr" ,
             scoreVal$current,"previcorrscore",scoreVal$previous))

      print(c("scorerank",score))
      cest = as.numeric(input$corrFIT) != score

      cest2 = as.numeric(input$corrFIT) != corrVal$current

      req(cest ==TRUE |cest ==FALSE )

      print(c("cest",cest))

      print(c("cest2",cest2))

      print(c("rankfit",score,input$corrFIT,corrVal$current))

   sq = rev(seq(0,1,0.001))


   sss = rev(seq(0, 1, length.out = nrow(mfit)))

   ss = round(sss, 4)


   if( is.null(corrVal$current)){

      updateSliderTextInput(session,"corrFIT",label =  NULL,
      choices = ss, selected = score)

    } else if(as.numeric(input$corrFIT) != score){

      updateSliderTextInput(session,"corrFIT",label =  NULL,
      choices = ss, selected = score)}

      corrVal$previous = corrVal$current

      corrVal$current  = score

  }, ignoreNULL = F, ignoreInit = F)

### this is for fdrRANK

output$mxFIT <- renderUI({
  tagList(

    box(title = "Optimizing GO enrichment:",
        tags$head(tags$style(HTML("<h5><b>label(font-family:Arial;
                                  font-weight:bold;font-size;20)<h5><b>"))),

        HTML("<h5><b>Finds the cofitness/cofitness rank that maximizes
             the number of GO enrichment terms for the top 1000 cofit
             strains for a given FDR value. Stepsize is the increment
             between ranks tested.<h5><b>"),

        HTML("<h5><b>Choose a stepsize and FDR threshold, then hit GO to run.<h5><b>"),

        HTML("<h5><b>Warning: calculation time for stepsize 1: ~ 3min.<h5><b>"),
        br(),
        fluidRow(

          column(width = 12,align = "center",

                 prettyRadioButtons(inputId = "stepsize",label = "stepsize",
                                    choices = c(1,5,20,50), outline = T, fill = F,thick = T,
                                    status = "primary", shape = "round", bigger = T,
                                    width = "100%",
                                    selected = 50, inline = T)),

          HTML("<h5><b>&nbsp&nbsp&nbsp&nbsp&nbspSet FDR threshold</b></h5>"),

          br(),
          column(width = 12,align = "center",

                 sliderInput("fdrRANK",label = "", min = 0, max = 0.5,
                             value = 0.2, step = 0.1)),


          column(width = 12,align = "left",
               actionButton("rankGO", "GO!", class = "btn-success")),

          column(width = 12, align = "left",
               tags$head(tags$style(HTML("<h5><b>label(font-family:Arial;
                                         font-weight:bold;font-size;20)<h5><b>"))),
               htmlOutput("finishTEXT"),
               htmlOutput("stepTEXT"),

               withSpinner(htmlOutput("mxCOFITvalue"),type = 5,
                           color = "#4d88ff", hide.ui = F, size = 0.5))),

        status = "primary", solidHeader = T,width="100%", height = 650))

})


gseaMAX <- reactiveValues(
  rank     = NULL,
  cofit    = NULL,
  goterms     = NULL)

observeEvent(outFit$dFit,{
  useShinyjs()
  click("rankGO")
  reset("scoreFIT")
  reset("fdrRANK")
},ignoreNULL = T,ignoreInit = F, once = T)

observeEvent(input$gene,{
  useShinyjs()
  reset("scoreFIT")
  reset("fdrRANK")

  shinyjs::hide(id = "finishTEXT")
  shinyjs::hide(id = "stepTEXT")
  shinyjs::hide(id = "mxCOFITvalue")

},ignoreNULL = T,ignoreInit = T)


observeEvent(input$fdrRANK,{

  useShinyjs()
  reset("scoreFIT")
  shinyjs::hide(id = "finishTEXT")
  shinyjs::hide(id = "stepTEXT")
  shinyjs::hide(id = "mxCOFITvalue")
  reset("corrFIT")
},ignoreNULL = T, ignoreInit = T)

observeEvent(input$stepsize,{
  useShinyjs()
  #click("rankGO")
  reset("scoreFIT")
  reset("fdrRANK")
  shinyjs::hide(id = "finishTEXT")
  shinyjs::hide(id = "stepTEXT")
  shinyjs::hide(id = "mxCOFITvalue")
},ignoreNULL = T, ignoreInit = T, once = F)


observeEvent(input$rankGO,{

  req(outMit$mFit)

  req(input$fdrSERV)

  gseaMAX$rank = NULL
  gseaMAX$cofit = NULL
  gseaMAX$goterms = NULL
  shinyjs::hide(id = "finishTEXT")
  shinyjs::hide(id = "stepTEXT")
  shinyjs::hide(id = "mxCOFITvalue")

  mfit = outMit$mFit

  FDR = input$fdrRANK

  mn = 1

  mx = 1000

  steps = pretty(c(mn,mx),n = 100)

  range = max(steps)-min(steps)

  step = range/(length(steps) -1)


  if(as.numeric(input$stepsize) == 1) {

    steps = 1:1000} else if(as.numeric(input$stepsize) == 5){

    steps = seq(0,1000,5)} else if(as.numeric(input$stepsize) == 20){

    steps = seq(0,1000,20)} else if(as.numeric(input$stepsize) == 50){

    steps = seq(0,1000,50)}

  dFIT = NULL

  for(i in 1:length(steps)) dFIT[[i]] = compRANK(mfit,coln = 2, rank  = steps[i])

  names(dFIT)=steps

  net = lapply(dFIT,runGOENRICH,fdrThresh = FDR, minSetSize = 5,
    curr_exp = "tst",bp_path =  "2023_January30_GOBP.RDS",
    go_path = "2023_January30_GOID_GOBP_SGD.txt")

  enrich = lapply(net,function(x) x = x[[1]])

  isnull = sapply(enrich,is.null)

  wisnull = which(isnull ==T)

  nrw = sapply(enrich,nrow)

  if(length(wisnull) > 0) isnull = lapply(nrw[wisnull], function(x) x = 0)

  nrw[wisnull] = isnull

  unrw = unlist(nrw)

  mx = which.max(unrw)

  nmx = as.numeric(names(mx))

  cofit = mfit[nmx,1]

  goterms = unrw[mx]

  if(length(nmx)>0) gseaMAX$rank = nmx

  if(length(nmx)>0) gseaMAX$cofit = cofit

  if(!is.null(unrw[mx])) gseaMAX$goterms = goterms

  shinyjs::show(id = "finishTEXT")

  shinyjs::show(id = "stepTEXT")

  shinyjs::show(id = "mxCOFITvalue")

  req(!is.null(gseaMAX$cofit))

  output$finishTEXT = renderText({

    paste("<h5><b>","Results:","</b></h5>")})

  output$stepTEXT = renderText({

    paste("<h5><b>","STEPSIZE = ",input$stepsize,";","FDR = ",input$fdrSERV,"</b></h5>")})

  output$mxCOFITvalue = renderText({

   req(!is.null(gseaMAX$cofit))

   paste("<h5><b>","COFITNESS = ",round(as.numeric(isolate(gseaMAX$cofit)),3),
      ";","RANK = ", isolate(gseaMAX$rank),"</b></h5>")})

},ignoreNULL = T,ignoreInit = T)


output$mxCOFITvalue = renderText({

  NULL

})

output$GOterms <- renderUI({

  htmlOutput("countGO")

})

observeEvent(gseaMAX$goterms,{

  output$countGO = renderText({paste("<h5><b>","no. GO terms =",
    gseaMAX$goterms,"</b></h5>")})

})

observeEvent(gseaMAX$rank,{

  req(outMit$mFit)

  req(!is.null(gseaMAX$rank))

  mfit = outMit$mFit

  mx = nrow(mfit)

  value = gseaMAX$rank

  scoreVal$previous = input$scoreFIT

  scoreVal$current = gseaMAX$rank

  updateSliderInput(session, "scoreFIT",label = NULL,min = 1,
    max = 1000, value = gseaMAX$rank, step = 1)

},ignoreNULL = T, ignoreInit = T)

observeEvent(rankFIT(),{

  req(outMit$mFit)

  req(!is.null(input$corrFIT))

  if(is.null(scoreVal$current)) scoreVal$current = rankFIT()

  if(rankFIT()!=scoreVal$current) {

    scoreVal$previous = scoreVal$current

    scoreVal$current  = rankFIT()}

  mfit = outMit$mFit

  mx = nrow(mfit)

  w = which(mfit[,2]==rankFIT())

  req(length(w)!=0)

  corel = mfit[w,1]

  if(length(corel)>1) corel = corel[1]

    updateSliderInput(session,"corrFIT",label = NULL,
        min = 0, max = 1.0,  value = corel, step = 0.1)

  corrVal$previous = corrVal$current

  corrVal$current = corel#}

}, ignoreNULL = T,ignoreInit = F)

###########################################################
###################### END COFITNESS ######################
###########################################################

##################################################
########### START GO ENRICHMENT BY RANK ##########
##################################################

goEnrich <- reactive({

     req(input$fdrSERV)

     req(outMit$mFit)

     req(input$scoreFIT)

     mfit = outMit$mFit

     gseaMAX$goterms = NULL

     num = is.numeric(input$scoreFIT)

     req(num == TRUE)

     wfit = which(mfit[,2] <= input$scoreFIT)

     validate(need(length(wfit)!=0, message = "No CoFitness scores above threshold"))

     req(length(wfit)!=0)

     qFIT = compRANK(mfit,coln = 2, rank  = input$scoreFIT)

     curr_exp = "cofit"

     FDR = input$fdrSERV

     network = runGOENRICH(fdrThresh = FDR, minSetSize = 5,

        curr_exp = curr_exp,score = qFIT, bp_path =  "2023_January30_GOBP.RDS",

        go_path =    "2023_January30_GOID_GOBP_SGD.txt")


     validate(need(!is.null(network$enrichInfo), message = "No GO enrichment,
           try relaxing the FDR or scorethreshold"))

    enrichInfo = network$enrichInfo

    req(!(is.null(enrichInfo)))

    if(!is.null(enrichInfo)) gseaMAX$goterms = nrow(enrichInfo)

    edgeMat = network$edgeMat

    if(!is.null(enrichInfo)) return(network)
  })

##################################################
############# END GO ENRICHMENT ##################
##################################################
visNet <- reactive({

    req(!is.null(goEnrich()$enrichInfo))

    enrich = goEnrich()$enrichInfo

    edge = goEnrich()$edgeMat

    vis = visSetup(enrichInfo = enrich,edgeMat = edge, fontsize = 20, fontface = "Courier")


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

             style = 'width: 500px; height = 31px; font-size: 18px;

             color: #000066;border: 3px solid #4d88ff;')) %>%

          visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);;}")} else {

           visNetwork(vis$nodes, vis$edges, width = "100%") %>%

           visNodes(shadow=list(enabled=T,size=25)) %>%

           visOptions(highlightNearest = list(enabled = T, degree = 5, hover = T),

          nodesIdSelection = list(enabled = TRUE, values = names,

          style = 'width: 500px; height = 31px; font-size: 18px;

          color: #000066;border: 3px solid #4d88ff;'),

          selectedBy = list(variable="FDR",

            style = 'width: 500px; height = 31px; font-size: 18px;

            color: #000066;border: 3px solid #4d88ff;')) %>%

        visIgraphLayout(type = "full") %>%

        visEvents(select = "function(nodes) {

                Shiny.onInputChange('current_node_id', nodes.nodes);;}")}

  })

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

          style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;')) %>%

        visEvents(select = "function(nodes) {

          Shiny.onInputChange('current_node_id', nodes.nodes);;}")} else {

          visNetwork(vis$nodes, vis$edges, width = "100%") %>%

          visNodes(shadow=list(enabled=T,size=25)) %>%

          visOptions(

            highlightNearest = list(enabled = T, degree = 5, hover = T),

            nodesIdSelection = list(enabled = TRUE, values = names,

            style = 'width: 500px; height = 31px; font-size: 18px;

            color: #000066;border: 3px solid #4d88ff;'),

            selectedBy = list(variable="FDR",

            style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;')) %>%

        visIgraphLayout(type = "full") %>%

        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);;}")}})

observeEvent(input$network_Proxy_selectedBy, {

  req(input$network_Proxy_selectedBy)

  n = visNet()$nodes

  w = which(n$FDR %in% as.numeric(input$network_Proxy_selectedBy))

  id = n$id[w]


  visNetworkProxy("network_Proxy") %>% visSelectNodes(id = id)

}, ignoreNULL = T, ignoreInit = F)

##################################################
################### START PROXY ##################
##################################################

######################################################################
######################### START BARPLOTS #############################
######################################################################

output$GoTable = renderUI({

  req(input$network_Proxy_selected)

  DT::dataTableOutput("gotermTable")

  })

output$gotermTable <- DT::renderDataTable({

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

    datatable(t,width=220,caption = htmltools::tags$caption(term,

    style = "caption-side: top; text-align: center; color:black;background:white;

    font-weight:bold;"),

    options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F,

      rowCallback = JS("function(row, data) {", "for (i = 1; i < data.length; i++) {",

      "if (data[i]>1000 | data[i]<1){", "$('td:eq('+i+')', row).html(data[i].toExponential(1));", "}", "}", "}")),

    height = 400,colnames = "") %>%

      formatStyle( target = "row", color = "black",backgroundColor = "white",

    columns = c(1,2),fontWeight = "bold")

  })

output$leadingEdge = renderUI({
    req(input$network_Proxy_selected)

    plotOutput("BarPlot",
          width = input$shiny_width/300,
          height = input$shiny_height/300)})

output$BarPlot <- renderPlot({

  req(input$network_Proxy_selected)

  vis = visNet()

  n = visNet()$nodes

  w = which(vis$nodes$id %in% c(input$network_Proxy_selected))

  req(length(w) > 0)

  n = vis$nodes[w,]

  req(nrow(n)!=0)

  s6 = geneBARPLOT(n$overlapGenes)

  d = s6

  if(nrow(d)> 10) d = d[1:10,]

  barplot(d$score, names.arg = d$gene,las=1, horiz=T,cex.axis = 0.8,

    cex.lab = 1.1, col="navy",xlab = "median fitness defect score")})

######################################################################
######################### END BARPLOTS ###############################
######################################################################

#################################################
########## START GO ENCHICHMENT TABLES ##########
#################################################

EnrichReact = reactive({

    req(!is.null(visNet()$nodes))

    enrich = visNet()$nodes

    req(length(nrow(enrich)) > 0)

    row <- input$EnrichTable_rows_selected

    out = OutEnrich()

    id = out$id[row]

    id})

observeEvent(EnrichReact(), {

  req(EnrichReact())

  id = EnrichReact()

  visNetworkProxy("network_Proxy") %>%

  visSelectNodes(id = id)

  }, ignoreNULL = T, ignoreInit = F)

OutEnrich = reactive({

  req(goEnrich()$enrichInfo)

  enrich = goEnrich()$enrichInfo

  enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")] =

    format(round(enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")],2),

    nsmall = 1,scientific = F)

  enrich[,c("P", "FDR")] =

    format(signif(enrich[,c("P", "FDR")],2),nsmall = 1,scientific = T)

    w = which(names(enrich) %in% c("formattedLabel","cluster","size",

      "filename","maxOverlapGeneScore"))

    enrich = enrich[,-w]

    enrich = enrich[,c("id",  "GOID","term","nGenes","nQuery",

     "nOverlap", "querySetFraction","geneSetFraction",

     "foldEnrichment","P","FDR","overlapGenes")]

    enrich})

###########################  GO ENRICHMENT TABLE  ###################################

output$EnrichTable = DT::renderDataTable({

  out = OutEnrich()

  out = out[,c("GOID","term","querySetFraction", "geneSetFraction",

    "foldEnrichment" , "FDR","overlapGenes" )]

  opts = list(pageLength = 5, autoWidth = F,scrollX = T)#,

  df =  DT::datatable(out, options = opts, rownames = F, escape = F, selection = "single") %>%

    formatStyle(c(1:7),fontWeight = 'bold')})

output$enrichdownload <- downloadHandler(

  filename = function() {

    paste0("enrich:",input$gene,"_", Sys.Date(), ".txt")
  },

  content = function(file) {

    write.table(OutEnrich(), file, row.names = F,sep="\t")})

#################################################
########## END GO ENCHICHMENT TABLES ############
#################################################

############################################################
############################################################
############# START HIPHOPTABLES AND PLOTS #################
############################################################
############################################################

output$hiptab = DT::renderDataTable({

  req(input$cmpSERV)

  mat = de3

  hip = geneAnno(mat,fdat = orth,cmp = input$cmpSERV, xvar=F)

  hip = hip[,c("gene","GENE","FD","descriptor")]

  hip$FD = format(round(hip$FD,2),nsmall = 1,scientific = F)

  hip

    },escape = F, options=list(pageLength=10, autoWidth = T,scrollX = TRUE,

      columnDefs = list(list(width = "80px", targets = c(-2,-3)))),rownames = F)

    outhip = reactive({

    req(input$cmpSERV)

    hip = geneAnno(mat = de3,fdat = orth,cmp = input$cmpSERV, xvar=F)

    hip = hip[,c("gene","GENE","FD","descriptor")]

    hip

  })

observeEvent(input$hiptab_rows_selected, {

  row <- input$hiptab_rows_selected

  df = outhip()

  updateTextInput(session, inputId = "gene", label = NULL, value = df$gene[row])

  newtab <- switch(input$tabs,

  "hiphop" = "genebyscreen",

  "genebyscreen" = "hiphop")

  updateTabItems(session, "tabs", newtab)})

output$downloadhip <- downloadHandler(

  filename = function() {

   paste0("HIP:",input$cmpSERV,  "_",Sys.Date(), ".txt")},

   content = function(file) {

  write.table(as.data.frame(outhip()), file, row.names = F,sep="\t",quote=F)})


output$efd  = renderPlot({

  req(input$cmpSERV)
  validate(
    need(input$cmpSERV,message =
           "please select condition"))
  cond1 = input$cmpSERV

  x = de3
  wx = which(colnames(x) %in% cond1)

  req(length(wx)!=0)
    req(input$site == "marjhip")
  colnames(x)[wx] = paste0("HIP|",colnames(x)[wx])
  coln = colnames(x)[wx]
  # coln   = str_wrap(coln,width=30,whitespace_only = F)
  colnames(x)[wx] = coln
  p10(x,wx,pch = 17,ylim = input$elim,sig = input$scorethreshSERV)
    })

essPlot = function(){
  validate(
    need(input$cmpSERV,message =
           "please select condition"))
  cond1 = input$cmpSERV
  if(input$site == "marjhip") {x = de3}

  wx = which(colnames(x) %in% cond1)

  colnames(x)[wx] = paste0("HIP|",colnames(x)[wx])



  p10(x,wx,pch = 17)


}
essPlot = function(){

  validate(need(input$cmpSERV,message = "please select condition"))

  cond1 = input$cmpSERV

  if(input$site == "marjhip") {x = de3}

  wx = which(colnames(x) %in% cond1)

  colnames(x)[wx] = paste0("HIP|",colnames(x)[wx])

  p10(x,wx,pch = 17)}

##########################

output$ExportessPlot <- downloadHandler(

  filename = function() {

  paste0("HIPplot:",input$cmpSERV,  "_",Sys.Date(), ".png")},

  content = function(file){png(file,

  width = (input$shiny_width/12)*6, height = 1000)

  essPlot()

    dev.off()})

output$hoptab = DT::renderDataTable({

  req(input$cmpSERV)

  if(input$site == "marjhop"){

    mat = d3[,phop$name]} else if(input$site == "marjhip"){

    mat = dn3} else {mat = delke}

    hop =  geneAnno(mat = mat,cmp = input$cmpSERV,fdat = orth, xvar = F)

    hop = hop[,c("gene","GENE","FD","descriptor")]

    hop$FD = format(round(hop$FD,2),nsmall = 1,scientific = F)

    hop

    },escape = F, options=list(pageLength=5, autoWidth = FALSE,scrollX = TRUE,

      columnDefs = list(list(width = "80px", targets = c(-2,-3)))),rownames = F)

outhop = reactive({

  req(input$cmpSERV)

  if(input$site == "marjhop"){

    mat = d3[,phop$name]} else if(input$site == "marjhip"){

    mat = dn3} else {mat = delke}

    hop =  geneAnno(mat = mat,cmp = input$cmpSERV,fdat = orth, xvar = F)

    hop = hop[,c("gene","GENE","FD","descriptor")]

    hop$FD = format(round(hop$FD,2),nsmall = 1,scientific = F)

    hop})

observeEvent(input$hoptab_rows_selected, {

  row <- input$hoptab_rows_selected

  df = outhop()

  updateTextInput(session, inputId = "gene", label = NULL, value = df$gene[row])})



output$downloadhop <- downloadHandler(

  filename = function() {

   paste0("HOP:",input$cmpSERV,  "_",Sys.Date(), ".txt")},

   content = function(file) {

  write.table(as.data.frame(outhip()), file, row.names = F,sep="\t",quote=F)})


##################################################
output$nfd = renderPlot({
  req(input$site)
  req(input$cmpSERV)
  validate(
    need(input$cmpSERV,message =
           "please select condition"))

  cond1 = input$cmpSERV

  if(input$site == "marjhip") {x = dn3}
  if(input$site == "marjhop") {x = d3}
  if(input$site == "elke") {x = delke}
  wx = which(colnames(x) %in% cond1)
  #print("wx",length(wx))
  colnames(x)[wx] = paste0("HOP|",colnames(x)[wx])

# if(input$site == "marjhip"){
#   # coln = colnames(x)[wx]
#   # coln   = str_wrap(coln,width=30,whitespace_only = F)
#   colnames(x)[wx] = coln
#    }
  req(length(wx) > 0)
  p10(x,wx,pch = 17, ylim = input$ylim,sig = input$scorethreshSERV)


})

##########################
nonPlot = function(){

  validate(need(input$cmpSERV,message = "please select condition"))

  cond1 = input$cmpSERV

  if(input$site == "marjhop"){

    x = d3[,phop$name]} else if(input$site == "marjhip"){

    x = dn3} else {x = delke}

    wx = which(colnames(x) %in% cond1)

    colnames(x)[wx] = paste0("HOP|",colnames(x)[wx])

    p10(x,wx,pch = 17)}

output$ExportnonPlot <- downloadHandler(

  filename = function() {

  paste0("HOPplot:",input$cmpSERV,  "_",Sys.Date(), ".png")},

  content = function(file){

  png(file, width = (input$shiny_width/12)*6, height = 1000)

  nonPlot()})

############################################################
############################################################
############# END HIPHOPTABLES AND PLOTS ###################
############################################################
############################################################

#######################################################
#######################################################
################## START Y-AXIS LIMITS ################
#######################################################
#######################################################

output$ylimits = renderUI({

   x = xinput()

   req(input$cmpSERV)

   coln = colnames(x)

   wcoln = which(coln %in% input$cmpSERV)
   req(length(wcoln) > 0)

   req(xinput())
   req(input$cmpSERV)

   x = xinput()

   wne = which(rownames(x) %in% noness$best_human_gene)

   we = which(rownames(x) %in% ess$best_human_gene)

   print(c("wewne", length(we) , length(wne)))

   if(length(wne) > 0) {

      ylims = c(floor(min(x[wne,wcoln,drop=F],na.rm = T)),
                ceiling(max(x[wne,wcoln,drop=F],na.rm = T)))
      }
     req(length(wne) > 0)

     if(length(we) > 28) {

     elims = c(floor(min(x[we,wcoln,drop=F],na.rm = T)),
               ceiling(max(x[we,wcoln,drop=F],na.rm = T)))
     }




tagList(



  if(length(we) > 28){
    column(width = 3,
  box(title= "HIP y-axis range:",

    fluidRow(
      column(width = 12,
      sliderInput("elim", label = NULL, min = elims[1],
        max = elims[2], value = c(-0.5,elims[2])))),

    fluidRow(
        column(width = 12,align = "center",
        actionLink(
          inputId = "eorig",
          label = "reset",
          style = "height: 100px;width: 150px;font-size:120%;text-align:center",
          size = "lg"))),
        status = "primary", solidHeader = T,width = "100%", height = 200) )

},


   if(length(ylims) == 2){
   column(width = 3,

    #column(1,
  box(title= "HOP y-axis range:",
       fluidRow(
         column(width=12,

       # Copy the line below to make a slider range
         sliderInput("ylim", label = NULL, min = ylims[1],
        max = ylims[2], value = c(-0.5,ylims[2])))),

        fluidRow(
        column(width = 12,align = "center",
        actionLink(
          inputId = "yorig",
          label = "reset",
          style = "height: 100px;width: 150px;font-size:120%;text-align:center",
          size = "lg"))),


         status = "primary", solidHeader = T,width = "100%", height = 200))

}
    )

})

################# Y AXES LIMITS ###################
###################################################
###################################################
 observeEvent(input$yorig,{

   x = xinput()
   coln = colnames(xinput())

   wcoln = which(coln %in% input$cmpSERV)

   wne = which(rownames(x) %in% noness$best_human_gene)

   if(length(wne) > 0) {
      xne = xinput()
      ylims = c(floor(min(xne[wne,wcoln,drop=F],na.rm = T)),
                ceiling(max(xne[wne,wcoln,drop=F],na.rm = T)))
      }

   updateNumericRangeInput(session, "ylim", "", value = ylims)

 },ignoreInit = F, ignoreNULL = T)


observeEvent(input$eorig,{

   cole = colnames(xinput())
   xe = xinput()
   wcole = which(cole %in% input$cmpSERV)

   we = which(rownames(xe)%in% ess$best_human_gene)

   if(length(we) > 28) {elims = c(floor(min(xe[we,wcole,drop=F],na.rm = T)),
                                 ceiling(max(xe[we,wcole,drop=F],na.rm = T)))
     }
   updateNumericRangeInput(session, "elim", "",value = elims)

})


################# Y AXES LIMITS ###################
###################################################
###################################################


observeEvent(input$resetscore, {
  updateNumericInput(session, "scorethreshSERV", "", value = 1.0, step =0.1, min = 0, max = 5.0)
},ignoreInit = F, ignoreNULL = T)



#######################################################
#######################################################
################## END Y-AXIS LIMITS ##################
#######################################################
#######################################################


#######################################################
#######################################################
################## START COINHIBITION #################
#######################################################
#######################################################

outcoInhib = reactive({

  req(input$cmpSERV)

  if(input$cmpSERV%in% phop$name) {

            xcov = xcoinh
            xpv =  xpval
            pdata = phop

      } else if(input$cmpSERV%in% pmarjhip$name){

            xcov = xcohh
            xpv = xpvhh
            pdata = pmarjhip

        } else {

            xcov = xelke
            xpv = xelkepv
            pdata = pelke
        }

      w = which(colnames(xcov) == input$cmpSERV)

      validate(need(length(w) == 1 ,message = "please enter a valid condition"))

      p = xpv[,w,drop=F]

      d = xcov[,w,drop = F]

      d = d[order(d[,1],decreasing = T),,drop=F]

      df = data.frame(screen = rownames(d)[1:nrow(d)],

      CoI = d[1:nrow(d),1],pvalue = p[rownames(d),1],stringsAsFactors = F)

      df$pvalue = formatC(df$pvalue,format = "e", digits = 2)

      w = which(df$CoI >= 0)

      if(length(w) > 0) df = df[w,]

      m = match(df$screen,pdata$name)

      df$CoI = format(round(df$CoI,2),nsmall = 1,scientific = F)

      df$CID = pdata$PCID[m]

      df$structure = pdata$image[m]

      w = which(names(df) == "CID")

      names(df)[w] = '<a href="https://pubchem.ncbi.nlm.nih.gov" target="_blank">CID</a>'

      df$mechanism = pdata$mechanism[m]

      wna = which(is.na(df$structure))


      d = df[-wna,]


      df

   })

####################################################
####################################################
####################################################
####################################################

output$coFitNess = DT::renderDataTable({

      req(outFit$dFit)

      outFit$dFit

    },escape = F,rownames=F,options=list(pageLength=10),selection = "single",server = F)

####################################################

observeEvent(input$coFitNess_rows_selected, {

      row = input$coFitNess_rows_selected

      req(outFit$dFit)

      xfit = outFit$dFit

      w = xfit$gene[row]

      validate(need(length(w) == 1 ,message = "please enter a valid gene"))

      updateTextInput(session, inputId = "gene", label = NULL, value = xfit$gene[row])

    })


output$coInhib = renderDataTable({
    xcoinhib = outcoInhib()

    options = list(
        dom = 'Bfrtip',
        pageLength = 1000,
        autoWidth = FALSE,
        scrollX = TRUE)

     DT::datatable( xcoinhib, options = options,
   escape = F,  rownames = F, selection = "single")

  })

observeEvent(input$coInhib_rows_selected, {

  row <- input$coInhib_rows_selected

  df = outcoInhib()

  print(c(nrow(df),"outcoInhib"))

  updateSelectizeInput(session,"cmpSERV",label = NULL,
  choices = df$screen, selected = df$screen[row][1],server = T)})


  output$downloadcoinhib <- downloadHandler(
    filename = function() {
      paste0("HOP:",
         input$cmpSERV,
             "_",Sys.Date(), ".txt")
    },
    content = function(file) {
      write.table(as.data.frame(outcoInhib()), file, row.names = F,sep="\t",quote=F)})

##################################################
################# END COFITNESS ##################
##################################################

##################################################
############### START outFit outMit ##############
##################################################

outFit = reactiveValues(dFit = NULL)

outMit = reactiveValues(mFit = NULL)

observeEvent(input$gene,{

  req(input$gene)

  req(input$site)

  if(input$site == "marjhip") {xfit = xfithh}
  if(input$site == "elke") {xfit = xefit}
  if(input$site == "marjhop") {xfit = xcofit}

  if(input$site == "marjhip") {xpv = xpithh}
  if(input$site == "elke") {xpv = xefitpv}
  if(input$site == "marjhop") {xpv = xpfit}

  ww = which(colnames(xfit) == toupper(input$gene))

  validate(need(length(ww) == 1 ,message = "please enter a valid gene"))

  if(length(ww) != 0) {
        p = xpv[,ww,drop=F]
        d = xfit[,ww,drop = F]
       }

  d = d[order(d[,1],decreasing = T),,drop=F]

  df = data.frame(gene = rownames(d), cofit = d[rownames(d),1],
    pvalue = p[rownames(d),1],stringsAsFactors = F)

  w = which(df$cofit >= 0)

  if(length(w) >0) df = df[w,]

  xfit = df

  xfit$cofit = as.numeric(xfit$cofit)

  xfit$pvalue = as.numeric(xfit$pvalue)

  xfit$pvalue = formatC(xfit$pvalue,format = "e", digits = 2)

  xfit$cofit = format(round(xfit$cofit, 6), nsmall = 2)

  mfit = as.matrix(xfit[,2:3])

  mfit[,2] = 1:nrow(mfit)

  colnames(mfit)[2]="rank"

  m = mfit

  sss = rev(seq(0, 1, length.out = nrow(mfit)))

  ss = round(sss, 4)

  m = cbind(m,ss)

  mm = apply(m,2,as.numeric)

  rownames(mm) = rownames(m)

  outFit$dFit = xfit

  outMit$mFit = mm

})

##################################################
############### END outFit outMit ################
##################################################

 }
  ##########################
shinyApp(ui, server)
