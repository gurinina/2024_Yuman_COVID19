
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
library(shinyjs)
library(shinyBS)
library(stringr)
library(shinyjqui)

options(rsconnect.max.bundle.size = 3145728000)

source(file.path("modules", "2022_March13_GOENRICH.R"))
source(file.path("modules" ,"2024_June19_hgeneAnnofunctions.R"))
source(file.path("modules" ,"2024_March22_functions.R"))
source(file.path("modules" ,"2024_July13_signature.R"))
source(file.path("modules" ,"2024_April11_visNetwork.R"))

xsig        <-      readRDS("24Jan8/x21.RDS")
orth        <-      readRDS("24Jan8/orthologs.RDS")
   w        <-      which(orth$hessential=="ess")
noness      <-      orth[-w,]
ess         <-      orth[w,]

phop        <-      readRDS("24Jan8/phop.RDS")
pelke       <-      readRDS("24Jan8/pelke.RDS")
pmarjhip    <-      readRDS("24Jan8/phip.RDS")

d3          <-      readRDS("24Jan8/d4.RDS")
dn3         <-      readRDS("24Jan8/dn4.RDS")
de3         <-      readRDS("24Jan8/de4.RDS")
delke       <-      readRDS("24Jan8/delke.RDS")
dint        <-      readRDS("24Jan8/dint.RDS")

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

##########################

header <-  dashboardHeader(title =  span("Yuman Chemogenomic Profiling  COVID19",
             style = "font-weight:bold;font-size:18px;
             text-align: right;"), titleWidth = 1200)

sidebar <- dashboardSidebar(width = 220,
  tags$style(HTML(".main-sidebar {
  background-color: #000066; color: #fff; font-weight:bold}
  .treeview-menu>li>a {font-weight:bold; color: #fff;}")),
  sidebarMenu(sidebarMenuOutput(("sidebarPanel"))))

body <-  dashboardBody(tags$style(HTML(".main-sidebar {font-size:16px;
  background-color:#000066;color:#fff;font-weight:bold}
 .treeview-menu>li>a{font-size:16px;font-weight:bold;color: #fff ;}")),

tags$style("#controls {
  background-color: #dcd0ff;opacity: 1;}
  #controls:hover{opacity: 1;}
  #econtrols {
  background-color: #dcd0ff;opacity: 1;};  }"),

tags$style("#controls {
  background-color: #dcd0ff;opacity: 1;}
  #econtrols:hover{opacity: 1;}
  #econtrols {
  background-color: #dcd0ff;opacity: 1;};  }"),

tags$head(tags$style(HTML("#point_einfo tr.selected
  {background-color:white},
  #point_info tr.selected {background-color:white}"))),

#makes pngs same size as window after rescaling
tags$script("$(document).on('shiny:connected', function(event) {
  var myWidth = $(window).width();
  Shiny.onInputChange('shiny_width',myWidth)});"),

useShinyjs(),

tags$script("$(document).on('shiny:connected', function(event) {
  var myHeight = $(window).height();
  Shiny.onInputChange('shiny_height',myHeight)});"),

tags$script(func <- JS('function(event, ui){return $(event.target).offset();}')),

tags$script(HTML("
  //Get mouse coordinates
  var mouseX, mouseY;$(document).mousemove(function(e) {
  #mouseX = e.pageX;
  #mouseY = e.pageY;}).mouseover();
  //Function to position draggable, place on current mouse coordinates
  Shiny.addCustomMessageHandler ('placeDraggable',function (message) {
  var element = $('#click_info').parent();
  #element.css({'top': mouseY + 'px', 'left' : mouseX + 'px'})});

  //Show or hide draggable
  Shiny.addCustomMessageHandler ('hideDraggable',function (message) {
  if(message.hide == true){
  $('#click_info').parent().hide();} else{
  $('#click_info').parent().show();}});

  //Get mouse coordinates
  var mouseX, mouseY;$(document).mousemove(function(e) {
  #emouseX = e.pageX;
  #emouseY = e.pageY;}).mouseover();
  //Function to position draggable, place on current mouse coordinates
  Shiny.addCustomMessageHandler ('placeDraggable',function (message) {
  var element = $('#click_info').parent();
  #element.css({'top': mouseY + 'px', 'left' : mouseX + 'px'})});

  //Show or hide draggable
  Shiny.addCustomMessageHandler ('hideDraggable',function (message) {
  if(message.hide == true){
  $('#eclick_info').parent().hide();} else{
  $('#eclick_info').parent().show();}});

  //Show or hide draggable
  Shiny.addCustomMessageHandler ('ehideDraggable',function (message) {
  if(message.hide == true){
  $('#eclick_info').parent().hide();} else{
  $('#eclick_info').parent().show();}}),")),

tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "chemogenomics.css")),

shinyjs::useShinyjs(),

tags$style(HTML(".selectize-input{white-space: nowrap;
  font-size: 14px; font-weight:bold!important;")),

#########################################################################
########################### START TABITEMS ##############################
#########################################################################

tabItems(

tabItem("hiphop", class = "active",

fluidRow(

  column(width = 3,

    box(title = "Select input datasets:",

      br(),

      prettyRadioButtons("site", label = NULL,
        choices = c(
          "2008 Erikson PLOS Genetics" = "elke",
          "2021 HIPHOP Marjan" = "marjhip",
          "2020:2024 HOP Marjan" = "marjhop"),

      outline = TRUE, fill = FALSE, thick = TRUE, shape = "square",
      bigger = FALSE, selected = "marjhop", inline = FALSE),

    status = "primary", solidHeader = TRUE, width = "100%", height = 160)),

  column(width = 9,

    box(title = "Select screen:",

  selectizeInput('cmpSERV', label = NULL,

choices = phop$name, multiple = FALSE,
  selected = "24Jan05:methyl-methanesulfonate-mms_800uM",

options = list('plugins' = list('autofill_disable'))),

  status = "primary", solidHeader = TRUE, width = "100%", height = 160))), # fluidRow ok

fluidRow(

  column(width = 3,

    box(align = "center", title = "Reset compound menu:",

      br(),

      actionButton("resetCmp", "Reset cmpMenu"),

    status = "primary", solidHeader = TRUE, width = "100%", height = 200)), # column ok

  column(width = 3,

    box(title = "Fitness score threshold:",

  sliderInput("scorethreshSERV", label = NULL, min = 0, value = 1.0,

step = 0.1, max = 5.0), # ok

  fluidRow(column(width = 12, align = "center",

actionLink(inputId = "resetScore", label = "reset",

style = "height: 100px;width: 150px;font-size:120%;
  text-align:center", size = "lg"))),

status = "primary", solidHeader = TRUE, width = "100%", height = 200)),

   uiOutput("ylimits")), # fluidRow ok

   uiOutput("hiphoppanel")), # tab ok

tabItem("signature",

fluidRow(sigNetworkModuleUI("sigNetworkModule1")),

fluidRow(

      box(title = "HOP profiles in response signature:",

DT::dataTableOutput("screenResp"),

  status = "primary", solidHeader = TRUE, width = "100%", height = 12))),

tabItem("goenrich",

fluidRow(visNetworkModuleUI("visNetwork1")))))


#####################################################################
########################### END TABITEMS ############################
#####################################################################

ui<-dashboardPage(header, sidebar, body, skin = "blue",
  setBackgroundColor(color = "#e6e6ff", shinydashboard = TRUE))

server <- function(input, output, session) {

observeEvent(input$site, {

  req(input$site)

if(input$site %in% c("elke", "marjhop")) {

output$sidebarPanel <-renderMenu({

  sidebarMenu(id = 'tabs',selected = "HOmozygous Profiles",

  menuItem("HOmozygous Profiles", tabName="hiphop",
    icon = icon("bullseye"), selected = 1),

  menuItem('GO enrichments', tabName = 'goenrich',
    icon = icon('dna')),

  menuItem('Signatures', tabName = 'signature',
    icon = icon('database')))})


} else {

output$sidebarPanel <-renderMenu({

sidebarMenu(id = "tabs", selected = HTML("HaploInsufficiency
    Profiles<br/>&nbsp&nbsp&nbsp&nbsp&nbsp&nbspHomozygous Profiles"),

  menuItem(HTML("HaploInsufficiency Profiles<br
    />&nbsp&nbsp&nbsp&nbsp&nbsp&nbspHomozygous Profiles"),

  tabName="hiphop",icon = icon("bullseye"), selected = 1),

  menuItem('GO enrichments', tabName = 'goenrich',
    icon = icon('dna')),

  menuItem('Signatures', tabName = 'signature',
    icon = icon('database')))})

}

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

  plotOutput("nfd", width = "90%", height = 700,

  click = clickOpts(id = "plot_click")),

  status = "primary", solidHeader = TRUE, width = 12),

jqui_draggable(absolutePanel(id = "controls", width = 600,
draggable = TRUE,uiOutput("click_info")),
options = list(cancel = ".shiny-input-container"))),# fluidRow ok

fluidRow(

  box(title="HOP genes: mouse click on row links to gene fitness profile;

  mouse click on genename links to SGD",

DT::dataTableOutput("hoptab"),

  status = "primary", solidHeader = TRUE,width = 12), # fluidRow ok

   box(title="download HOP:",

     column(width = 6,

       downloadButton("downloadhop", "HOPdownload")),

     column(width = 6,

       downloadButton("ExportnonPlot", "HOPfitness plot")),

   status = "primary", solidHeader = TRUE,width = 12)),

fluidRow(

  box(title="Compound information: CoInhibitory (CoI) value,
    pvalue, PubChem ID (CID) link, Mechanism of Action (MoA) link,
    indication and mechanism",

    h5(strong("Mouse click on row to view CoInhibitory (CoI)
      Homozygous Profile (HOP) and corresponding GO
      enrichments in next tab")),

    br(),

    DT::dataTableOutput('coInhib'),

    status = "primary", solidHeader = TRUE,width = 12)), # fluidRow ok

fluidRow(align = "center",

  box(title="download CoInhibition:",

  downloadButton('downloadcoinhib', 'CoInhibition'),

  status = "primary", solidHeader = TRUE,width = 3))) # tabItem ok

} else if(length(wboth) > 0){

tabItem("HIPHOP",

fluidRow(

  box(title = "HIP fitness profile:",

    plotOutput("efd", width = "90%",height = 500,

    click = clickOpts(id="plot_eclick")),

    status = "primary", solidHeader = TRUE, width = 6),

  box(title = "HOP fitness profile:
  mouse click on points to view gene description",

    plotOutput("nfd", width = "90%", height = 500,

    click = clickOpts(id="plot_click")),

    status = "primary", solidHeader = TRUE, width = 6),

jqui_draggable(absolutePanel(id = "controls", width = 600,
draggable = TRUE,uiOutput("click_info")),
options = list(cancel = ".shiny-input-container"))), # fluidRow oK

fluidRow(

  box(title="HIP genes: mouse click on row links to gene fitness profile",

    DT::dataTableOutput("hiptab", width = "100%"),

  status = "primary", solidHeader = TRUE,width = 6),

  box(title="HOP genes: mouse click on genename links to SGD",

    DT::dataTableOutput("hoptab",width = "100%"),

  status = "primary", solidHeader = TRUE, width = 6)), # fluidRow oK

fluidRow(

  box(title="download HIP:",

    column(width = 6,

      downloadButton("downloadhip", "HIPdownload")),

    column(width = 6,

      downloadButton("ExportessPlot", "HIPfitness plot")),

  status = "primary", solidHeader = TRUE,width = 6),

   box(title="download HOP:",

     column(width = 6,

      downloadButton("downloadhop", "HOPdownload")),

     column(width = 6,

       downloadButton("ExportnonPlot", "HOPfitness plot")),

  status = "primary", solidHeader = TRUE,width = 6)), # fluidRow oK

fluidRow(

  box(title="Compound information: CoInhibitory (CoI) value,
  pvalue, PubChem ID (CID) link, Mechanism of Action (MoA) link,
  indication and mechanism",

    h5(strong("Mouse click on row to view CoInhibitory
      (CoI) Homozygous Profile (HOP) and corresponding GO enrichments in next tab")),

    br(),

    DT::dataTableOutput('coInhib'),

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

    DT::dataTableOutput("hoptab"),

  status = "primary", solidHeader = TRUE,width = 12), # fluidRow ok

   box(title="download HOP:",

     column(width = 6,

    downloadButton("downloadhop", "HOPdownload")),

    column(width = 6,

    downloadButton("ExportnonPlot", "HOPfitness plot")),

  status = "primary", solidHeader = TRUE,width = 12)),

fluidRow(

  box(title = "Compound information: CoInhibitory
    (CoI) value,,pvalue, PubChem ID (CID) link,
    Mechanism of Action (MoA) link,
    indication and mechanism",

    h5(strong("Mouse click on row to view CoInhibitory (CoI)
    Homozygous Profile (HOP) and corresponding
    GO enrichments in next tab")),

    br(),

    DT::dataTableOutput('coInhib'),

  status = "primary", solidHeader = TRUE,width = 12)), # fluidRow ok

fluidRow(align = "center",

  box(title="download CoInhibition:",

    downloadButton('downloadcoinhib', 'CoInhibition'),

  status = "primary", solidHeader = TRUE,width = 3))) # tabItem ok

  }})

##############################################################################
########################### END HIPHOPPANEL ##################################
##############################################################################

##############################################################################
##########################  START CALL MODULES ##############################@
##############################################################################

tabSEND = reactiveValues(

  tab = NULL)

xinput = reactive({

  req(input$site)

  site = input$site

if(site == "marjhop"){

  xinput = d3[,phop$name]

  } else if(site == "marjhip"){

  xinput = dint

  } else {

  xinput = delke}

})

pdata = reactive({

  req(input$site)

  site = input$site

if(site == "marjhop"){

  pdata = phop

  } else if(site == "elke"){

  pdata = pmarjhip

  } else {

  pdata = pelke}

})

respSEND = reactiveValues(resp = NULL)

threshSEND = reactiveValues(thresh = NULL)

cmpSEND = reactiveValues(cmp = NULL)

observeEvent(input$tabs, {

  tabSEND$tab = input$tabs
  cmpSEND$cmp = input$cmpSERV

}, ignoreInit = FALSE, ignoreNULL = T)

returnedMOD   = callModule(visNetworkModule,'visNetwork1',xinput = xinput,
  cmpInput    = reactive(cmpSEND$cmp),
  threshInput = reactive(threshSEND$thresh),
  siteInput   = reactive(input$site),
  tabsInput   = reactive(tabSEND$tab))

cmpRETURN = reactive({
  print(paste("cmpRETURNED",returnedMOD[[1]]()))
  returnedMOD[[1]]()})

threshRETURN = reactive({
  print(paste("threshRETURNED",returnedMOD[[2]]()))
  returnedMOD[[2]]()})

siteRETURN = reactive({
  print(paste("siteRETURNED",returnedMOD[[3]]()))
  returnedMOD[[3]]()})

returnedSIG <- callModule(sigNetworkModule, 'sigNetworkModule1',
  xRespDat = reactive(xsig),
  xRespInput = reactive(respSEND$resp),
  inputTab  = reactive(tabSEND$tab))

sigRETURN = reactive({
  print(paste("sigRETURN", returnedSIG[[1]]()))
  returnedSIG[[1]]()})

tabRETURN = reactive({
 print(paste("tabservRETURN", returnedSIG[[2]]()))
 returnedSIG[[2]]()})

###################################################
################# END CALL MODULES ################
###################################################

###################################################
############### START OBSERVEEVENT ################
###################################################

observeEvent(cmpRETURN(), {

req(cmpRETURN())

tabNOThop <- input$tabs != "hiphop"

threshMOD <- threshRETURN()

cmpMOD <- cmpRETURN()

req(cmpMOD)

req(cmpMOD != "")

req(!is.null(cmpMOD))

req(!is.null(input$cmpSERV))

cmpDIFF <- cmpMOD != input$cmpSERV

req(tabNOThop == TRUE | tabNOThop == FALSE)

 if(cmpMOD %in% phop$name){

     choices = phop$name} else if(cmpMOD %in% pmarjhip$name){

     choices = pmarjhip$name} else if(cmpMOD %in% pelke$name){

     choices =  pelke$name}

  if(cmpMOD %in% phop$name){

     site = "marjhop"} else if(cmpMOD %in% pmarjhip$name){

     site = "marjhip"} else if(cmpMOD %in% pelke$name){

     site = "elke"}

req(site)

req(site != "")

req(!is.null(site))

updatePrettyRadioButtons(session,"site", label = NULL,
  c("2008 Erikson PLOS Genetics" = "elke",
    "2021 HIPHOP Marjan" = "marjhip",
    "2020:2024 HOP Marjan" = "marjhop"), selected = site)

if (cmpMOD != input$cmpSERV) {
updateSelectizeInput(session, "cmpSERV", label = NULL,
choices = choices, selected = cmpRETURN(), server = TRUE)
}

}, ignoreInit = FALSE, ignoreNULL = TRUE)

observeEvent(threshRETURN(), {

  req(threshRETURN())

  tabNOThop = input$tabs != "hiphop"

  threshMOD = threshRETURN()

  req(threshMOD)

  req(threshMOD!="")

  req(!is.null(input$cmpSERV))

  threshDIFF = threshMOD!=input$scorethreshSERV

  if(threshMOD!=input$scorethreshSERV){

   updateSliderInput(session, "scorethreshSERV", "",

   value = threshMOD, step = 0.1, min = 0, max = 5.0)
 }

}, ignoreInit = FALSE, ignoreNULL = T)

observeEvent(input$cmpSERV, {

req(input$cmpSERV)

req(!is.null(input$tabs))

req(input$tabs!= "")

updateSliderInput(session, "scorethreshSERV", label = NULL,
  value = 1.0, step = 0.1, min = 0, max = 5.0)

cmpSEND$cmp = input$cmpSERV

threshSEND$thresh = input$scorethreshSERV

tabSEND$tabs = input$tabs

}, ignoreInit = FALSE, ignoreNULL = T)

observeEvent(input$scorethreshSERV,{

  req(input$cmpSERV)

  cmpSEND$cmp = input$cmpSERV

  threshSEND$thresh = input$scorethreshSERV

  tabSEND$tabs = input$tabs

}, ignoreInit = FALSE, ignoreNULL = T)

observeEvent(input$site, {

req(input$site)

site = input$site

req(site)

if(site == "marjhop") {

choices = phop$name} else if(site == "marjhip"){

choices = pmarjhip$name} else if(site == "elke"){

choices = pelke$name}

if(site == "marjhop") {

selected = "24Jan05:methyl-methanesulfonate-mms_800uM"} else if(site == "marjhip"){

selected = "23May09:doxorubicin-hydrochloride_13uM"} else if(site == "elke"){

selected =  "sertraline hydrochloride"}

updateSelectizeInput(session,"cmpSERV", label = NULL,

choices = choices, selected = selected, server = T)

}, ignoreNULL = TRUE, ignoreInit = FALSE)

#######################################################
################## END OBSERVEEVENTS ##################
#######################################################

#######################################################
############# START SIGNATURE EVENTS ##################
#######################################################

output$screenResp = DT::renderDataTable({

  w = which(phop$signature %in% sigRETURN())

  validate(need(length(w) > 0 , message =
  "please enter a valid signature"))

  df = data.frame(

  screen = phop$name[w],

  mechanism = phop$mechanism[w],

  PCID = phop$PCID[w],

  compound = phop$cmp[w],

  stringsAsFactors = FALSE)

  opts = list(

  pageLength = 25,

  autoWidth = FALSE,

  scrollX = FALSE,

  columnDefs = list(

  list(className = 'dt-left', targets = c(0, 1, 2, 3)),

  list(width = c('100px'), targets = c(2)),

  list(width = c('500px'), targets = c(1)),

  list(width = c('300px'), targets = 0)))

  df = df %>% dplyr::arrange(compound)

  df =  DT::datatable(df,

  options = opts, rownames = FALSE, escape = FALSE,
  selection = "single") %>%

  formatStyle(c(1:4), fontWeight = 'bold', fontSize = '14px')

})

observeEvent(input$screenResp_rows_selected, {

  w = which(phop$signature %in% sigRETURN())

  validate(need(length(w) > 0 , message = "please enter a valid signature"))

  d = phop$name[w]

  df = data.frame(

  screen = phop$name[w],

  mechanism = phop$mechanism[w],

  PCID = phop$PCID[w],

  compound = phop$cmp[w],

  stringsAsFactors = FALSE)

  df = df %>% dplyr::arrange(compound)

  row = input$screenResp_rows_selected

  w = which(df$screen %in% df$screen[row])

  validate(need(length(w) > 0 , message = "please enter a valid signature"))

  selected = df$screen[row]

  cmpSEND$cmp = selected

  updateSelectizeInput(session, "cmpSERV", label = NULL,

choices = phop$name, selected = selected, server = T)

  newtab <- switch(input$tabs,

   "signature" = "hiphop",

   "hiphop" = "signature")

  updateTabItems(session, "tabs", newtab)

})

output$targethip <- DT::renderDataTable({

whna = which(!is.na(phop$pcid))
cids = phop$pcid[whna]
pdata = phop[whna,]
image_paths <- sapply(cids, function(cid) {
paste0("hopimg/", cid, ".png")
 })
  # Use the helper function to generate the img tags
img_tags <- sapply(image_paths, function(path) {
  generate_img_tag(path)
})

w = which(pdata$name %in% input$cmpSERV)

d = pdata[w, c(
"name",
"cmp",
"signature",
"PCID",
"images",
"mechanism"
)]

  names(d)[5] = "structure"

  cids = pdata$pcid[w]

  d$structure = img_tags[w]

  PCID = '<a href="https://pubchem.ncbi.nlm.nih.gov" target="_blank">PCID</a>'

  names(d)[1] = "screen"

  names(d)[4] = PCID

  names(d)[2] = "compound"

opts = list(dom = 'Bfrtip', paging = FALSE, target = "cell", searching = FALSE,
  info = FALSE, autowidth = TRUE, scrollX = TRUE, ordering = FALSE,
  columnDefs = list(
  list(className = 'dt-left', targets = c(0,1,4)),
  list(className = 'dt-center', targets = c(2,3)),
  list(width = c('200px'), targets = c(0,2)),
  list(width = c('100px'), targets = c(1,3)),
  list(width = c('400px'), targets = c(5))))

datatable(d, options = opts, escape = FALSE,
  class = 'table-bordered stripe table-condensed',
  rownames = FALSE)

})

observeEvent(input$targethip_rows_selected, {

whna = which(!is.na(phop$pcid))
cids = phop$pcid[whna]
pdata = phop[whna,]
w = which(pdata$name %in% input$cmpSERV)

row = input$targethip_rows_selected

respSEND$resp = pdata$signature[w][row]

newtab <- switch(input$tabs,
 "hiphop" = "signature",
 "signature" = "hiphop")

updateTabItems(session, "tabs", newtab)

})

#######################################################
############# END SIGNATURE EVENTS ####################
#######################################################

############## RESET COMPOUND #############

observeEvent(input$resetCmp,{

  site = input$site

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
   selected = "sertraline hydrochloride",server = T)
  }


  },ignoreInit = FALSE, ignoreNULL = T)

##################################################
################ START MOUSE-OVERS ###############
##################################################

show_click = reactiveVal(NULL)
print_click = reactiveVal(NULL)
output$click_info <- renderUI(show_click())
output$point_info <- DT::renderDataTable({

  df = print_click()

  opts = list(dom = 'Bfrtip', paging = FALSE, target = "cell", searching = FALSE,
info = FALSE, autowidth = TRUE, scrollX = TRUE, ordering = FALSE,
columnDefs = list(
  list(width = c('10%'), targets = c(0, 1, 2)),
  list(width = c('70%'), targets = 3)))

  df =  DT::datatable(df, options = opts,rownames = FALSE, escape = FALSE, selection = "single") %>%
  formatStyle(c(1:4),fontWeight = 'bold', fontSize = '12px',target="row")

})

observeEvent(input$plot_click, {

   welke = which(colnames(delke)%in% input$cmpSERV)
   whop = which(colnames(d3) %in% input$cmpSERV)
   wboth = which(colnames(dn3) %in% input$cmpSERV)

  if(length(whop) > 0) {x = d3} else if(length(wboth) > 0) {x = dn3} else {x = delke}

  req(input$cmpSERV)
  w1 = which(colnames(x) == input$cmpSERV)

  req(length(w1) > 0)

  validate(need(length(w1) == 1 ,message = "please enter a valid compound"))


  hop <- geneAnno(mat = x, cmp = colnames(x)[w1], fdat = fdat, sgdlink = T)

  pclick <- nearPoints(df = hop,  xvar = "gene",
                       yvar = "FD", coordinfo = input$plot_click,threshold = 10)

  hideTooltip <- function(hide){
  session$sendCustomMessage(type = 'hideDraggable', message = list('hide'= hide))
  }

  g = grep("FD",names(pclick))

  if(length(g) > 0) pclick[,g] = format(round(pclick[,g],2),nsmall = 1,scientific = FALSE)

  g = which(names(pclick) %in% c("xvar","gene"))

  if(length(g) > 0) pclick = pclick[,-g]

  if( nrow(pclick) == 0 ) {
  show_click(NULL)

  hideTooltip(TRUE) # Hide tooltip if there's no info to show
  return()

  } else {

  session$sendCustomMessage(type = 'placeDraggable', message = list())
  show_click(tagList(

  {DT::dataTableOutput("point_info", width = "100%")}))

  print_click({pclick})

  }

})

show_eclick = reactiveVal(NULL)
print_eclick = reactiveVal(NULL)
output$eclick_info <- renderUI(show_eclick())
output$point_einfo <- DT::renderDataTable({

  df = print_eclick()

 opts = list(dom = 'Bfrtip', paging = FALSE, target = "cell", searching = FALSE,
   info = FALSE, autowidth = TRUE, scrollX = TRUE, ordering = FALSE,
   columnDefs = list(
  list(width = c('10%'), targets = c(0, 1, 2)),
  list(width = c('70%'), targets = 3)))

  df =  DT::datatable(df, options = opts,rownames = FALSE, escape = FALSE, selection = "single") %>%
formatStyle(c(1:4),fontWeight = 'bold', fontSize = '12px',target="row")

})

observeEvent(input$plot_eclick, {

  x = de3

  req(input$cmpSERV)

  w1 = which(colnames(x) == input$cmpSERV)

  req(length(w1) > 0)

  validate(need(length(w1) == 1 ,message = "please enter a valid compound"))

  hideTooltip <- function( hide ){
  session$sendCustomMessage(type = 'ehideDraggable', message = list('hide'=hide))
  }

  hop <- geneAnno(mat = x, cmp = colnames(x)[w1],
    fdat = fdat, sgdlink = T)

  peclick <- nearPoints(df = hop,  xvar = "gene", yvar = "FD", coordinfo = input$plot_eclick,threshold = 5)

  g = grep("FD",names(peclick))

  if(length(g) > 0) peclick[,g] = format(round(peclick[,g],2),nsmall = 1,scientific = FALSE)

  g = which(names(peclick) %in% c("xvar","gene"))

  if(length(g) > 0) peclick = peclick[,-g]

  if( nrow(peclick) == 0 ) {
show_eclick(NULL)

hideTooltip(TRUE) # Hide tooltip if there's no info to show
return()

 } else {

  session$sendCustomMessage(type = 'eplaceDraggable', message = list())
  show_click(tagList(

  {DT::dataTableOutput("point_info", width = "100%")}))

  print_click({peclick})

  }

})



##################################################
################### START COFITNESS ##############
##################################################

output$coFitNess = DT::renderDataTable({

  outFit$dFit

},escape = FALSE,rownames= FALSE,options =
  list(pageLength=10), selection = "single", server = FALSE)

observeEvent(input$coFitNess_rows_selected, {

  row = input$coFitNess_rows_selected

  xfit = outFit$dFit

   w = xfit$gene[row]

   validate(need(length(w) == 1 , message = "please enter a valid gene"))

   updateTextInput(session, inputId = "gene",
     label = NULL, value = xfit$gene[row])

})
####################################################

output$downloadcofitness <- downloadHandler(
   filename = function() {
  paste0("CoFit:",input$gene, "_",Sys.Date(), ".txt")
  },
  content = function(file) {

  write.table(as.data.frame(outFit$dFit),
    file, row.names = FALSE,sep="\t",quote=FALSE)
  }
)

##################################################
###################### END COFITNESS #############
##################################################

#######################################################
#######################################################
################## START COINHIBITION #################
#######################################################
#######################################################

outcoInhib = reactive({

  # Create a data frame with CIDs and corresponding image paths

  req(input$cmpSERV)

  if(input$cmpSERV%in% phop$name) {
    whna = which(!is.na(phop$pcid))
    cids = phop$pcid[whna]
    xcov  = xcoinh
    xpv   =  xpval
    pdata = phop[whna,]
    image_paths <- sapply(cids, function(cid) {
    paste0("hopimg/", cid, ".png")
     })
  # Use the helper function to generate the img tags
    img_tags <- sapply(image_paths, function(path) {
      generate_img_tag(path)
    })

} else if(input$cmpSERV%in% pmarjhip$name){
    wina = which(!is.na(pmarjhip$pcid))
    cids = pmarjhip$pcid[wina]
    xcov  = xcohh
    xpv   = xpvhh
    pdata = pmarjhip[wina,]
    image_paths <- sapply(cids, function(cid) {
    paste0("hipimg/", cid, ".png")
  })
    # Use the helper function to generate the img tags
    img_tags <- sapply(image_paths, function(path) {
      generate_img_tag(path)
    })
} else {
    cids = pelke$pcid
    xcov  = xelke
    xpv   = xelkepv
    pdata = pelke
    image_paths <- sapply(cids, function(cid) {
    paste0("elkeimg/", cid, ".png")
  })
    # Use the helper function to generate the img tags
    img_tags <- sapply(image_paths, function(path) {
      generate_img_tag(path)
    })
}

w = which(colnames(xcov) == input$cmpSERV)

validate(need(length(w) == 1 ,message =
 "please enter a valid condition"))

  p = xpv[,w,drop=F]
  d = xcov[,w,drop = F]

  d = d[order(d[,1],decreasing = T),,drop=F]

  df = data.frame(screen = rownames(d)[1:nrow(d)],
  CoI = d[1:nrow(d),1],pvalue = p[rownames(d),1],
  stringsAsFactors = FALSE)

  df$pvalue = formatC(df$pvalue,format = "e", digits = 2)

  w = which(df$CoI >= 0)

  if(length(w) > 0) df = df[w,]

  m = match(df$screen,pdata$name)

  df$CoI = format(round(df$CoI,2),nsmall = 1,scientific = FALSE)

  df$CID = pdata$PCID[m]
  cids = pdata$pcid[m]
  df$structure = img_tags[m]

  w = which(names(df) == "CID")

  names(df)[w] = '<a href="https://pubchem.ncbi.nlm.nih.gov" target="_blank">CID</a>'

  df$mechanism = pdata$mechanism[m]
  df

   })

output$coInhib = DT::renderDataTable({

xcoinhib = outcoInhib()

opts = list(
dom = 'Bfrtip',
pageLength = 10,
autoWidth = FALSE,
scrollX = TRUE)


 DT::datatable( xcoinhib, options = opts,
 escape = FALSE,  rownames = FALSE, selection = "single")

  })

observeEvent(input$coInhib_rows_selected, {

  row <- input$coInhib_rows_selected

  df = outcoInhib()

  updateSelectizeInput(session,"cmpSERV",label = NULL,
    choices = df$screen, selected = df$screen[row][1],server = T)

})

output$downloadcoinhib <- downloadHandler(
filename = function() {
  paste0("HOP:",
 input$cmpSERV,
 "_",Sys.Date(), ".txt")
},
content = function(file) {
  write.table(as.data.frame(outcoInhib()), file, row.names = FALSE,sep="\t",quote=FALSE)
})

##################################################
################# END COINHIBITION ###############
##################################################

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

observeEvent(input$network_Proxy_selected, {

req(input$network_Proxy_selected)

n = visNet()$nodes

w = which(n$id %in% as.numeric(input$network_Proxy_selected))

id = n$id[w]

visNetworkProxy("network_Proxy") %>%

  visSelectNodes(id = id)

}, ignoreNULL = T, ignoreInit = FALSE)


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

######################################################################
######################### START BARPLOTS #############################
######################################################################

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

datatable(t, width=220,caption = htmltools::tags$caption(term,
  style = "caption-side: top; text-align: center; color:black;b
  background:white;font-weight:bold;"),

  options=list(paging= FALSE, scrollY = FALSE, dom= "t",scroller = FALSE,
    searching = FALSE,ordering = FALSE, rowCallback = JS(
    "function(row, data) {",
    "for (i = 1; i < data.length; i++) {",
    "if (data[i]>1000 | data[i]<1){",
    "$('td:eq('+i+')', row).html(data[i].toExponential(1));","}","}","}")),

  height = 400, colnames = "") %>%

  formatStyle( target = "row", color = "black",backgroundColor = "white",
  columns = c(1,2),fontWeight = "bold")

  })

output$LeadingEdge = renderUI({

req(input$network_Proxy_selected)

plotOutput("BarPlot",
  width = 300,
  height = 300)
  })

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

  barplot(d$score,names.arg = d$gene,las = 1, horiz = TRUE, cex.axis = 0.8,
    cex.lab = 1.1, col ="navy",xlab = "fitness defect score")

   })

######################################################################
######################### END BARPLOTS ###############################
######################################################################

####################################################
################ START HIPHOP TABLES ###############
####################################################

output$hiptab = DT::renderDataTable({

  req(input$cmpSERV)

  mat = de3

  w1 = which(colnames(mat) == input$cmpSERV)

  hop <- hgeneAnno(mat = mat, cmp = colnames(mat)[w1], fdat = orth)

  hop <- hop[, c("FD", "gene", "GeneCards", "descriptor")]

  opts = list(pageLength=10, autoWidth = TRUE,scrollX = TRUE,
    columnDefs = list(list(width = "80px", targets = c(-2,-3))))


  hip

  },escape = FALSE,

  options=list(pageLength=5, autoWidth = FALSE,scrollX = TRUE,
 columnDefs = list(list(width = "80px", targets = c(-2,-3)))),rownames = FALSE)

###
###
###
outhip = reactive({

req(input$cmpSERV)

  hop <- hgeneAnno(mat = de3, cmp = input$cmpSERV, fdat = orth)

  hop <- hop[, c("FD", "gene", "GeneCards", "descriptor")]

  m = match(hop$gene,orth$human_gene)

  hop$yeast_homolog = orth$yeast_gene[m]
  hop$SGD = orth$yeast_raw[m]
  hop$yeast_descriptor = orth$yeast_anno[m]
  hop$GeneCards = orth$human_raw[m]
  hop

  })

output$downloadhip <- downloadHandler(
  filename = function() {
    paste0("HIP:",input$cmpSERV,  "_",Sys.Date(), ".txt")
  },
  content = function(file) {
    write.table(as.data.frame(outhip()), file, row.names = FALSE,
    sep = "\t",quote = FALSE)
  }
  )

output$efd  = renderPlot({

  req(input$cmpSERV)
  validate(need(input$cmpSERV,message =
   "please select condition"))

  cond1 = input$cmpSERV

  x = de3

  wx = which(colnames(x) %in% cond1)

  req(length(wx)!=0)
  req(input$site == "marjhip")
  colnames(x)[wx] = paste0("HIP|",colnames(x)[wx])
  coln = colnames(x)[wx]
  coln   = str_wrap(coln, width=30,whitespace_only = FALSE)
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
##########################
output$ExportessPlot <- downloadHandler(
  # file name
  filename = function() {
  paste0("HIPplot:",input$cmpSERV,  "_",Sys.Date(), ".png")
  },
  content = function(file){
  # create plot
  png(file,
  width = (input$shiny_width/12)*6,
  height = 1000)

  essPlot()
  dev.off()

  })

output$hoptab = DT::renderDataTable({

req(input$cmpSERV)

if(input$site == "marjhop"){
  mat = d3[,phop$name]} else if(input$site == "marjhip"){
  mat = dn3} else {
  mat = delke}

  hop <- hgeneAnno(mat = mat, cmp = input$cmpSERV, fdat = orth)

  hop <- hop[, c("FD", "gene", "GeneCards", "descriptor")]

  hop$FD <- format(round(hop$FD, 2), nsmall = 1, scientific = FALSE)

  hop

}, escape = FALSE, options = list(pageLength = 10, autoWidth = FALSE,
scrollX = TRUE, columnDefs = list(list(width = "80px",
targets = c(-2, -3)))), rownames = FALSE)

outhop = reactive({

req(input$cmpSERV)

if(input$site == "marjhop"){mat = d3[,phop$name]} else if(
input$site == "marjhip"){mat = dn3} else (mat = delke)


  hop <- hgeneAnno(mat = x, cmp = input$cmpSERV, fdat = orth)

  hop <- hop[, c("FD", "gene", "GeneCards", "descriptor")]

  m = match(hop$gene,orth$human_gene)

  hop$yeast_homolog = orth$yeast_gene[m]
  hop$SGD = orth$yeast_raw[m]
  hop$yeast_descriptor = orth$yeast_anno[m]
  hop$GeneCards = orth$human_raw[m]
  hop

  })

observeEvent(input$hoptab_rows_selected, {

  row <- input$hoptab_rows_selected

  df = outhop()

  updateTextInput(session, inputId = "gene", label = NULL, value = df$gene[row])

})

output$downloadhop <- downloadHandler(
  filename = function() {
    paste0("HOP:",input$cmpSERV,  "_",Sys.Date(), ".txt")
  },
  content = function(file) {
    write.table(as.data.frame(outhop()), file, row.names = FALSE,
                sep = "\t",quote = FALSE)
  }
)

output$nfd = renderPlot({

  req(input$site)
  req(input$cmpSERV)

  validate(need(input$cmpSERV, message = "please select condition"))

  cond1 = input$cmpSERV

  if(input$site == "marjhip") {x = dn3}
  if(input$site == "marjhop") {x = d3}
  if(input$site == "elke") {x = delke}

  wx = which(colnames(x) %in% cond1)

  colnames(x)[wx] = paste0("HOP|",colnames(x)[wx])

  if(input$site == "marjhip"){
    coln = colnames(x)[wx]
    coln   = str_wrap(coln, width=30,whitespace_only = FALSE)
    colnames(x)[wx] = coln
  }

  req(length(wx) > 0)
  p10(x,wx,pch = 17, ylim = input$ylim,sig = input$scorethreshSERV)

})
##########################
nonPlot = function(){
  validate(
  need(input$cmpSERV,message =
   "please select condition"))
  cond1 = input$cmpSERV


if(input$site == "marjhop"){
  x = d3[,phop$name]} else if(input$site == "marjhip"){
  x = dn3} else {
  x = delke}

  wx = which(colnames(x) %in% cond1)

  colnames(x)[wx] = paste0("HOP|",colnames(x)[wx])

  p10(x,wx,pch = 17)

}
##########################
output$ExportnonPlot <- downloadHandler(

  filename = function() {
  paste0("HOPplot:",input$cmpSERV,  "_",Sys.Date(), ".png")
  },
  content = function(file){
  # create plot
  png(file,
  width = (input$shiny_width/12)*6,
  height = 1000)

  nonPlot()
  dev.off()

  })

#######################################################
#######################################################
############# END HOPTABLES  AND PLOTS ################
#######################################################
#######################################################

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

wne = which(rownames(x) %in% noness$sgd_gene)

we = which(rownames(x) %in% ess$sgd_gene)

if(length(wne) > 0) {

  ylims = c(floor(min(x[wne,wcoln,drop=F],na.rm = T)),
    ceiling(max(x[wne,wcoln,drop=F],na.rm = T)))
}

req(length(wne)>0)

if(length(we) > 0) {

  elims = c(floor(min(x[we,wcoln,drop=F],na.rm = T)),
    ceiling(max(x[we,wcoln,drop=F],na.rm = T)))

}

tagList(

if(length(we) > 0) {

  column(width = 3,

    box(title= "HIP y-axis range:",

      fluidRow(column(width = 12,
        sliderInput("elim", label = NULL, min = elims[1],
        max = elims[2], value = c(-0.5,elims[2])))),

      fluidRow(column(width = 12, align = "center",

        actionLink(
          inputId = "eorig",
          label = "reset",
          style = "height: 100px;width: 150px;font-size:120%;
          text-align:center", size = "lg"))),

    status = "primary", solidHeader = TRUE, width = "100%", height = 200))

},

if(length(ylims) == 2){

  column(width = 3,

    box(title= "HOP y-axis range:",

      fluidRow(column(width = 12,
        sliderInput("ylim", label = NULL, min = ylims[1],
        max = ylims[2], value = c(-0.5, ylims[2])))),

      fluidRow(column(width = 12,align = "center",

        actionLink(
          inputId = "yorig",
          label = "reset",
          style = "height: 100px;width: 150px;font-size:120%;
          text-align:center", size = "lg"))),

    status = "primary", solidHeader = TRUE, width = "100%", height = 200))

}
  )

})

observeEvent(input$yorig,{

x = xinput()
coln = colnames(xinput())

wcoln = which(coln %in% input$cmpSERV)

wne = which(rownames(x) %in% noness$sgd_gene)

if(length(wne) > 0) {

xne = xinput()
ylims = c(floor(min(xne[wne,wcoln,drop=F],na.rm = T)),
  ceiling(max(xne[wne,wcoln,drop=F],na.rm = T)))

}


updateNumericRangeInput(session, "ylim", "", value = ylims)

 },ignoreInit = FALSE, ignoreNULL = T)


observeEvent(input$eorig,{

  cole = colnames(xinput())
  xe = xinput()
  wcole = which(cole %in% input$cmpSERV)

  we = which(rownames(xe)%in% ess$sgd_gene)

  if(length(we) > 0) {elims = c(floor(min(xe[we,wcole,drop=F],na.rm = T)),
    ceiling(max(xe[we,wcole,drop=F],na.rm = T)))}

   updateNumericRangeInput(session, "elim", "",value = elims)

})

observeEvent(input$resetscore, {

updateSliiderInput(session, "scorethreshSERV", "", value = 1.0,
  step = 0.1,min = 0, max = 5.0)


},ignoreInit = FALSE, ignoreNULL = T)

#######################################################
#######################################################
################## END Y-AXIS LIMITS ##################
#######################################################
#######################################################

 }
##########################
shinyApp(ui, server)
