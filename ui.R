
library(DT)
library(shinythemes)
library(shinyBS)
library(shinyWidgets)
library(data.table)


key.fl <- fread("key.txt", data.table = F)

footerTagList <- list(
  tags$footer(id = "myFooter",
              shiny::includeHTML("www/footer.html")
  )
)

shinyUI(
  fluidPage(
    titlePanel(
      title = div(
        img(src = "headerN.png"),
        span("funRiceGenes:", style = "font-size:36px;color:white;"),
        span(
          "A comprehensive database of functionally characterized rice genes", style = "font-size:28px;color:white;"
        ), style = "background-color:#33CCCC;margin-left: -15px;margin-right: -15px;margin-top: -20px;margin-bottom: -10px;"
      ), windowTitle = "Welcome to funRiceGenes!"
    ), 
    
    includeCSS("www/footer.css"),
    
    navbarPage(title="", windowTitle="Welcome to funRiceGenes!",
               footer = footerTagList,
               
               ## gene panel
               tabPanel(HTML("<strong style='font-size:18px'>Gene</strong>"),
                        icon = icon("cube", class = NULL, lib = "font-awesome"),
                        tabsetPanel(id = "infor_gene",
                                    tabPanel(HTML("<strong style='font-size:18px'>Search</strong>"),
                                             fluidRow(
                                               column(12,
                                                      radioButtons("query", 
                                                                   label=h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Search by a single gene symbol or a single genomic locus</font>'),
                                                                            bsButton("q1", label="", icon=icon("question"), style="info", size="small")),
                                                                   choices=c("MSU Locus", "RAPdb Locus", "Gene Symbol"), inline=TRUE),
                                                      
                                                      bsPopover("q1", "Fill in the table cells (case insensitive) to fetch information on cloned rice genes.", trigger = "focus"),
                                                      uiOutput("inText"),
                                                      tabsetPanel(
                                                        tabPanel('Information', DT::dataTableOutput("mytable1")),
                                                        tabPanel('Reference', DT::dataTableOutput("mytable2")),
                                                        tabPanel('Accession', DT::dataTableOutput("mytable3")),
                                                        tabPanel('Expression', DT::dataTableOutput("mytable4")),
                                                        tabPanel('Keyword', DT::dataTableOutput("mytable5")),
                                                        tabPanel('Connection', DT::dataTableOutput("mytable6")),
                                                        tabPanel('RiceNet', DT::dataTableOutput("mytable13"))
                                                      )
                                               )
                                             )
                                    ),
                                    
                                    tabPanel(HTML("<strong style='font-size:18px'>Submit</strong>"),
                                             fluidRow(
                                               column(6, h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Submit a new gene or submit new information for an existing gene</font>'),
                                                            bsButton("qq1", label="", icon=icon("question"), style="info", size="small")),
                                                      bsPopover("qq1","Note: To submit a new publication, fill in the Pubmed ID cell. To submit a new publication for an existing gene, fill in the Gene symbol and Pubmed ID cells. To submit a new gene, fill in all the four cells.",
                                                                trigger = "focus"),
                                                      wellPanel(
                                                        column(12, textInput('symsub7', strong("Gene symbol"),value="")),
                                                        column(12, textInput('msusub7', strong("MSU genomic locus"),value="")),
                                                        column(12, textInput('rapsub7', strong("RAPdb genomic locus"),value="")),
                                                        column(12, textInput('pubmed7', strong("Pubmed ID"),value="")),
                                                        column(12, passwordInput('key7', strong("Password"),value="")),
                                                        actionButton("submit7", strong("Submit"),style = "color: white; background-color: black"),
                                                        actionButton("clear1", strong("Reset"),style = "color: white; background-color: black")
                                                      )
                                               ),
                                               column(6,
                                                      h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Edit existing gene information</font>'),
                                                         bsButton("qq2", label="", icon=icon("question"), style="info", size="small") 
                                                      ),
                                                      bsPopover("qq2", "Note: Fill in all the table cells to edit either of the three items: symbol, MSU and RAPdb locus.", trigger = "focus"),
                                                      wellPanel(
                                                        column(12, textInput('oldsym', strong("Old gene symbol"),value="")),
                                                        column(12, textInput('newsym', strong("New gene symbol"),value="")),
                                                        column(12, textInput('oldmsu', strong("Old MSU locus"),value="")),
                                                        column(12, textInput('newmsu', strong("New MSU locus"),value="")),
                                                        column(12, textInput('oldrap', strong("Old RAPdb locus"),value="")),
                                                        column(12, textInput('newrap', strong("New RAPdb locus"),value="")),
                                                        column(12, passwordInput('key6', strong("Password"),value="")),
                                                        actionButton("submit6", strong("Submit"),style = "color: white;background-color: black" ),
                                                        actionButton("clear2", strong("Reset"),style = "color: white;background-color: black")
                                                      )
                                               )
                                             )
                                    )
                        )
               ),
               
               ## gene family panel
               tabPanel(HTML("<strong style='font-size:18px'>Gene Family</strong>"),
                        icon = icon("cubes", class = NULL, lib = "font-awesome"),
                        
                        tabsetPanel(id = "infor_genefamily",
                                    tabPanel(HTML("<strong style='font-size:18px'>Search</strong>"),
                                             fluidRow(
                                               column(12,
                                                      radioButtons("queryfam", h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Search by a single gene symbol or a single genomic locus</font>'),
                                                                                  bsButton("q2", label="", icon=icon("question"), style="info", size="small")),
                                                                   choices=c("MSU Locus", "RAPdb Locus", "Gene Symbol"), inline=TRUE),
                                                      
                                                      bsPopover("q2", "Fill in the table cells (case insensitive) to fetch information on rice genes.", trigger = "focus"),
                                                      uiOutput("inTextfam"),
                                                      tabsetPanel(
                                                        tabPanel(strong('Information'), DT::dataTableOutput("mytable8")),
                                                        tabPanel('Reference', DT::dataTableOutput("mytable9"))
                                                      )
                                               )
                                             )
                                    ),
                                    
                                    tabPanel(HTML("<strong style='font-size:18px'>Submit</strong>"),
                                             fluidRow(
                                               column(12,
                                                      h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Submit a new gene family</font>'),
                                                         bsButton("qq3", label="", icon=icon("question"), style="info", size="small")  
                                                      ),
                                                      bsPopover("qq3", "Note: The file Gene Family info should contain 5 columns with names Accession, Symbol, MSU,RAPdb and Name. Different columns should be separated by \t. And the 1st row should be the column names.",
                                                                trigger = "focus"),
                                                      wellPanel(
                                                        column(4, textInput('pubmed10', strong("Pubmed ID"),value="")),
                                                        column(4, passwordInput('key10', strong("Password"),value="")),  
                                                        
                                                        column(4, fileInput('genfamin', strong("Gene Family info"),
                                                                            accept=c(".txt"))),
                                                        actionButton("submit10", strong("Submit"),style = "color: white;background-color: black"),
                                                        actionButton("clear6", strong("Reset"),style = "color: white;background-color: black")
                                                      )
                                               )
                                             )
                                    )
                        )
               ),
               
               ## keyword panel
               tabPanel(HTML("<strong style='font-size:18px'>Keyword</strong>"),
                        icon = icon("key", class = NULL, lib = "font-awesome"),
                        
                        tabsetPanel(id = "infor_key",
                                    tabPanel(HTML("<strong style='font-size:18px'>Search</strong>"),
                                             fluidRow(
                                               column(4,
                                                      shinyWidgets::multiInput(
                                                        inputId = "keyword",
                                                        label = tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Search with predefined keywords</font>')),
                                                        choices = NULL, width="100%", 
                                                        choiceNames = unique(key.fl$Keyword),
                                                        choiceValues = unique(key.fl$Keyword),
                                                        selected = "nitrogen",
                                                        options = list(
                                                          enable_search = TRUE, limit = 1,
                                                          non_selected_header = "Choose from:",
                                                          selected_header = "You have selected:"
                                                        )
                                                      )
                                               ),
                                               
                                               column(8,
                                                      tabPanel('Information', DT::dataTableOutput("mytable7"))
                                               )
                                             )
                                    ),
                                    
                                    tabPanel(HTML("<strong style='font-size:18px'>Submit</strong>"),
                                             fluidRow(
                                               column(12,
                                                      h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Submit new keywords on genes</font>')),
                                                      
                                                      wellPanel(
                                                        column(2, textInput('symsub4', strong("Gene symbol"),value="")),
                                                        column(2, textInput('keysub4', strong("Keyword"),value="")),
                                                        column(2, textInput('tilsub4', strong("Title"),value="")),
                                                        column(2, textInput('evisub4', strong("Evidence"),value="")),
                                                        column(4, passwordInput('key4', strong("Password"),value="")),
                                                        actionButton("submit4", strong("Submit"),style = "color: white;background-color: black"),
                                                        actionButton("clear3", strong("Reset"),style = "color: white;background-color: black")
                                                      )
                                               )
                                             )
                                    )
                        )
               ),
               
               ## publication panel
               tabPanel(HTML("<strong style='font-size:18px'>Publication</strong>"),
                        icon = icon("book", class = NULL, lib = "font-awesome"),
                        
                        tabsetPanel(id = "infor_public",
                                    tabPanel( HTML("<strong style='font-size:18px'>Search</strong>"),
                                              fluidRow(
                                                column(12,
                                                       h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Search with any word concerning rice functional genomic studies</font>')),
                                                       textInput("publication", label = "",                              
                                                                 value = "heading date"),
                                                       tabsetPanel(
                                                         tabPanel('Result', DT::dataTableOutput("mytable11"))
                                                       )
                                                )
                                              )
                                    ),
                                    
                                    tabPanel(HTML("<strong style='font-size:18px'>Submit</strong>"),
                                             fluidRow(
                                               column(12,
                                                      h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Submit a new publication</font>')),
                                                      
                                                      wellPanel(
                                                        column(2, textInput('symsub2', strong("Gene symbol"),value="")),
                                                        column(1, textInput('tilsub2', strong("Title"),value="")),
                                                        column(1, textInput('yearsub2', strong("Year"),value="")),
                                                        column(2, textInput('jousub2', strong("Journal"),value="")), 
                                                        column(2, textInput('afisub2', strong("Affiliation"),value="")),
                                                        column(2, textInput('abssub2', strong("Abstract"),value="")),
                                                        column(2, passwordInput('key2', strong("Password"),value="")),
                                                        actionButton("submit2", strong("Submit"),style = "color: white;background-color: black"),
                                                        actionButton("clear4", strong("Reset"),style = "color: white;background-color: black")
                                                      )
                                               )
                                             )
                                    )
                        )
               ),
               
               ## ID conversion panel
               tabPanel(HTML("<strong style='font-size:18px'>ID Conversion</strong>"),
                        icon = icon("exchange", class = NULL, lib = "font-awesome"),
                        tabsetPanel(id = "infor_idc",
                                    tabPanel(HTML("<strong style='font-size:18px'>Convert between MSU and RAPdb gene IDs</strong>"),
                                             fluidRow(
                                               column(12,
                                                      radioButtons("queryconv", h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Convert between MSU genomic locus and RAPdb genomic locus</font>'),
                                                                                   bsButton("q3", label="", icon=icon("question"), style="info", size="small")),
                                                                   choices=c("RAPdb to MSU", "MSU to RAPdb"), inline=TRUE),
                                                      
                                                      bsPopover("q3", "You can submit multiple IDs delimited by space.", trigger = "focus"),
                                                      uiOutput("inTextconv"),
                                                      
                                                      tabsetPanel(
                                                        tabPanel(strong('Result'), DT::dataTableOutput("mytable10"))
                                                      )
                                               )
                                             )
                                    ),
                                    
                                    tabPanel(HTML("<strong style='font-size:18px'>Convert between indica and japonica gene IDs</strong>"),
                                             fluidRow(
                                               column(12,
                                                      radioButtons("IJconv", h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Convert between indica and japonica gene IDs</font>'),
                                                                                bsButton("q31", label="", icon=icon("question"), style="info", size="small")),
                                                                   choices=c("MSU Nipponbare", "RAPdb Nipponbare", "Minghui 63", "Zhenshan 97"), inline=TRUE),
                                                      
                                                      bsPopover("q31", "You can submit multiple IDs delimited by space.", trigger = "focus"),
                                                      uiOutput("inIJconv"),
                                                      
                                                      tabsetPanel(
                                                        tabPanel(strong('Result'), DT::dataTableOutput("mytable12"))
                                                      )
                                               )
                                             )
                                    )
                        )
               ),
               
               ## Bulk download panel
               tabPanel(HTML("<strong style='font-size:18px'>Download</strong>"),
                        icon = icon("download", class = NULL, lib = "font-awesome"),
                        sidebarPanel(
                          selectInput("selextractdata", h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Extract  data</font>'),
                                                           bsButton("bs00", label="", icon=icon("question"), style="info", size="small")
                          ), c("Extract data using MSU genomic locus" = "1","Extract data using RAP genomic locus" = "2"), "1"),
                          bsPopover("bs00", 'Select extract data type',"Select how to extract the data",trigger = "focus"),
                          conditionalPanel(condition="input.selextractdata == '1'",
                                           textAreaInput("msuarea", label = NULL, width="600px", resize="vertical", height="200px", 
                                                         placeholder = "One locus in one row",value = "LOC_Os07g49460\nLOC_Os06g40780\nLOC_Os05g06660")
                          ),
                          conditionalPanel(condition="input.selextractdata == '2'",
                                           textAreaInput("raparea", label = NULL, width="600px", resize="vertical", height="200px", 
                                                         placeholder = "One locus in one row",value = "Os02g0677300\nOs10g0528100\nOs06g0610350")
                                           
                          ),
                          actionButton("submit_ID", strong("Submit"), style = "color: white;background-color: black"),
                          actionButton("example_ID", strong("Load Example"), style = "color: white;background-color: black"),
                          actionButton("Clear_ID", strong("Reset"), style = "color: white;background-color: black")
                        ),
                        
                        mainPanel(
                          fluidRow(conditionalPanel(condition="input.selextractdata == '1'",
                                                    tabsetPanel(id = "infor_down"  ,       
                                                                tabPanel(strong('Gene Information'), 
                                                                         column(12, DT::dataTableOutput("dMsuInfo"))),
                                                                tabPanel(strong('Keyword'),
                                                                         column(12, DT::dataTableOutput("dMsuKey"))),
                                                                tabPanel(strong('Publication'),
                                                                         column(12,  DT::dataTableOutput("dMsuPub"))),
                                                                tabPanel( strong('RiceNet'),
                                                                          column(12, DT::dataTableOutput("dMsuRiceNet")))
                                                    )),
                                   conditionalPanel(condition="input.selextractdata == '2'",
                                                    tabsetPanel(id = "infor_downor"  ,       
                                                                tabPanel(strong('Gene Information'),
                                                                         column(12,  DT::dataTableOutput("dRapInfo"))
                                                                ),
                                                                tabPanel(strong('Keyword'),
                                                                         column(12,  DT::dataTableOutput("dRapKey"))
                                                                ),
                                                                tabPanel(strong('Publication'),
                                                                         column(12,  DT::dataTableOutput("dRapPub"))
                                                                )
                                                    )
                                   )
                          )
                        )
               ),
               
               ## submit panel
               tabPanel(HTML("<strong style='font-size:18px'>Submit</strong>"),
                        icon = icon("upload", class = NULL, lib = "font-awesome"),
                        fixedRow(
                          column(6,
                                 h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Submit a new Genbank accession</font>')),
                                 
                                 wellPanel(
                                   column(12, textInput('symsub3', strong("Gene symbol"),value="")),
                                   column(12, textInput('accsub3', strong("Accession"),value="")),
                                   column(12, passwordInput('key3', strong("Password"),value="")),
                                   actionButton("submit3", strong("Submit"),style = "color: white;background-color: black"),
                                   actionButton("clear5", strong("Reset"),style = "color: white;background-color: black")
                                 )
                          ),
                          column(6,
                                 h4(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Submit a new connection between genes</font>')),
                                 
                                 wellPanel(
                                   column(12, textInput('symsub5', strong("Gene symbol 1"),value="")),
                                   column(12, textInput('sym2sub5', strong("Gene symbol 2"),value="")),
                                   column(12, textInput('tilsub5', strong("Title"),value="")),
                                   column(12, textInput('evisub5', strong("Evidence"),value="")),
                                   column(12, passwordInput('key5', strong("Password"),value="")),
                                   actionButton("submit5", strong("Submit"),style = "color: white;background-color: black"),
                                   actionButton("clear7", strong("Reset"),style = "color: white;background-color: black")
                                 )
                          )
                        )
               )
    ), right=5, left=5
  )
)
