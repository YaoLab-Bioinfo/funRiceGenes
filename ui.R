                                                  
library(shinythemes)
library(shinyBS)

shinyUI(fluidPage(theme=shinytheme("darkly"),
                  tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});')), includeScript("google-analytics.js")),
                  
  fluidRow(
    absolutePanel(
    navbarPage(title="funRiceGenes", windowTitle="Welcome to funRiceGenes!",
      footer=p(HTML("<div align='center'>

                 <a href='http://croplab.hzau.edu.cn' target='_blank'><img src='croplab.png' width='130' height='130'></a>
                 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                 <a href='http://funricegenes.github.io' target='_blank'><img src='rice.png' width='130' height='130'></a>

              </div>")),
      
    ## gene panel
	  tabPanel(HTML("<span>Gene</span>"),
	           radioButtons("query", label=h4("* Query with a gene symbol or a genomic locus",
	                                          bsButton("q1", label="", icon=icon("question"), style="info", size="small")), 
            choices=c("MSU Locus", "RAPdb Locus", "Gene Symbol"), inline=TRUE),
            
            bsPopover("q1", "Fill in the table cells (case insensitive) to fetch information on cloned rice genes.",
                      trigger = "focus"),

      uiOutput("inText"),

      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable1")),
        tabPanel('Reference', dataTableOutput("mytable2")),
        tabPanel('Accession', dataTableOutput("mytable3")),
        tabPanel('Expression', dataTableOutput("mytable4")),
        tabPanel('Keyword', dataTableOutput("mytable5")),
        tabPanel('Connection', dataTableOutput("mytable6")),
        tabPanel('RiceNet', dataTableOutput("mytable13"))
      ),
      
      
      h4("* Submit a new Gene or add new information for an existing gene"),
      
      wellPanel(
        column(2, textInput('symsub7', strong("Gene symbol"),value="")),
        column(2, textInput('msusub7', strong("MSU genomic locus"),value="")),
        column(2, textInput('rapsub7', strong("RAPdb genomic locus"),value="")),
        column(2, textInput('pubmed7', strong("Pubmed ID"),value="")),
        column(4, textInput('key7', strong("Password"),value="")),
        actionButton("submit7", strong("Submit")),
        actionButton("clear1", strong("Clear")),
        helpText(HTML("<span style='color:slateblue'>Note: To submit a new publication, fill in the 'Pubmed ID' cell. 
                To submit a new publication for an existing gene, fill in the 'Gene symbol' and 'Pubmed ID' cells. 
                 To submit a new gene, fill in all the four cells.</span>"))
      ),
      
      p(HTML("<table><tr><td>
             <h4>* Edit existing gene information</h4>
			 </td></tr>
	                                      </table>
       ")
      ),
      wellPanel(
        column(2, textInput('oldsym', strong("Old Gene symbol"),value="")),
        column(2, textInput('newsym', strong("New"),value="")),
        column(2, textInput('oldmsu', strong("Old MSU locus"),value="")),
        column(2, textInput('newmsu', strong("New"),value="")),
        column(2, textInput('oldrap', strong("Old RAPdb locus"),value="")),
        column(1, textInput('newrap', strong("New"),value="")),
        column(1, textInput('key6', strong("Password"),value="")),
        actionButton("submit6", strong("Submit")),
        actionButton("clear2", strong("Clear")),
        helpText(HTML("<p style='color:slateblue'>Note: Fill in all the table cells to edit either of the three items: symbol, MSU and RAPdb locus.</p>"))
      )
	  ),

	  ## gene family panel
	  tabPanel(HTML("<span>GeneFamily</span>"),
	           radioButtons("queryfam", h4("* Query with a gene symbol or a genomic locus",
	                                          bsButton("q2", label="", icon=icon("question"), style="info", size="small")),
      choices=c("MSU Locus", "RAPdb Locus", "Gene Symbol"), inline=TRUE),
            
      bsPopover("q2", "Fill in the table cells (case insensitive) to fetch information on rice genes.",
                trigger = "focus"),

      uiOutput("inTextfam"),

      tabsetPanel(
        tabPanel(strong('Information'), dataTableOutput("mytable8")),
        tabPanel('Reference', dataTableOutput("mytable9"))
      ),
      
      p(HTML("<table><tr><td>
             <h4>* Submit a new gene family
             </h4></td></tr>
             </table>
       ")
      ),
      
      wellPanel(
        column(4, textInput('pubmed10', strong("Pubmed ID"),value="")),
        column(4, textInput('key10', strong("Password"),value="")),
        column(4, fileInput('genfamin', strong("Gene Family info"),
                            accept=c(".txt"))),
        actionButton("submit10", strong("Submit")),
        actionButton("clear6", strong("Clear")),
        helpText(HTML("<p style='color:slateblue'>Note: The file 'Gene Family info' should contain 5 columns with names 'Accession',
             'Symbol', 'MSU', 'RAPdb' and 'Name'. Different columns should be separated by '\\t'. And the 1st row should
                 be the column names.</p>"))
      )
	  ),
      
	  ## keyword panel
	  tabPanel(HTML("<span>Keyword</span>"),
      textInput("keyword", label = HTML("<span style='white-space: nowrap'>
                                        <h4>* Query with keywords characterizing agronomic trait of rice listed on 
                                        <a href='https://funricegenes.github.io/tags/' target='_blank'>this page</a></h4></span>"), 
                value = "heading date"),
      
      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable7"))
      ),
      
      p(HTML("<table><tr><td>
             <h4>* Submit new keywords on genes
             </h4></td></tr>
	                                      </table>
       ")
      ),
      wellPanel(
        column(2, textInput('symsub4', strong("Gene symbol"),value="")),
        column(2, textInput('keysub4', strong("Keyword"),value="")),
        column(2, textInput('tilsub4', strong("Title"),value="")),
        column(2, textInput('evisub4', strong("Evidence"),value="")),
        column(4, textInput('key4', strong("Password"),value="")),
        actionButton("submit4", strong("Submit")),
        actionButton("clear3", strong("Clear"))
      )
	  ),
    
	  ## publication panel
    tabPanel(HTML("<span>Publication</span>"),
      textInput("publication", label = HTML("<span style='white-space: nowrap'>
                                            <h4>* Query with any word concerning rice functional genomic studies
                                            </h4></span>"), 	
               value = "heading date"),
      tabsetPanel(
        tabPanel('Result', dataTableOutput("mytable11"))
      ),
      
      p(HTML("<table><tr><td>
             <h4>* Submit a new publication
             </h4></td></tr>
	                                      </table>
        ")
      ),
      wellPanel(
        column(2, textInput('symsub2', strong("Gene symbol"),value="")),
        column(1, textInput('tilsub2', strong("Title"),value="")),
        column(1, textInput('yearsub2', strong("Year"),value="")),
        column(2, textInput('jousub2', strong("Journal"),value="")), 
        column(2, textInput('afisub2', strong("Affiliation"),value="")),
        column(2, textInput('abssub2', strong("Abstract"),value="")),
        column(2, textInput('key2', strong("Password"),value="")),
        actionButton("submit2", strong("Submit")),
        actionButton("clear4", strong("Clear"))
      )
	  ),

	  ## ID conversion panel
     tabPanel(HTML("<span>IDConversion</span>"),
              radioButtons("queryconv", h4("* Convert between MSU genomic locus and RAPdb genomic locus",
                                           bsButton("q3", label="", icon=icon("question"), style="info", size="small")),
                 choices=c("RAPdb to MSU", "MSU to RAPdb"), inline=TRUE),
              
              bsPopover("q3", "You can submit multiple IDs delimited by space.", trigger = "focus"),
      uiOutput("inTextconv"),

      tabsetPanel(
          tabPanel(strong('Result'), dataTableOutput("mytable10"))
      ),
      
      radioButtons("IJconv", h4("* Convert between ", em("indica"), " and ", em("japonica"), " gene IDs",
                                   bsButton("q31", label="", icon=icon("question"), style="info", size="small")),
                   choices=c("MSU Nipponbare", "RAPdb Nipponbare", "Minghui 63", "Zhenshan 97"), inline=TRUE),
      
      bsPopover("q31", "You can submit multiple IDs delimited by space.", trigger = "focus"),
      
      uiOutput("inIJconv"),
      
      tabsetPanel(
        tabPanel(strong('Result'), dataTableOutput("mytable12"))
      )

	  ),

	  ## Bulk download panel
	  tabPanel(HTML("<span>Download</span>"),
	           textAreaInput("msuarea", HTML("<span style='white-space: nowrap'><h4>* Extract data using MSU genomic locus</h4></span>
	                                          "), width="400px", resize="vertical", height="200px", 
	                         placeholder = "One locus in one row"),
	           
	           downloadButton("dMsuInfo", "Download locus information"),
	           downloadButton("dMsuKey", "Download Keywords data"),
	           downloadButton("dMsuPub", "Download literatures"),
	           downloadButton("dMsuRiceNet", "Download RiceNet data"),
	           
	           textAreaInput("raparea", HTML("<span style='white-space: nowrap'><h4>* Extract data using RAP genomic locus</h4></span>
	                                          "), width="400px", resize="vertical", height="200px", 
	                         placeholder = "One locus in one row"),
	           
	           downloadButton("dRapInfo", "Download locus information"),
	           downloadButton("dRapKey", "Download Keywords data"),
	           downloadButton("dRapPub", "Download literatures")
	           
	     ),
	  
	  ## submit panel
	  tabPanel(HTML("<span>Submit</span>"),
      
      p(HTML("<table><tr><td>
             <h4>* Submit a new Genbank accession</h4>
	      </td></tr></table>
        ")
      ),
      wellPanel(
        column(4, textInput('symsub3', strong("Gene symbol"),value="")),
        column(4, textInput('accsub3', strong("Accession"),value="")),
        column(4, textInput('key3', strong("Password"),value="")),
        actionButton("submit3", strong("Submit")),
        actionButton("clear5", strong("Clear"))
      ),
      
      p(HTML("<table><tr><td>
             <h4>* Submit a new connection between genes
             </h4></td></tr>
	                                      </table>
       ")
      ),
      wellPanel(
        column(2, textInput('symsub5', strong("Gene symbol 1"),value="")),
        column(2, textInput('sym2sub5', strong("Gene symbol 2"),value="")),
        column(2, textInput('tilsub5', strong("Title"),value="")),
        column(2, textInput('evisub5', strong("Evidence"),value="")),
        column(4, textInput('key5', strong("Password"),value="")),
        actionButton("submit5", strong("Submit")),
        actionButton("clear7", strong("Clear"))
      ),

      p(HTML("<table><tr><td>
        <h4>* Submit new phenotype and expression figures
        </h4> </td></tr>
	                                      </table>
        ")
      ),
      wellPanel(
        column(2, textInput('symsub8', strong("Gene symbol"),value="")),
        column(2, textInput('key8', strong("Password"),value="")),
        column(4, fileInput('phenofig', strong("Phenotype Figure"),
                            accept=c(".tif", ".png", ".tiff"))),
        column(4, fileInput('expfig', strong("Expression Figure"),
                            accept=c(".tif", ".png", ".tiff"))),
        actionButton("submit8", strong("Submit")),
        actionButton("clear8", strong("Clear"))
      )
      
	  )
    ),
  right=5, left=5
  )
  ) 
))

