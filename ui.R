
shinyUI(fluidPage(theme="mystyle.css",
                  tags$head(includeScript("google-analytics.js")),
                  
  fluidRow(
    absolutePanel(
    navbarPage(HTML("<span style='font-family: Comic Sans MS;color:white;'>RicENcode</span>
                    <span style='font-family: Comic Sans MS;color:white; font-size: 50%'>
                       The knowledge of cloned rice genes lost in the 
                       information of rice functional genomic studies
                    </span>
                    "),
      footer=p(HTML("<div align='center'>

                 <a href='http://croplab.hzau.edu.cn' target='_blank'><img src='croplab.png' width='130' height='130'></a>
                 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                 <a href='http://venyao.github.io/RICENCODE/' target='_blank'><img src='rice.png' width='130' height='130'></a>
                 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                 <a href='https://github.com/venyao/RICENCODE' target='_blank'><img src='github.png' width='130' height='130'></a>
                 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                 <a href='https://www.researchgate.net/profile/Wen_Yao' target='_blank'><img src='RG.jpg' width='130' height='130'></a>

              </div>")),
      tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>About</span>"),
      p(HTML("
			 <h3><span style='font-family: Comic Sans MS'>Introdution</span></h3>
             <div align='justify' style='font-family: Comic Sans MS'>
                   More than 2204 cloned rice genes and 303 gene families 
                   comprised of 3381 genes were collected. For each gene, various information 
                   concerning the gene symbol, the genomic locus in the 
                   Nipponbare reference genome, the Genbank accession number, 
                   the title, journal, year, author affiliation, abstract of 
                   related publications and the text mining result of the 
                   publications were provided.
             </div>
			 
			 <h3>How to query this database?</h3>
			 <div align='justify'>
                   This database was designed as a Shiny application and was deployed in the Cloud. 
                   You can query this database <a href='http://ricencode.ncpgr.cn' target='_blank'>HERE</a>.
                   You can query this database using the genomic locus assigned by the <a href='http://rice.plantbiology.msu.edu/' target='_blank'>Rice Genome Annotation Project</a>, 
                   or the genomic locus assigned by the <a href='http://rapdb.dna.affrc.go.jp/' target='_blank'>Rice Annotation Project</a>, or the gene symbol. 
                   You can retrieve the basic information of a gene, the publications on a gene, the agronomic traits associated with a gene, the connections between genes from this database.
                   
				   <br>
                   You can also query this database by downloading the whole database to your local computer. In this way, you need <a href='http://www.rstudio.com/' target='_blank'>RStudio</a> and the <a href='http://shiny.rstudio.com/' target='_blank'>Shiny</a> package installed on your computer (See the <a href='http://venyao.github.io/RICENCODE/assets/ricencode-intro.pdf' target='_blank'>manual</a>).
				   
				   <br>
				   Finally, this database were also provided to users as <a href='http://venyao.github.io/RICENCODE/' target='_blank'>static web pages</a>. You can query this website through in-site search using the searching box on the top right of this website.
             </div>
			 
			 <h3>How to contribute to this database?</h3>
			 <div align='justify'>
			       You can also contribute to this database by submitting information on newly cloned rice gene or new publications to this database. To do this, you need to download the database and have <a href='http://www.rstudio.com/' target='_blank'>RStudio</a> and the <a href='http://shiny.rstudio.com/' target='_blank'>Shiny</a> package installed on your computer(See the <a href='http://venyao.github.io/RICENCODE/assets/ricencode-intro.pdf' target='_blank'>manual</a>).
			 </div>
			 
			 <h3>Further information</h3>
			 <div align='justify'>
			      Further information concerning the details on querying and contributing to this database can be found in the <a href='http://venyao.github.io/RICENCODE/assets/ricencode-intro.pdf' target='_blank'>manual</a>.
			 </div>
			 
			 <h3>Contact us</h3>
			 <div align='justify'>
			      If you encounter any problems or have any suggestions concerning this database, please send email to ywhzau at gmail.com.
			 </div>
			 
             ")),
      br()
	  ),
      
	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Gene</span>"),
	           radioButtons("query", h4("* Query with gene symbol or genomic locus 
                              assigned by the Rice Genome Annotation Project"), 
            choices=c("MSU Locus", "RAPdb Locus", "Gene Symbol"), inline=TRUE),

      uiOutput("inText"),

      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable1"), style = "color: white; background-color: black"),
        tabPanel('Reference', dataTableOutput("mytable2"), style = "color: white; background-color: black"),
        tabPanel('Accession', dataTableOutput("mytable3"), style = "color: white; background-color: black"),
        tabPanel('Expression', dataTableOutput("mytable4"), style = "color: white; background-color: black"),
        tabPanel('Keyword', dataTableOutput("mytable5"), style = "color: white; background-color: black"),
        tabPanel('Connection', dataTableOutput("mytable6"), style = "color: white; background-color: black")
      ),
      
      p(HTML("
             <h4>* Submit new Gene or add new information for existing genes
             </h4>
        ")
      ),
      wellPanel(style = "background-color: #00b271",
        column(3, textInput('symsub7', strong("Gene symbol"),value="")),
        column(3, textInput('msusub7', strong("MSU genomic locus"),value="")),
        column(3, textInput('rapsub7', strong("RAPdb genomic locus"),value="")),
        column(3, textInput('pubmed7', strong("Pubmed ID"),value="")),
        actionButton("submit7", strong("Submit")),
        actionButton("clear1", strong("Clear"))
      ),
      
      p(HTML("
             <h4>* Edit existing gene information
             </h4>
       ")
      ),
      wellPanel(style = "background-color: #b45b3e",
        column(2, textInput('oldsym', strong("Old Gene symbol"),value="")),
        column(2, textInput('newsym', strong("New Gene Symbol"),value="")),
        column(2, textInput('oldmsu', strong("Old MSU locus"),value="")),
        column(2, textInput('newmsu', strong("New MSU locus"),value="")),
        column(2, textInput('oldrap', strong("Old RAPdb locus"),value="")),
        column(2, textInput('newrap', strong("New RAPdb locus"),value="")),
        actionButton("submit6", strong("Submit")),
        actionButton("clear2", strong("Clear"))
      )
	  ),

	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Gene Family</span>"),
	           radioButtons("queryfam", h4("* Query gene family using 
                                    gene symbol or genomic locus"), 
            choices=c("MSU Locus", "RAPdb Locus", "Gene Symbol"), inline=TRUE),

      uiOutput("inTextfam"),

      tabsetPanel(
        tabPanel(strong('Information'), dataTableOutput("mytable8"), style = "color: white; background-color: black"),
        tabPanel('Reference', dataTableOutput("mytable9"), style = "color: white; background-color: black")
      ),
      
      p(HTML("
             <h4>* Submit new gene family
             </h4>
       ")
      ),
      
      wellPanel(style = "background-color: #336699",
        column(6, textInput('pubmed10', strong("Pubmed ID"),value="")),
        #column(4, textInput('famname', strong("Family Name"),value="")),
        column(6, fileInput('genfamin', strong("Gene Family info"),
                            accept=c(".txt"))),
        actionButton("submit10", strong("Submit"))
      )
      
	  ),
      
	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Keyword</span>"),
	  p(HTML("<div align='justify'>
                         Keywords in this database were used to describe phenotypic trait or biological process. As these keywords were collected from the publications, you can only use the keywords listed on <a href='http://venyao.github.io/RICENCODE/tags.html' target='_blank'>this web page</a> to query this database. To use any keyword you like to query this database, you can use the 'Publication' panel on the navigation bar.
             </div>"
	  )),
      textInput("keyword", label = h4("* Query with keyword to characterize agronomic trait"), 
                value = "heading date"),
      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable7"), style = "color: white; background-color: black")
      ),
      
      br(),
      p(HTML("
             <h4>* Submit new keyword to this database
             </h4>
       ")
      ),
      wellPanel(style = "background-color: #336699",
        column(3, textInput('symsub4', strong("Gene symbol"),value="")),
        column(3, textInput('keysub4', strong("Keyword"),value="")),
        column(3, textInput('tilsub4', strong("Title"),value="")),
        column(3, textInput('evisub4', strong("Evidence"),value="")),
        actionButton("submit4", strong("Submit")),
        actionButton("clear3", strong("Clear"))
      )
	  ),
    
    tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Publication</span>"),
      textInput("publication", label = h4("* Query the title or abstract of collected publications"), 
               value = "heading date"),
      tabsetPanel(
        tabPanel('Result', dataTableOutput("mytable11"), style = "color: white; background-color: black")
      ),

      br(),
      
      p(HTML("
             <h4>* Submit new publication to this database
             </h4>
        ")
      ),
      wellPanel(style = "background-color: #8080c0",
        column(2, textInput('symsub2', strong("Gene symbol"),value="")),
        column(2, textInput('tilsub2', strong("Title"),value="")),
        column(2, textInput('yearsub2', strong("Year"),value="")),
        column(2, textInput('jousub2', strong("Journal"),value="")), 
        column(2, textInput('afisub2', strong("Affiliation"),value="")),
        column(2, textInput('abssub2', strong("Abstract"),value="")),
        actionButton("submit2", strong("Submit")),
        actionButton("clear4", strong("Clear"))
      )
	  ),

     tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>ID Conversion</span>"),
              radioButtons("queryconv", h4("* Convert ID of MSU genomic locus
                                  and RAPdb genomic locus"), 
                 choices=c("RAPdb to MSU", "MSU to RAPdb"), inline=TRUE),

      uiOutput("inTextconv"),

      tabsetPanel(
          tabPanel(strong('Result'), dataTableOutput("mytable10"), style = "color: white; background-color: black")
      )
	  ),

	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Submit</span>"),
      
      p(HTML("
             <h4>* Submit new Genbank accession to this database
             </h4>
        ")
      ),
      wellPanel(style = "background-color: #00b271",
        column(6, textInput('symsub3', strong("Gene symbol"),value="")),
        column(6, textInput('accsub3', strong("Accession"),value="")),
        actionButton("submit3", strong("Submit")),
        actionButton("clear5", strong("Clear"))
      ),

      p(HTML("
             <h4>* Submit new Expression information to this database
             </h4>
        ")
      ),
      wellPanel(style = "background-color: #b45b3e",
        column(3, textInput('symsub9', strong("Gene symbol"),value="")),
        column(3, textInput('exp9', strong("Expression"),value="")),
        column(3, textInput('ove9', strong("Overexpression"),value="")),
        column(3, textInput('rnai9', strong("RNAi"),value="")),
        actionButton("submit9", strong("Submit")),
        actionButton("clear6", strong("Clear"))
      ),
      
      p(HTML("
             <h4>* Submit new connection to this database
             </h4>
       ")
      ),
      wellPanel(style = "background-color: #336699",
        column(3, textInput('symsub5', strong("Gene symbol 1"),value="")),
        column(3, textInput('sym2sub5', strong("Gene symbol 2"),value="")),
        column(3, textInput('tilsub5', strong("Title"),value="")),
        column(3, textInput('evisub5', strong("Evidence"),value="")),
        actionButton("submit5", strong("Submit")),
        actionButton("clear7", strong("Clear"))
      ),

      p(HTML("
        <h4>* Submit new phenotype and expression figures to this database
        </h4>
        ")
      ),
      wellPanel(style = "background-color: #8080c0",
        column(4, textInput('symsub8', strong("Gene symbol"),value="")),
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



