
shinyUI(fluidPage(theme="mystyle.css",
                  tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});')), includeScript("google-analytics.js")),
                  
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
                 <a href='http://ricencode.github.io' target='_blank'><img src='rice.png' width='130' height='130'></a>
                 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
                 <a href='https://github.com/venyao/RICENCODE' target='_blank'><img src='github.png' width='130' height='130'></a>

              </div>")),
      
	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Gene</span>"),
	           radioButtons("query", HTML("<table style='background-color: #DDF3FF'><tr><td><h4>* Query with gene symbol or genomic locus</h4></td>
<td>
<div class='help-tip'>
	<p>Fill in the table cells (case insensitive) to fetch information on cloned rice genes.</p>
	                                      </div></td></tr>
	                                      </table>
	                                      "), 
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
      
      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
             <h4>* Submit new Gene or add new information for existing genes</h4>
             </td><td>
<div class='help-tip'>
	<p>To only submit a new publication, fill in the 'Pubmed ID' cell. 
     To submit a new publication for an existing gene, fill in the 'Gene symbol' and 'Pubmed ID' cells.
     To submit a new gene, fill in all the four cells.
</p>
	                                      </div></td></tr>
	                                      </table>
        ")
      ),
      wellPanel(style = "background-color: #00b271",
        column(2, textInput('symsub7', strong("Gene symbol"),value="")),
        column(2, textInput('msusub7', strong("MSU genomic locus"),value="")),
        column(2, textInput('rapsub7', strong("RAPdb genomic locus"),value="")),
        column(2, textInput('pubmed7', strong("Pubmed ID"),value="")),
        column(4, textInput('key7', strong("Password"),value="")),
        actionButton("submit7", strong("Submit")),
        actionButton("clear1", strong("Clear"))
      ),
      
      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
             <h4>* Edit existing gene information</h4>
			 </td><td>
<div class='help-tip'>
	<p>Fill in all the table cells to edit either of the three items: symbol, MSU and RAPdb locus.
</p>
	                                      </div></td></tr>
	                                      </table>
       ")
      ),
      wellPanel(style = "background-color: #b45b3e",
        column(2, textInput('oldsym', strong("Old Gene symbol"),value="")),
        column(2, textInput('newsym', strong("New"),value="")),
        column(2, textInput('oldmsu', strong("Old MSU locus"),value="")),
        column(2, textInput('newmsu', strong("New"),value="")),
        column(2, textInput('oldrap', strong("Old RAPdb locus"),value="")),
        column(1, textInput('newrap', strong("New"),value="")),
        column(1, textInput('key6', strong("Password"),value="")),
        actionButton("submit6", strong("Submit")),
        actionButton("clear2", strong("Clear"))
      )
	  ),

	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>GeneFamily</span>"),
	           radioButtons("queryfam", HTML("<table style='background-color: #DDF3FF'><tr><td><h4>* Query with gene symbol or genomic locus</h4></td>
<td>
<div class='help-tip'>
	<p>Fill in the table cells (case insensitive) to fetch information on rice genes.</p>
	                                      </div></td></tr>
	                                      </table>
	                                      "),
            choices=c("MSU Locus", "RAPdb Locus", "Gene Symbol"), inline=TRUE),

      uiOutput("inTextfam"),

      tabsetPanel(
        tabPanel(strong('Information'), dataTableOutput("mytable8"), style = "color: white; background-color: black"),
        tabPanel('Reference', dataTableOutput("mytable9"), style = "color: white; background-color: black")
      ),
      
      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
             <h4>* Submit new gene family
             </h4></td><td>
             <div class='help-tip'>
             <p>The file 'Gene Family info' should contain 5 columns with names 'Accession',
             'Symbol', 'MSU', 'RAPdb' and 'Name'. Different columns should be separated by '\\t'. And the 1st row should
             be the column names.
             </p>
             </div></td></tr>
             </table>
       ")
      ),
      
      wellPanel(style = "background-color: #336699",
        column(4, textInput('pubmed10', strong("Pubmed ID"),value="")),
        column(4, textInput('key10', strong("Password"),value="")),
        column(4, fileInput('genfamin', strong("Gene Family info"),
                            accept=c(".txt"))),
        actionButton("submit10", strong("Submit")),
        actionButton("clear6", strong("Clear"))
      )
      
	  ),
      
	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Keyword</span>"),
      textInput("keyword", label = HTML("<table style='background-color: #DDF3FF'><tr><td>
                                        <h4>* Query with keyword characterizing agronomic trait of rice</h4>
                                        </td><td>
<div class='help-tip'>
	<p>Keywords in this database were used to describe phenotypic trait or biological process. 
As these keywords were collected from publications, 
you can only use the keywords listed on <a href='http://ricencode.github.io/tags/' target='_blank'>
this web page</a> to query this database. 
To use any keyword you like to query this database, 
you can use the 'Publication' panel on the navigation bar.
</p>
	                                      </div></td></tr>
	                                      </table>
                                        "), 
                value = "heading date"),
      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable7"), style = "color: white; background-color: black")
      ),
      
      br(),
      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
             <h4>* Submit new keywords on genes
             </h4></td><td>
<div class='help-tip'>
	<p>Key information on genes extracted from abstracts.
</p>
	                                      </div></td></tr>
	                                      </table>
       ")
      ),
      wellPanel(style = "background-color: #336699",
        column(2, textInput('symsub4', strong("Gene symbol"),value="")),
        column(2, textInput('keysub4', strong("Keyword"),value="")),
        column(2, textInput('tilsub4', strong("Title"),value="")),
        column(2, textInput('evisub4', strong("Evidence"),value="")),
        column(4, textInput('key4', strong("Password"),value="")),
        actionButton("submit4", strong("Submit")),
        actionButton("clear3", strong("Clear"))
      )
	  ),
    
    tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Publication</span>"),
      textInput("publication", label = HTML("<table style='background-color: #DDF3FF'><tr><td>
                                        <h4>* Query with any word concerning rice functional genomic studies</h4>
                                        </td><td>
<div class='help-tip'>
	<p>Any words you can think of.
</p>
	                                      </div></td></tr>
	                                      </table>
                                        "), 	
               value = "heading date"),
      tabsetPanel(
        tabPanel('Result', dataTableOutput("mytable11"), style = "color: white; background-color: black")
      ),

      br(),
      
      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
             <h4>* Submit new publication
             </h4></td><td>
<div class='help-tip'>
	<p>Submit publications not archived in Pubmed.
</p>
	                                      </div></td></tr>
	                                      </table>
        ")
      ),
      wellPanel(style = "background-color: #8080c0",
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

     tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>IDConversion</span>"),
              radioButtons("queryconv", HTML("<table style='background-color: #DDF3FF'><tr><td><h4>* Convert between MSU genomic locus
                                  and RAPdb genomic locus</h4></td>
<td>
<div class='help-tip'>
	<p>Conversion.</p>
	                                      </div></td></tr>
	                                      </table>
	                                      "),
                 choices=c("RAPdb to MSU", "MSU to RAPdb"), inline=TRUE),

      uiOutput("inTextconv"),

      tabsetPanel(
          tabPanel(strong('Result'), dataTableOutput("mytable10"), style = "color: white; background-color: black")
      )
	  ),

	  tabPanel(HTML("<span style='font-family: Comic Sans MS;color:white'>Submit</span>"),
      
      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
             <h4>* Submit new Genbank accession</h4>
	      </td><td>
<div class='help-tip'>
	<p>One item at a time.
</p>
	                                      </div></td></tr>
	                                      </table>
        ")
      ),
      wellPanel(style = "background-color: #00b271",
        column(4, textInput('symsub3', strong("Gene symbol"),value="")),
        column(4, textInput('accsub3', strong("Accession"),value="")),
        column(4, textInput('key3', strong("Password"),value="")),
        actionButton("submit3", strong("Submit")),
        actionButton("clear5", strong("Clear"))
      ),
      
      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
             <h4>* Submit new connection between genes
             </h4></td><td>
<div class='help-tip'>
	<p>One item at a time.
</p>
	                                      </div></td></tr>
	                                      </table>
       ")
      ),
      wellPanel(style = "background-color: #336699",
        column(2, textInput('symsub5', strong("Gene symbol 1"),value="")),
        column(2, textInput('sym2sub5', strong("Gene symbol 2"),value="")),
        column(2, textInput('tilsub5', strong("Title"),value="")),
        column(2, textInput('evisub5', strong("Evidence"),value="")),
        column(4, textInput('key5', strong("Password"),value="")),
        actionButton("submit5", strong("Submit")),
        actionButton("clear7", strong("Clear"))
      ),

      p(HTML("<table style='background-color: #DDF3FF'><tr><td>
        <h4>* Submit new phenotype and expression figures
        </h4> </td><td>
<div class='help-tip'>
	<p>One item at a time.
</p>
	                                      </div></td></tr>
	                                      </table>
        ")
      ),
      wellPanel(style = "background-color: #8080c0",
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



