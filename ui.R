
shinyUI(fluidPage(
  
  fluidRow(
    absolutePanel(
      br(),
      p(HTML("<div align='center'; style='font-size:500%'><font face='Monotype Corsiva'; color='green3';>RicENcode</font></div>")),
#      p(HTML("<div align='center'><img src='E:/GIT/RICENCODE/ricencode.jpg'></div>")),
      br(),
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      textInput("msu", label = h4("MSU genomic locus:"), 
                value = "LOC_Os07g15770"),
      tabsetPanel(
        tabPanel(strong('Information'), dataTableOutput("mytable1")),
        tabPanel('Reference', dataTableOutput("mytable2")),
        tabPanel('Accession', dataTableOutput("mytable3")),
        tabPanel('Text-mining', dataTableOutput("mytable4")),
        tabPanel('Keyword', dataTableOutput("mytable5")),
        tabPanel('Connection', dataTableOutput("mytable6"))
      ),
      
      br(),
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      textInput("rap", label = h4("RAPdb genomic locus:"), 
                value = "Os05g0158500"),
      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable7")),
        tabPanel('Reference', dataTableOutput("mytable8")),
        tabPanel('Accession', dataTableOutput("mytable9")),
        tabPanel('Text-mining', dataTableOutput("mytable10")),
        tabPanel('Keyword', dataTableOutput("mytable11")),
        tabPanel('Connection', dataTableOutput("mytable12"))
      ), 
      
      br(),
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      textInput("symbol", label = h4("Gene symbol:"), 
                value = "Moc1"),
      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable13")),
        tabPanel('Reference', dataTableOutput("mytable14")),
        tabPanel('Accession', dataTableOutput("mytable15")),
        tabPanel('Text-mining', dataTableOutput("mytable16")),
        tabPanel('Keyword', dataTableOutput("mytable17")),
        tabPanel('Connection', dataTableOutput("mytable18"))
      ),
      
      br(),
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      textInput("keyword", label = h4("Keyword:"), 
                value = "heading date"),
      tabsetPanel(
        tabPanel('Information', dataTableOutput("mytable19"))
      ),
      
      br(),
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      helpText(h5("Submit new gene.")),
      wellPanel(
        column(4, textInput('symsub1', strong("Gene symbol"),value="")),
        column(4, textInput('msusub', strong("MSU genomic locus"),value="")),
        column(4, textInput('rapsub', strong("RAPdb genomic locus"),value="")),
        actionButton("submit1", strong("Submit"))
      ),
      
      textOutput("mytext20"),
      
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      helpText(h5("Submit new publication.")),
      wellPanel(
        column(2, textInput('symsub2', strong("Gene symbol"),value="")),
        column(2, textInput('tilsub', strong("Title"),value="")),
        column(2, textInput('yearsub', strong("Year"),value="")),
        column(2, textInput('jousub', strong("Journal"),value="")), 
        column(2, textInput('afisub', strong("Affiliation"),value="")),
        column(2, textInput('abssub', strong("Abstract"),value="")),
        actionButton("submit2", strong("Submit"))
      ),
      
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      helpText(h5("Submit new Genbank accession.")),
      wellPanel(
        column(6, textInput('symsub3', strong("Gene symbol"),value="")),
        column(6, textInput('accsub', strong("Accession"),value="")),
        actionButton("submit3", strong("Submit"))
      ),
      
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      helpText(h5("Submit new Keyword.")),
      wellPanel(
        column(3, textInput('symsub4', strong("Gene symbol"),value="")),
        column(3, textInput('keysub', strong("Keyword"),value="")),
        column(3, textInput('tilsub4', strong("Title"),value="")),
        column(3, textInput('evisub4', strong("Evidence"),value="")),
        actionButton("submit4", strong("Submit"))
      ),
      
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
                       blue;'></div></b>")),
      
      helpText(h5("Submit new Connection.")),
      wellPanel(
        column(3, textInput('symsub5', strong("Gene symbol 1"),value="")),
        column(3, textInput('sym2sub5', strong("Gene symbol 2"),value="")),
        column(3, textInput('tilsub5', strong("Title"),value="")),
        column(3, textInput('evisub5', strong("Evidence"),value="")),
        actionButton("submit5", strong("Submit"))
      ),
      
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
             blue;'></div></b>")),
      
      helpText(h4("Update the database.")),
      wellPanel(
        actionButton("submit6", strong("Submit"))
      ),
      
      p(HTML("<b><div style='background-color:#FADDF2;border:1px solid
             blue;'></div></b>")),
      
      p(HTML("<div align='center'>
                 <a href='http://www.croplab.org'>National Key Laboratory of Crop Genetic Improvement</a>
                 <br>
                 <a href='http://www.ncpgr.cn'>National Center of Plant Gene Research (Wuhan)</a>
              </div>")),
      br(),
      
      right=5, left=10
    )
  )
  
))



