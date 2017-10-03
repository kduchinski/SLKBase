library(shiny)
library(DBI)
library(dplyr)
library(pathview)


shinyServer(function(input, output){
  
  SUMLines_DB <-  dbConnect(RMySQL::MySQL(),
                            username = [USERNAME],
                            password = [PASSWORD],
                            host = [INSTANCE IP],
                            port = [PORT],
                            dbname = "SUMLines_DB"
  )
  tables = dbListTables(SUMLines_DB)
  dbDisconnect(SUMLines_DB)
  sumLineList = c()
  i = 1
  while(i < length(tables)){
    sumLineList = append(sumLineList, strsplit(tables[i], '_')[[1]][1])
    i = i + 3
  }

  output$sumlineSelect <- renderUI({
    selectInput("sumline", "Cell line:", choices = sumLineList, selected = sumLineList[1])
  })
  output$pathway <- renderImage({
    if(is.null(input$sumline)){return(NULL)}
    if(is.null(input$pathway)){return(NULL)}
    SUMLines_DB <-  dbConnect(RMySQL::MySQL(),
                              username = "root",
                              password = "sumlines",
                              host = "35.190.141.134",
                              port = 3306,
                              dbname = "SUMLines_DB"
    )
    on.exit(dbDisconnect(SUMLines_DB))
    rs.data = dbSendQuery(SUMLines_DB, paste0('select ids,quantlog from ', input$sumline, '_CellectaData'))
    df.data = fetch(rs.data, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    rs.mut = dbSendQuery(SUMLines_DB, paste0('select * from ', input$sumline, '_Mut_COSMIC'))
    df.mut = fetch(rs.mut, n=-1)
    df.data$COSMIC = df.data[,2]
    df.data$COSMIC[df.data[,1] %in% df.mut[,2]] = NA
    
    head(df.data)
    df.mat = as.matrix(df.data[,2:3]/max(df.data[,2]))
    
    rownames(df.mat) = id2eg(df.data[,1])[,2]
    colnames(df.mat) = c("QuantLog", "COSMIC")

    keggID = strsplit(input$pathway, '_')[[1]][2]
    outfile <- paste0("./hsa", keggID, ".pathview.multi.png")
    pathview(gene.data = df.mat[,1:2], pathway.id = keggID, species = "hsa", 
             limit=list(gene = 1, cpd = 1),  kegg.dir = tempdir())
    list(src = outfile,
         contentType = 'image/png',
         width = 1000,
         height = 1200,
         alt = "loading...")
  })
  output$geneInfo <- renderDataTable({
    if(is.null(input$sumline)){return(NULL)}
    SUMLines_DB <-  dbConnect(RMySQL::MySQL(),
                              username = "root",
                              password = "sumlines",
                              host = "35.190.141.134",
                              port = 3306,
                              dbname = "SUMLines_DB"
    )
    on.exit(dbDisconnect(SUMLines_DB))
    rs.cell = dbSendQuery(SUMLines_DB, paste0('select ids,quantlog,quantlogrank from ', input$sumline, '_CellectaData where ids = "', input$gene, '"'))
    df.cell = fetch(rs.cell, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    rs.amp = dbSendQuery(SUMLines_DB, paste0('select `log Fold Change`,`DNA Amp.` from ', input$sumline, '_amplificationData where Symbol = "', input$gene, '"'))
    df.amp = fetch(rs.amp, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    rs.mut = dbSendQuery(SUMLines_DB, paste0('select geneMutation,Occurences_in_COSMIC from ', input$sumline, '_Mut_COSMIC where gene = "', input$gene, '"'))
    df.mut = fetch(rs.mut, n=-1)
    if(nrow(df.cell) == 0){return(NULL)}
    df = df.cell
    if(nrow(df.amp) == 0){
      df$log_Fold_Change = c("NA")
      df$DNA_Amp = c("NA")
    }else{
      df = cbind(df.cell, df.amp)
    }
    if(nrow(df.mut) == 0){
      df$geneMutation = c("NA")
      df$Occurences_in_COSMIC = c("NA")
    }else{
      df = cbind(df,df.mut)
    }
    df
  }, options = c(paging = F, searching = F))
  
})
