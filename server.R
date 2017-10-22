library(shiny)
library(DBI)
library(dplyr)
library(pathview)


shinyServer(function(input, output){
  
  # SUMLines_DB <-  dbConnect(RMySQL::MySQL(),
  #                           username = "root",
  #                           password = "sumlines",
  #                           host = "35.190.141.134",
  #                           port = 3306,
  #                           dbname = "SUMLines_DB"
  # )
  # tables = dbListTables(SUMLines_DB)
  # tables = tables[2:length(tables)]
  # dbDisconnect(SUMLines_DB)
  # sumLineList = c()
  # i = 1
  # while(i < length(tables)){
  #   sumLineList = append(sumLineList, strsplit(tables[i], '_')[[1]][1])
  #   i = i + 3
  # }
  sumLineList = c("SUM149", "SUM185", "SUM190", "SUM225", "SUM229", "SUM44", "SUM52")
  
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
    
    rs.fc = dbSendQuery(SUMLines_DB, paste0('select Symbol,', input$sumline, ' from ExpressionDataSUMLinesAllGenes'))
    df.fc = fetch(rs.fc, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    
    
    if(input$enrich == "Gene Expression"){
      # rs.fc = dbSendQuery(SUMLines_DB, paste0('select Symbol,', input$sumline, ' from ExpressionDataSUMLinesAllGenes'))
      # df.fc = fetch(rs.fc, n=-1)
      # dbClearResult(dbListResults(SUMLines_DB)[[1]])
      df.data = df.fc[(df.fc[,2] > 1 | df.fc[,2] < -1),]
      cutoff = ceiling(max(df.data[,2]))
    }else{
      rs.mut = dbSendQuery(SUMLines_DB, paste0('select gene from ', input$sumline, '_Mut_COSMIC'))
      df.mut = fetch(rs.mut, n=-1)
      dbClearResult(dbListResults(SUMLines_DB)[[1]])
      rs.amp = dbSendQuery(SUMLines_DB, paste0('select Symbol from ', input$sumline, '_amplificationData'))
      df.amp = fetch(rs.amp, n=-1)
      dbClearResult(dbListResults(SUMLines_DB)[[1]])
      
      df.overexp = df.fc[df.fc[,2] > 1, 1]
      rs.data = dbSendQuery(SUMLines_DB, paste0('select ids,quantlog from ', input$sumline, '_CellectaData where quantloghit = 1'))
      df.data = fetch(rs.data, n=-1)
      dbClearResult(dbListResults(SUMLines_DB)[[1]])
      if(input$sumline == "SUM185"){ cutoff = max(df.data[,2])}
      else{ cutoff = quantile(df.data[,2])[4] + 1.5*IQR(df.data[,2])}
    }
    
    df.data$Buffer1 = df.data[,2]
    df.data$Buffer2 = df.data[,2]
    df.data$Buffer3 = df.data[,2]
    df.data$Mutated = df.data[,2]
    df.data$Amplified = df.data[,2]
    df.data$Overexpressed = df.data[,2]
    
    if(input$enrich == "shRNA Screen"){
    df.data$Mutated[df.data[,1] %in% df.mut[,1]] = -1*cutoff
    df.data$Amplified[df.data[,1] %in% df.amp[,1]] = 0
    df.data$Overexpressed[df.data[,1] %in% df.overexp] = NA
    }
    
    df.mat = as.matrix(df.data[,-1])
    
    rownames(df.mat) = id2eg(df.data[,1])[,2]
    
    keggID = strsplit(input$pathway, '_')[[1]][2]
    outfile <- paste0("./hsa", keggID, ".pathview.multi.png")
    pathview(gene.data = df.mat, pathway.id = keggID, species = "hsa", 
             limit=list(gene = cutoff, cpd = 1),  kegg.dir = tempdir())
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
    rs.amp = dbSendQuery(SUMLines_DB, paste0('select `DNA Amp.` from ', input$sumline, '_amplificationData where Symbol = "', input$gene, '"'))
    df.amp = fetch(rs.amp, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    rs.fc = dbSendQuery(SUMLines_DB, paste0('select ', input$sumline, ' from ExpressionDataSUMLinesAllGenes where Symbol = "', input$gene, '"'))
    df.fc = fetch(rs.fc, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    rs.mut = dbSendQuery(SUMLines_DB, paste0('select geneMutation,Occurences_in_COSMIC from ', input$sumline, '_Mut_COSMIC where gene = "', input$gene, '"'))
    df.mut = fetch(rs.mut, n=-1)
    if(nrow(df.cell) == 0){return(NULL)}
    df = df.cell
    if(nrow(df.amp) == 0){
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
    if(nrow(df.fc) == 0){
      df$log_Fold_Change = c("NA")
    }else{
      df$log_Fold_Change = df.fc[,1]
    }
    df
  }, options = c(paging = F, searching = F))
  
})