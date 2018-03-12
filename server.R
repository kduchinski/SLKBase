library(shiny)
library(DBI)
library(dplyr)
library(pathview)
library(org.Hs.eg.db)

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
  sumLineList = c("SUM149", "SUM185", "SUM190", "SUM225", "SUM229", "SUM44", "SUM52", "SUM159")
  
  output$sumlineSelect <- renderUI({
    selectInput("sumline", "Cell line:", choices = sumLineList)
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
    
    rs.fc = dbSendQuery(SUMLines_DB, paste0('select Symbol,', input$sumline, ',EntrezId from ExpressionDataSUMLinesAllGenes'))
    df.fc = fetch(rs.fc, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    
    
    if(input$enrich == "Gene Expression"){
      df.data = df.fc[(df.fc[,2] > 1 | df.fc[,2] < -1),]
      if(nrow(df.data) == 0){return(
        list(src = "./www/not_available_message.png",
             contentType = 'image/png',
             width = 722,
             height = 166,
             alt = "loading...")
      )}
      cutoff = ceiling(max(df.data[,2]))
      rownames(df.data) = as.character(seq(1, nrow(df.data)))
    }else{
      rs.mut = dbSendQuery(SUMLines_DB, paste0('select gene from ', input$sumline, '_Mut_COSMIC'))
      df.mut = fetch(rs.mut, n=-1)
      dbClearResult(dbListResults(SUMLines_DB)[[1]])
      rs.amp = dbSendQuery(SUMLines_DB, paste0('select Symbol from ', input$sumline, '_amplificationData'))
      df.amp = fetch(rs.amp, n=-1)
      dbClearResult(dbListResults(SUMLines_DB)[[1]])
      
      df.overexp = df.fc[df.fc[,2] > 1, 1]
      rs.data = dbSendQuery(SUMLines_DB, paste0('select ', input$sumline, '_CellectaData.ids,', input$sumline, 
                                                '_CellectaData.quantlog,ExpressionDataSUMLinesAllGenes.EntrezId from ', input$sumline, 
                                                '_CellectaData left join ExpressionDataSUMLinesAllGenes on ', input$sumline, 
                                                '_CellectaData.ids = ExpressionDataSUMLinesAllGenes.Symbol where quantloghit = 1'))
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
      
      rows.noID = as.numeric(rownames(df.data[is.na(df.data$EntrezId) & df.data[,1] != "",]))
      try.match = pathview::id2eg(df.data[rows.noID,1])[,2]
      df.data[rows.noID,3] = try.match
    }
    
    df.mat = as.matrix(df.data[,-c(1,3)])
    
    
    rownames(df.mat) = as.character(df.data[,3])
    
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
  
  output$pathwayGenes <- renderDataTable({
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
    
    if(input$enrich == "Gene Expression"){
      rs = dbSendQuery(SUMLines_DB, paste0('select ids,EntrezId,quantlog,quantlogrank,`DNA Amp.`,', input$sumline,
                                           ',geneMutation,Occurences_in_COSMIC from ', input$sumline,
                                           '_CellectaData left join ExpressionDataSUMLinesAllGenes on ', input$sumline,
                                           '_CellectaData.ids = ExpressionDataSUMLinesAllGenes.Symbol left join ', input$sumline,
                                           '_amplificationData on ids = ', input$sumline, '_amplificationData.Symbol left join ',
                                           input$sumline, '_Mut_COSMIC on ', input$sumline, '_CellectaData.ids = ', input$sumline,
                                           '_Mut_COSMIC.gene where (', input$sumline, ' > 1 OR ', input$sumline, ' < -1)'))
    }else{
      rs = dbSendQuery(SUMLines_DB, paste0('select ids,EntrezId,quantlog,quantlogrank,`DNA Amp.`,', input$sumline,
                                           ',geneMutation,Occurences_in_COSMIC from ', input$sumline,
                                           '_CellectaData left join ExpressionDataSUMLinesAllGenes on ', input$sumline,
                                           '_CellectaData.ids = ExpressionDataSUMLinesAllGenes.Symbol left join ', input$sumline,
                                           '_amplificationData on ids = ', input$sumline, '_amplificationData.Symbol left join ',
                                           input$sumline, '_Mut_COSMIC on ', input$sumline, '_CellectaData.ids = ', input$sumline,
                                           '_Mut_COSMIC.gene where quantloghit = 1'))
    }
    
    df = fetch(rs, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    df.unmapped = df[is.na(df$EntrezId),]
    colnames(df)[6] = "foldChange"
    
    rows.noID = as.numeric(rownames(df[is.na(df$EntrezId),]))
    try.match = id2eg(df[rows.noID,1])[,2]
    df[rows.noID,2] = try.match
    
    keggID = strsplit(input$pathway, '_')[[1]][2]
    kegg = org.Hs.egPATH2EG
    mapped = mappedkeys(kegg)
    kegg2 = as.list(kegg[mapped])
    genes = unlist(kegg2[keggID])
    if(is.null(genes)){
      df = df[!is.na(df$EntrezId) & !duplicated(df$EntrezId),]
      rownames(df) = df$EntrezId
      temp = tempfile()
      download.file(paste0("http://rest.kegg.jp/get/hsa", keggID, "/kgml"), destfile = temp)
      node.data = node.info(temp)
      node.type=c("gene","enzyme", "compound", "ortholog")
      sel.idx=node.data$type %in% node.type
      nna.idx=!is.na(node.data$x+node.data$y+node.data$width+node.data$height)
      sel.idx=sel.idx & nna.idx
      node.data=lapply(node.data, "[", sel.idx)
      
      plot.data.gene=node.map(df[,3:4], node.data, node.types="gene", node.sum="sum", entrez.gnodes=FALSE)
      plot.data.gene=plot.data.gene[plot.data.gene$all.mapped != "",]
      
      genes = unlist(strsplit(as.character(plot.data.gene$all.mapped), ","))
    }

    df = df[df$EntrezId %in% genes,]
    df[,-2]
    #if(nrow(df) == 0){return(NULL)}
  })
  
  
})