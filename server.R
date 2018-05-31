library(shiny)
library(DBI)
library(dplyr)
library(pathview)
library(org.Hs.eg.db)

shinyServer(function(input, output){
  
  
  sumLineList = c("SUM149", "SUM185", "SUM190", "SUM225", "SUM229", "SUM44", 
                  "SUM52", "SUM159", "SUM102", "SUM1315", "MCF10A", "MCF7",
                  "MCF7_LTED")
  
  # Drop-down menu of SUM lines to choose from
  # List could alternatively be pulled from DBIConnection::dbListTables()
  output$sumlineSelect <- renderUI({
    selectInput("sumline", "Cell line:", choices = sumLineList)
  })
  
  # Reactively create enriched KEGG pathway
  output$pathway <- renderImage({
    
    # Prevent possible errors before default inputs load
    if(is.null(input$sumline)){return(NULL)}
    if(is.null(input$pathway)){return(NULL)}
    
    # Connection to Google Cloud SQL server
    SUMLines_DB <-  dbConnect(RMySQL::MySQL(),
                              username = [USERNAME],
                              password = [PASSWORD],
                              host = [IP],
                              port = [PORT],
                              dbname = "SLKBase"
    )
    # Prevent exceeding maximum 16 connections
    on.exit(dbDisconnect(SUMLines_DB))
    
    rs = dbSendQuery(SUMLines_DB, paste0('select * from ', input$sumline))
    df = fetch(rs, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    
    # Take only shRNA screen hits or genes with > 1 expression fold change, depending on setting
    if(input$enrich == "Gene Expression"){
      df.data = df[(df$LogFoldChange > 1 | df$LogFoldChange < -1),c("EntrezID", "LogFoldChange")]
      df.data = df.data[!is.na(df.data$LogFoldChange),]
      if(nrow(df.data) == 0){return(
        list(src = "./www/not_available_message.png",
             contentType = 'image/png',
             width = 722,
             height = 166,
             alt = "loading...")
      )}
      cutoff = ceiling(max(df.data$LogFoldChange))
      rownames(df.data) = as.character(seq(1, nrow(df.data)))
    }else{
      df.data = df[df$ScreenHit == 1, c("EntrezID", "QuantLog")]
      # Removing outliers necessary to scale every SUM line except SUM185
      if(input$sumline == "SUM185"){ cutoff = max(df.data$QuantLog)}
      else{ cutoff = quantile(df.data$QuantLog)[4] + 1.5*IQR(df.data$QuantLog)}
    }
    
    # Create space for color-coding (shRNA screen option only)
    # This color-coding will appropriate pathview's ability to represent
    # multiple datasets/replicates per gene
    df.data$Buffer1 = df.data[,2]
    df.data$Buffer2 = df.data[,2]
    df.data$Buffer3 = df.data[,2]
    df.data$Mutated = df.data[,2]
    df.data$Amplified = df.data[,2]
    df.data$Overexpressed = df.data[,2]
    
    if(input$enrich == "shRNA Screen"){
      # Red band for mutated genes
      df.data$Mutated[df$ScreenHit == 1 && df$OccurencesInCOSMIC > 0] = -1*cutoff
      # Grey band for copy-number amplified genes
      df.data$Amplified[df$ScreenHit == 1 && !is.na(df$DnaAmp)] = 0
      # White band for overexpressed genes
      df.data$Overexpressed[df$ScreenHit == 1 && df$LogFoldChange > 1] = NA
    }
    
    # Remove Entrez IDs and make them row names
    df.mat = as.matrix(df.data[,-1])
    
    rownames(df.mat) = as.character(df.data$EntrezID)
    
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
  
  # Table representing all info for a single input gene
  output$geneInfo <- renderDataTable({
    
    # Table only generated if there is a valid input gene
    if(is.null(input$sumline)){return(NULL)}
    SUMLines_DB <-  dbConnect(RMySQL::MySQL(),
                              username = [USERNAME],
                              password = [PASSWORD],
                              host = [IP],
                              port = [PORT],
                              dbname = "SLKBase"
    )
    on.exit(dbDisconnect(SUMLines_DB))
    rs = dbSendQuery(SUMLines_DB, paste0('select * from ', input$sumline, ' where GeneSymbol = "', input$gene, '"'))
    df = fetch(rs, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    
    if(nrow(df) == 0){return(NULL)}
    df
  }, options = c(paging = F, searching = F))
  
  # Table with info for all shRNA screen hits or genes with > 1 expression fold change in the KEGG pathway, depending on setting
  output$pathwayGenes <- renderDataTable({
    if(is.null(input$sumline)){return(NULL)}
    if(is.null(input$pathway)){return(NULL)}
    SUMLines_DB <-  dbConnect(RMySQL::MySQL(),
                              username = [USERNAME],
                              password = [PASSWORD],
                              host = [IP],
                              port = [PORT],
                              dbname = "SLKBase"
    )
    on.exit(dbDisconnect(SUMLines_DB))
    
    if(input$enrich == "Gene Expression"){
      rs = dbSendQuery(SUMLines_DB, paste0('select * from ', input$sumline, ' where LogFoldChange > 1 OR LogFoldChange < -1'))
    }else{
      rs = dbSendQuery(SUMLines_DB, paste0('select * from ', input$sumline, ' where ScreenHit = 1'))
    }
    
    df = fetch(rs, n=-1)
    dbClearResult(dbListResults(SUMLines_DB)[[1]])
    
    # Get list of all genes in the input KEGG pathway
    keggID = strsplit(input$pathway, '_')[[1]][2]
    kegg = org.Hs.egPATH2EG
    mapped = mappedkeys(kegg)
    kegg2 = as.list(kegg[mapped])
    genes = unlist(kegg2[keggID])
    
    # The above method is easy but not reliable for all pathways
    # The below method gets genes for the rest of the pathways by extracting the gene IDs
    # from the KEGG pathway "KGML" file
    if(is.null(genes)){
      df = df[!is.na(df$EntrezID) & !duplicated(df$EntrezID),]
      rownames(df) = df$EntrezID
      temp = tempfile()
      download.file(paste0("http://rest.kegg.jp/get/hsa", keggID, "/kgml"), destfile = temp)
      node.data = node.info(temp)
      node.type=c("gene","enzyme", "compound", "ortholog")
      sel.idx=node.data$type %in% node.type
      nna.idx=!is.na(node.data$x+node.data$y+node.data$width+node.data$height)
      sel.idx=sel.idx & nna.idx
      node.data=lapply(node.data, "[", sel.idx)
      
      plot.data.gene=node.map(df[,c("QuantLog", "QuantLogRank")], node.data, node.types="gene", 
                              node.sum="sum", entrez.gnodes=FALSE)
      plot.data.gene=plot.data.gene[plot.data.gene$all.mapped != "",]
      
      genes = unlist(strsplit(as.character(plot.data.gene$all.mapped), ","))
    }
    
    df = df[df$EntrezID %in% genes,]
    df[,-1]
  })
  
  
})