server <- function(input, output) {

options(shiny.usecairo=T)
  
  plotInput <- function(){
    MSET <- input$MSET
    PIE.TYPE <- input$PieType
    QVAL <- as.numeric(input$qval)/100
    FDR <- as.numeric(input$fdr)/100
    LFC <- as.numeric(input$fc)
    res <- tmodLimmaTest(efit, genes = "SYMBOL", sort.by = "msd", tmodFunc = tmodCERNOtest, mset = MSET)
    pie <- tmodLimmaDecideTests(efit, lfc.thr = LFC, pval.thr = FDR, genes="SYMBOL")
    tmodPanelPlot(res, pie=pie, pie.style = PIE.TYPE, filter.rows.pval = QVAL, pval.thr = QVAL,
                  filter.empty.cols =T, filter.empty.rows = T, legend.style = "auto",
                  col.labels.style = "bottom")
  }
  
  output$distPlot <- renderPlot({
      plotInput()
  }, height=1400, width=1000)
    
  output$export <- downloadHandler(
    filename = function() { paste("Tmod plot Moduleset-", input$MSET, " Plot-",input$PieType, " FDR-", input$fdr, " LFC-", input$fc, ".pdf", sep="") },
    content = function(file) {
      # ggsave(file, plot = grid.draw(plotInput()), device = "pdf", width=20,height = 20)
      cairo_pdf(filename = file,
                width = 18, height = 20)
      plotInput()
      dev.off()
  })
  
  plotInput1 <- function(){
    plot(cmd1, col=as.character(module.colours), main="MDS plot",
         xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
  }
  plotInput2 <- function(){
    plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05, main="Consensus gene dendrogram and module colors")
  }
  plotInput3 <- function(){
    par(mar = c(9, 9, 3, 3))
    labeledHeatmap(Matrix = modTraitCor, xLabels = names(datvar),
                   yLabels = names(MEs), ySymbols = names(MEs),
                   colorLabels =FALSE,colors=blueWhiteRed(100),textMatrix=textMatrix,
                   setStdMargins = FALSE, zlim = c(-1,1),
                   #cex.text = 1.25,cex.lab.y = 1.15,cex.lab.x = 1.75,
                   main = paste("Module-trait relationships"),xLabelsAngle = 45) 
  }
  plotInput4 <- function(){
    midpoint <-0
    base_size <- 20
    ggplot(df.heatmap, aes(Comparison, pathway)) +
      geom_tile(aes(fill = NES), colour = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = midpoint, na.value = "gray50") +
      theme_grey(base_size = base_size) + labs(x = "",  y = "") +
      scale_x_discrete(expand = c(0, 0), position = "top") +
      scale_y_discrete(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 35, hjust = 0, vjust = 1), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank(), panel.border = element_blank(), panel.background = element_rect(fill = "white",color = 'white'), legend.key.size = unit(1, "cm"), legend.text=element_text(size=rel(1)))
  }
  
  output$distPlot1 <- renderPlot({
    if (input$MES=="MDS") {plotInput1()}
    if (input$MES=="MDRO") {plotInput2()}
    if (input$MES=="MT") {plotInput3()}
    if (input$MES=="BTM") {plotInput4()}
  }, height=800, width=800)
  
  output$export1 <- downloadHandler(
    filename = function() { paste("WGCNA-", input$MES, " Plot", ".pdf", sep="") },
    content = function(file) {
      # ggsave(file, plot = grid.draw(plotInput()), device = "pdf", width=20,height = 20)
      cairo_pdf(filename = file,
                width = 25, height = 20)
      if (input$MES=="MDS") {plotInput1()}
      if (input$MES=="MDRO") {plotInput2()}
      if (input$MES=="MT") {plotInput3()}
      if (input$MES=="BTM") {plotInput4()}
      dev.off()
    })
  
  datasetInput1 <- reactive({
    DEGFDR <- as.numeric(input$degfdr)/100
    DEGLFC <- as.numeric(input$degfc)
    if(input$genesymbolsonly == TRUE) { 
      alltabs <- alltabs[alltabs$SYMBOL!= "",]}
    if(input$duplicated == TRUE) {
      alltabs <- alltabs}
    if(input$duplicated == TRUE & input$genesymbolsonly == TRUE) {
      alltabs <- alltabs[alltabs$SYMBOL!= "",]}
    alltabs[alltabs$adj.P.Val < DEGFDR & abs(alltabs$logFC) > DEGLFC,] 
  })
  
  output$exportTable1 <- downloadHandler(
    filename = function() { paste0("All DEGs ", "FDR-", input$degfdr, " LFC-", input$degfc, ".csv", sep="") },
    content = function(file) {
      write.csv(datasetInput1(),file, row.names = F)
    }
  )
  
  output$table2 <- DT::renderDataTable({
    DEGFDR <- as.numeric(input$degfdr)/100
    DEGLFC <- as.numeric(input$degfc)
    if(input$genesymbolsonly == TRUE) { 
       tableS2 <- tableS2[tableS2$SYMBOL!= "",]}
    if(input$duplicated == TRUE) { 
      tableS2 <- tableS2[!duplicated(tableS2$SYMBOL),]}
    if(input$duplicated == TRUE & input$genesymbolsonly == TRUE) { 
      tableS2 <- tableS2[tableS2$SYMBOL!= "" & !duplicated(tableS2$SYMBOL),]}
    tableS2[tableS2$adj.P.Val < DEGFDR & abs(tableS2$logFC) > DEGLFC,]
  })
  
  datasetInput2 <- reactive({
    DEGFDR <- as.numeric(input$degfdr)/100
    DEGLFC <- as.numeric(input$degfc)
    if(input$genesymbolsonly == TRUE) { 
      tableS2 <- tableS2[tableS2$SYMBOL!= "",]}
    if(input$duplicated == TRUE) { 
      tableS2 <- tableS2[!duplicated(tableS2$SYMBOL),]}
    if(input$duplicated == TRUE & input$genesymbolsonly == TRUE) { 
      tableS2 <- tableS2[tableS2$SYMBOL!= "" & !duplicated(tableS2$SYMBOL),]}
    tableS2[tableS2$adj.P.Val < DEGFDR & abs(tableS2$logFC) > DEGLFC,]
  })
  
  output$exportTable2 <- downloadHandler(
    filename = function() { paste0("CMvsCC DEGs ", "FDR-", input$degfdr, " LFC-", input$degfc, ".csv", sep="") },
    content = function(file) {
      write.csv(datasetInput2(),file, row.names = F)
    }
  )
  
  output$table3 <- DT::renderDataTable({
    DEGFDR <- as.numeric(input$degfdr)/100
    DEGLFC <- as.numeric(input$degfc)
    if(input$genesymbolsonly == TRUE) { 
      tableS3 <- tableS3[tableS3$SYMBOL!= "",]}
    if(input$duplicated == TRUE) { 
      tableS3 <- tableS3[!duplicated(tableS3$SYMBOL),]}
    if(input$duplicated == TRUE & input$genesymbolsonly == TRUE) { 
      tableS3 <- tableS3[tableS3$SYMBOL!= "" & !duplicated(tableS3$SYMBOL),]} 
    tableS3[tableS3$adj.P.Val < DEGFDR & abs(tableS3$logFC) > DEGLFC,]
  })
  
  datasetInput3 <- reactive({
    DEGFDR <- as.numeric(input$degfdr)/100
    DEGLFC <- as.numeric(input$degfc)
    if(input$genesymbolsonly == TRUE) { 
      tableS3 <- tableS3[tableS3$SYMBOL!= "",]}
    if(input$duplicated == TRUE) { 
      tableS3 <- tableS3[!duplicated(tableS3$SYMBOL),]}
    if(input$duplicated == TRUE & input$genesymbolsonly == TRUE) { 
      tableS3 <- tableS3[tableS3$SYMBOL!= "" & !duplicated(tableS3$SYMBOL),]} 
    tableS3[tableS3$adj.P.Val < DEGFDR & abs(tableS3$logFC) > DEGLFC,]
  })
  
  output$exportTable3 <- downloadHandler(
    filename = function() { paste0("SMAvsCC DEGs ", "FDR-", input$degfdr, " LFC-", input$degfc, ".csv", sep="") },
    content = function(file) {
      write.csv(datasetInput3(),file, row.names = F)
    }
  )
  
  output$table4 <- DT::renderDataTable({
    DEGFDR <- as.numeric(input$degfdr)/100
    DEGLFC <- as.numeric(input$degfc)
    if(input$genesymbolsonly == TRUE) { 
      tableS4 <- tableS4[tableS4$SYMBOL!= "",]}
    if(input$duplicated == TRUE) { 
      tableS4 <- tableS4[!duplicated(tableS4$SYMBOL),]}
    if(input$duplicated == TRUE & input$genesymbolsonly == TRUE) { 
      tableS4 <- tableS4[tableS4$SYMBOL!= "" & !duplicated(tableS4$SYMBOL),]}  
        tableS4[tableS4$adj.P.Val < DEGFDR & abs(tableS4$logFC) > DEGLFC,]
  })
  
  datasetInput4 <- reactive({
    DEGFDR <- as.numeric(input$degfdr)/100
    DEGLFC <- as.numeric(input$degfc)
    if(input$genesymbolsonly == TRUE) { 
      tableS4 <- tableS4[tableS4$SYMBOL!= "",]}
    if(input$duplicated == TRUE) { 
      tableS4 <- tableS4[!duplicated(tableS4$SYMBOL),]}
    if(input$duplicated == TRUE & input$genesymbolsonly == TRUE) { 
      tableS4 <- tableS4[tableS4$SYMBOL!= "" & !duplicated(tableS4$SYMBOL),]}  
    tableS4[tableS4$adj.P.Val < DEGFDR & abs(tableS4$logFC) > DEGLFC,]
  })
  
  output$exportTable4 <- downloadHandler(
    filename = function() { paste0("CMvsSMA DEGs ", "FDR-", input$degfdr, " LFC-", input$degfc, ".csv", sep="") },
    content = function(file) {
      write.csv(datasetInput4(),file, row.names = F)
    }
  )
  
}


# Run the application 
shinyApp(ui = ui, server = server)
