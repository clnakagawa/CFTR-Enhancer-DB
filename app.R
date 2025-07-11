library(shinydashboard)
library(bslib)
library(DBI)
library(RSQLite)
library(DT)

tab <- function(...) {
  shiny::tabPanel(..., class = "p-3 border border-top-0 rounded-bottom")
}

# Evidence options keyed by table (for easy manual editing per association)
evidenceChoices <- list(
  Activity = c("ATAC", "DNAse", "Histone"),
  Contact = c("3C", "4C", "5C", "Hi-C"),
  Target = c("Reporter Assay", "ChIP", "ChIP-seq", "EMSA", "Supershift")
)

# Tabs info for viewing tables
tabsInfo <- list(
  list(title = "Regions", tableId = "Region", sql = "SELECT * FROM Region"),
  list(title = "Region Contacts", tableId = "Contact", sql = "SELECT r1.name AS region1, r2.name AS region2,
       c.name AS cell_type, c.tissue AS tissue, s.doi AS doi, evidence AS experiment, activity
       FROM Contact JOIN Region AS r1 ON Contact.rid1 = r1.rid
       JOIN Region AS r2 ON Contact.rid2 = r2.rid
       JOIN Cell AS c ON Contact.cid = c.cid
       JOIN Source AS s ON Contact.sid = s.sid;"),
  list(title = "Region Targeting", tableId = "Target", sql = "SELECT r.name as region, c.name as cell_type,
       c.tissue AS tissue, f.name as factor, s.doi as doi, evidence AS experiment, activity
       FROM Target JOIN Region AS r ON Target.rid = r.rid
       JOIN Cell AS c ON Target.cid = c.cid
       JOIN Factor AS f ON Target.fid = f.fid
       JOIN Source AS s ON Target.sid = s.sid;"),
  list(title = "Region Activity", tableId = "Activity", sql = "SELECT r.name AS region, c.name AS cell_type,
       c.tissue AS tissue, s.doi AS doi, evidence AS experiment, activity
       FROM Activity JOIN Region AS r ON Activity.rid = r.rid
       JOIN Cell AS c ON Activity.cid = c.cid
       JOIN Source AS s ON Activity.sid = s.sid;")
)

ui <- bootstrapPage(
  tags$div(
    style = "background-color:#028BBF; padding: 15px; width: 100%;",
    tags$h1("CFTR Regulatory Region Browser", style = "margin-bottom: 5px; color: white;"),
    tags$p("Interface for searching known information about CFTR enhancer, chromosome interactions, and TF targets",
           style = "margin-top: 0; font-size: 16px; color: white;")
  ),
  
  tags$div(
    style = "margin: 20px 40px;",
    tags$style(HTML("
      .nav.nav-tabs {
        margin-bottom: 20px;
      }
      .dataTables_wrapper .dataTables_info,
      .dataTables_wrapper .dataTables_paginate {
        margin-top: 40px;  
      }
      .dataTables_scrollBody {
        overflow-y: hidden!important;
       }
    ")),
    
    do.call(tabsetPanel, c(
      lapply(tabsInfo, function(tabDef) {
        tabPanel(
          tabDef$title,
          div(
            style = "padding-bottom: 80px;",
            DTOutput(outputId = paste0("table", tabDef$tableId))
          )
        )
      }),
      list(
        tabPanel("Upload",
                 h4("Add New Record"),
                 selectInput("formSelector", "Select table to add data to:",
                             choices = c("Region", "Cell", "Factor", "Activity", "Contact", "Target")),
                 uiOutput("formUi"),
                 verbatimTextOutput("submitStatus")
        )
      )
    ))
  )
)



server <- function(input, output, session) {
  conn <- dbConnect(RSQLite::SQLite(), "enhancerDB.sqlite")
  
  # refresh triggers
  refreshTriggers <- list(
    Region = reactiveVal(0),
    Contact = reactiveVal(0),
    Target = reactiveVal(0),
    Activity = reactiveVal(0)
  )
  
  # Load data tables for viewing
  for (tabDef in tabsInfo) {
    local({
      tabCopy <- tabDef
      output[[paste0("table", tabCopy$tableId)]] <- renderDT({
        refreshTriggers[[tabCopy$tableId]]()
        datatable(
          dbGetQuery(conn, tabCopy$sql),
          filter = "top",
          options = list(pageLength = 10)
        )
      })
    })
  }
  
  # Helper to insert generic records with auto-increment IDs
  insertRecord <- function(table, idCol, fields, values, outputId) {
    idQuery <- sprintf("SELECT IFNULL(MAX(%s), 0) + 1 AS nextId FROM %s", idCol, table)
    newId <- dbGetQuery(conn, idQuery)$nextId
    
    fieldList <- paste(c(idCol, fields), collapse = ", ")
    placeholderList <- paste(rep("?", length(fields) + 1), collapse = ", ")
    sql <- sprintf("INSERT INTO %s (%s) VALUES (%s)", table, fieldList, placeholderList)
    
    tryCatch({
      dbExecute(conn, sql, params = c(newId, values))
      output[[outputId]] <- renderText(paste(table, "added with ID:", newId))
      TRUE
    }, error = function(e) {
      output[[outputId]] <- renderText(paste("Failed to insert into", table, ":", e$message))
      FALSE
    })
  }
  
  # Dynamic UI for form
  output$formUi <- renderUI({
    switch(input$formSelector,
           "Region" = tagList(
             textInput("regionName", "Region Name"),
             selectInput("regionType", "Region Type", c("enhancer", "silencer", "promoter", "other")),
             selectInput("regionChr", "Chromosome", paste0("chr", c(1:22, "X", "Y"))),
             numericInput("regionStart", "Start", value = 1, min = 1),
             numericInput("regionEnd", "End", value = 1000, min = 1),
             actionButton("submitRegion", "Add Region")
           ),
           "Cell" = tagList(
             textInput("cellName", "Cell Name"),
             textInput("cellTissue", "Tissue"),
             actionButton("submitCell", "Add Cell")
           ),
           "Factor" = tagList(
             textInput("factorName", "Factor Name"),
             actionButton("submitFactor", "Add Factor")
           ),
           "Activity" = {
             regions <- dbGetQuery(conn, "SELECT rid, name FROM Region ORDER BY name")
             cells <- dbGetQuery(conn, "SELECT cid, name FROM Cell ORDER BY name")
             sources <- dbGetQuery(conn, "SELECT sid, doi FROM Source ORDER BY doi")
             
             fluidRow(
               column(6,
                      div(style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 10px; border-radius: 5px;",
                          tags$b("Region Data"),
                          selectInput("activityRegion", NULL,
                                      choices = c(setNames(as.character(regions$rid), regions$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.activityRegion == 'other'",
                            textInput("activityRegionName", "New Region Name"),
                            selectInput("activityRegionType", "Region Type", c("enhancer", "silencer", "promoter", "other")),
                            selectInput("activityRegionChr", "Chromosome", paste0("chr", c(1:22, "X", "Y"))),
                            numericInput("activityRegionStart", "Start", value = 1, min = 1),
                            numericInput("activityRegionEnd", "End", value = 1000, min = 1)
                          )
                      )
               ),
               column(6,
                      div(style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 10px; border-radius: 5px;",
                          tags$b("Cell Data"),
                          selectInput("activityCell", NULL,
                                      choices = c(setNames(as.character(cells$cid), cells$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.activityCell == 'other'",
                            textInput("activityCellName", "New Cell Name"),
                            textInput("activityCellTissue", "Tissue")
                          )
                      )
               ),
               column(12,
                      div(style = "border: 1px solid #ccc; padding: 15px; border-radius: 5px;",
                          tags$b("Activity Details"),
                          textInput("activitySourceDoi", "Source DOI (required)"),
                          selectInput("activityEvidence", "Evidence",
                                      choices = evidenceChoices$Activity),
                          selectInput("activityActivity", "Activity", choices = c("Yes", "No", "Unsure")),
                          actionButton("submitActivity", "Add Activity")
                      )
               )
             )
           },
           # Contact form UI snippet (replace old one)
           "Contact" = {
             regions <- dbGetQuery(conn, "SELECT rid, name FROM Region ORDER BY name")
             cells <- dbGetQuery(conn, "SELECT cid, name FROM Cell ORDER BY name")
             sources <- dbGetQuery(conn, "SELECT sid, doi FROM Source ORDER BY doi")
             
             fluidRow(
               column(6,
                      div(style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 10px; border-radius: 5px;",
                          tags$b("Region 1 Data"),
                          selectInput("contactRegion1", NULL,
                                      choices = c(setNames(as.character(regions$rid), regions$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.contactRegion1 == 'other'",
                            textInput("contactRegion1Name", "New Region 1 Name"),
                            selectInput("contactRegion1Type", "Region 1 Type", c("enhancer", "silencer", "promoter", "other")),
                            selectInput("contactRegion1Chr", "Chromosome", paste0("chr", c(1:22, "X", "Y"))),
                            numericInput("contactRegion1Start", "Start", value = 1, min = 1),
                            numericInput("contactRegion1End", "End", value = 1000, min = 1)
                          )
                      )
               ),
               column(6,
                      div(style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 10px; border-radius: 5px;",
                          tags$b("Region 2 Data"),
                          selectInput("contactRegion2", NULL,
                                      choices = c(setNames(as.character(regions$rid), regions$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.contactRegion2 == 'other'",
                            textInput("contactRegion2Name", "New Region 2 Name"),
                            selectInput("contactRegion2Type", "Region 2 Type", c("enhancer", "silencer", "promoter", "other")),
                            selectInput("contactRegion2Chr", "Chromosome", paste0("chr", c(1:22, "X", "Y"))),
                            numericInput("contactRegion2Start", "Start", value = 1, min = 1),
                            numericInput("contactRegion2End", "End", value = 1000, min = 1)
                          )
                      )
               ),
               column(12,
                      div(style = "border: 1px solid #ccc; padding: 15px; border-radius: 5px;",
                          tags$b("Cell Data"),
                          selectInput("contactCell", NULL,
                                      choices = c(setNames(as.character(cells$cid), cells$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.contactCell == 'other'",
                            textInput("contactCellName", "New Cell Name"),
                            textInput("contactCellTissue", "Tissue")
                          )
                      )
               ),
               column(12,
                      div(style = "border: 1px solid #ccc; padding: 15px; border-radius: 5px;",
                          tags$b("Contact Details"),
                          textInput("contactSourceDoi", "Source DOI (required)"),
                          selectInput("contactEvidence", "Evidence", choices = evidenceChoices$Contact),
                          selectInput("contactActivity", "Activity", choice = c("Yes", "No", "Unsure")),
                          actionButton("submitContact", "Add Contact")
                      )
               )
             )
           },
           # Target form UI snippet (replace old one)
           "Target" = {
             regions <- dbGetQuery(conn, "SELECT rid, name FROM Region ORDER BY name")
             factors <- dbGetQuery(conn, "SELECT fid, name FROM Factor ORDER BY name")
             cells <- dbGetQuery(conn, "SELECT cid, name FROM Cell ORDER BY name")
             sources <- dbGetQuery(conn, "SELECT sid, doi FROM Source ORDER BY doi")
             
             fluidRow(
               column(4,
                      div(style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 10px; border-radius: 5px;",
                          tags$b("Region Data"),
                          selectInput("targetRegion", NULL,
                                      choices = c(setNames(as.character(regions$rid), regions$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.targetRegion == 'other'",
                            textInput("targetRegionName", "New Region Name"),
                            selectInput("targetRegionType", "Region Type", c("enhancer", "silencer", "promoter", "other")),
                            selectInput("targetRegionChr", "Chromosome", paste0("chr", c(1:22, "X", "Y"))),
                            numericInput("targetRegionStart", "Start", value = 1, min = 1),
                            numericInput("targetRegionEnd", "End", value = 1000, min = 1)
                          )
                      )
               ),
               column(4,
                      div(style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 10px; border-radius: 5px;",
                          tags$b("Factor Data"),
                          selectInput("targetFactor", NULL,
                                      choices = c(setNames(as.character(factors$fid), factors$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.targetFactor == 'other'",
                            textInput("targetFactorName", "New Factor Name")
                          )
                      )
               ),
               column(4,
                      div(style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 10px; border-radius: 5px;",
                          tags$b("Cell Data"),
                          selectInput("targetCell", NULL,
                                      choices = c(setNames(as.character(cells$cid), cells$name), "Other" = "other")),
                          conditionalPanel(
                            condition = "input.targetCell == 'other'",
                            textInput("targetCellName", "New Cell Name"),
                            textInput("targetCellTissue", "Tissue")
                          )
                      )
               ),
               column(12,
                      div(style = "border: 1px solid #ccc; padding: 15px; border-radius: 5px;",
                          tags$b("Target Details"),
                          textInput("targetSourceDoi", "Source DOI (required)"),
                          selectInput("targetEvidence", "Evidence", choices = evidenceChoices$Target),
                          selectInput("targetActivity", "Activity", choices = c("Yes", "No", "Unsure")),
                          actionButton("submitTarget", "Add Target")
                      )
               )
             )
           }
    )
  })
  
  # Region insert handler
  observeEvent(input$submitRegion, {
    req(input$regionName, input$regionType, input$regionChr,
        input$regionStart, input$regionEnd)
    
    if (input$regionStart >= input$regionEnd) {
      output$submitStatus <- renderText("Error: start must be less than end.")
      return()
    }
    
    insertRecord(
      table = "Region",
      idCol = "rid",
      fields = c("chr", "start", "end", "type", "name"),
      values = list(
        input$regionChr,
        input$regionStart,
        input$regionEnd,
        input$regionType,
        input$regionName
      ),
      outputId = "submitStatus"
    )
    refreshTriggers$Region(refreshTriggers$Region() + 1)
  })
  
  # Cell insert handler
  observeEvent(input$submitCell, {
    req(input$cellName, input$cellTissue)
    insertRecord(
      table = "Cell",
      idCol = "cid",
      fields = c("name", "tissue"),
      values = list(input$cellName, input$cellTissue),
      outputId = "submitStatus"
    )
  })
  
  # Factor insert handler
  observeEvent(input$submitFactor, {
    req(input$factorName)
    insertRecord(
      table = "Factor",
      idCol = "fid",
      fields = c("name"),
      values = list(input$factorName),
      outputId = "submitStatus"
    )
  })
  
  # Activity insert handler
  observeEvent(input$submitActivity, {
    if (is.null(input$activitySourceDoi) || input$activitySourceDoi == "") {
      output$submitStatus <- renderText("Error: Source DOI is required.")
      return()
    }
    
    getOrCreateSource <- function(doi) {
      existing <- dbGetQuery(conn, "SELECT sid FROM Source WHERE doi = ?", params = list(doi))
      if (nrow(existing) == 1) {
        return(existing$sid)
      } else {
        newSid <- dbGetQuery(conn, "SELECT IFNULL(MAX(sid), 0) + 1 AS nextId FROM Source")$nextId
        dbExecute(conn, "INSERT INTO Source (sid, doi) VALUES (?, ?)", params = list(newSid, doi))
        return(newSid)
      }
    }
    
    # Handle Region
    if (input$activityRegion == "other") {
      req(input$activityRegionName, input$activityRegionType, input$activityRegionChr,
          input$activityRegionStart, input$activityRegionEnd)
      if (input$activityRegionStart >= input$activityRegionEnd) {
        output$submitStatus <- renderText("Error: Region start must be less than end.")
        return()
      }
      newRid <- dbGetQuery(conn, "SELECT IFNULL(MAX(rid), 0) + 1 AS nextId FROM Region")$nextId
      dbExecute(conn, "INSERT INTO Region (rid, chr, start, end, type, name) VALUES (?, ?, ?, ?, ?, ?)",
                params = list(newRid, input$activityRegionChr, input$activityRegionStart,
                              input$activityRegionEnd, input$activityRegionType, input$activityRegionName))
      rid <- newRid
    } else {
      rid <- as.integer(input$activityRegion)
    }
    
    # Handle Cell
    if (input$activityCell == "other") {
      req(input$activityCellName, input$activityCellTissue)
      newCid <- dbGetQuery(conn, "SELECT IFNULL(MAX(cid), 0) + 1 AS nextId FROM Cell")$nextId
      dbExecute(conn, "INSERT INTO Cell (cid, name, tissue) VALUES (?, ?, ?)",
                params = list(newCid, input$activityCellName, input$activityCellTissue))
      cid <- newCid
    } else {
      cid <- as.integer(input$activityCell)
    }
    
    sid <- getOrCreateSource(input$activitySourceDoi)
    
    tryCatch({
      dbExecute(conn, "INSERT INTO Activity (rid, cid, sid, evidence, activity) VALUES (?, ?, ?, ?, ?)",
                params = list(rid, cid, sid, input$activityEvidence, input$activityActivity))
      output$submitStatus <- renderText("Activity record added successfully.")
      refreshTriggers$Activity(refreshTriggers$Activity() + 1)
    }, error = function(e) {
      output$submitStatus <- renderText(paste("Failed to insert Activity:", e$message))
    })
  })
  
  # Contact insert handler
  observeEvent(input$submitContact, {
    if (is.null(input$contactSourceDoi) || input$contactSourceDoi == "") {
      output$submitStatus <- renderText("Error: Source DOI is required.")
      return()
    }
    
    getOrCreateSource <- function(doi) {
      existing <- dbGetQuery(conn, "SELECT sid FROM Source WHERE doi = ?", params = list(doi))
      if (nrow(existing) == 1) {
        return(existing$sid)
      } else {
        newSid <- dbGetQuery(conn, "SELECT IFNULL(MAX(sid), 0) + 1 AS nextId FROM Source")$nextId
        dbExecute(conn, "INSERT INTO Source (sid, doi) VALUES (?, ?)", params = list(newSid, doi))
        return(newSid)
      }
    }
    
    # Helper to handle region inputs (1 or 2)
    handleRegion <- function(prefix) {
      regionInput <- input[[paste0("contactRegion", prefix)]]
      if (regionInput == "other") {
        req(input[[paste0("contactRegion", prefix, "Name")]],
            input[[paste0("contactRegion", prefix, "Type")]],
            input[[paste0("contactRegion", prefix, "Chr")]],
            input[[paste0("contactRegion", prefix, "Start")]],
            input[[paste0("contactRegion", prefix, "End")]])
        if (input[[paste0("contactRegion", prefix, "Start")]] >= input[[paste0("contactRegion", prefix, "End")]]) {
          stop(sprintf("Error: Region %s start must be less than end.", prefix))
        }
        newRid <- dbGetQuery(conn, "SELECT IFNULL(MAX(rid), 0) + 1 AS nextId FROM Region")$nextId
        dbExecute(conn, "INSERT INTO Region (rid, chr, start, end, type, name) VALUES (?, ?, ?, ?, ?, ?)",
                  params = list(
                    newRid,
                    input[[paste0("contactRegion", prefix, "Chr")]],
                    input[[paste0("contactRegion", prefix, "Start")]],
                    input[[paste0("contactRegion", prefix, "End")]],
                    input[[paste0("contactRegion", prefix, "Type")]],
                    input[[paste0("contactRegion", prefix, "Name")]]
                  ))
        return(newRid)
      } else {
        return(as.integer(regionInput))
      }
    }
    
    # Handle both regions
    rid1 <- NULL
    rid2 <- NULL
    tryCatch({
      rid1 <- handleRegion("1")
      rid2 <- handleRegion("2")
    }, error = function(e) {
      output$submitStatus <- renderText(e$message)
      return()
    })
    
    # Handle Cell
    if (input$contactCell == "other") {
      req(input$contactCellName, input$contactCellTissue)
      newCid <- dbGetQuery(conn, "SELECT IFNULL(MAX(cid), 0) + 1 AS nextId FROM Cell")$nextId
      dbExecute(conn, "INSERT INTO Cell (cid, name, tissue) VALUES (?, ?, ?)",
                params = list(newCid, input$contactCellName, input$contactCellTissue))
      cid <- newCid
    } else {
      cid <- as.integer(input$contactCell)
    }
    
    sid <- getOrCreateSource(input$contactSourceDoi)
    
    tryCatch({
      dbExecute(conn, "INSERT INTO Contact (rid1, rid2, cid, sid, evidence, activity) VALUES (?, ?, ?, ?, ?, ?)",
                params = list(rid1, rid2, cid, sid, input$contactEvidence, input$contactActivity))
      output$submitStatus <- renderText("Contact record added successfully.")
      refreshTriggers$Contact(refreshTriggers$Contact() + 1)
    }, error = function(e) {
      output$submitStatus <- renderText(paste("Failed to insert Contact:", e$message))
    })
  })
  
  # Target insert handler
  observeEvent(input$submitTarget, {
    if (is.null(input$targetSourceDoi) || input$targetSourceDoi == "") {
      output$submitStatus <- renderText("Error: Source DOI is required.")
      return()
    }
    
    getOrCreateSource <- function(doi) {
      existing <- dbGetQuery(conn, "SELECT sid FROM Source WHERE doi = ?", params = list(doi))
      if (nrow(existing) == 1) {
        return(existing$sid)
      } else {
        newSid <- dbGetQuery(conn, "SELECT IFNULL(MAX(sid), 0) + 1 AS nextId FROM Source")$nextId
        dbExecute(conn, "INSERT INTO Source (sid, doi) VALUES (?, ?)", params = list(newSid, doi))
        return(newSid)
      }
    }
    
    # Handle Region
    if (input$targetRegion == "other") {
      req(input$targetRegionName, input$targetRegionType, input$targetRegionChr,
          input$targetRegionStart, input$targetRegionEnd)
      if (input$targetRegionStart >= input$targetRegionEnd) {
        output$submitStatus <- renderText("Error: Region start must be less than end.")
        return()
      }
      newRid <- dbGetQuery(conn, "SELECT IFNULL(MAX(rid), 0) + 1 AS nextId FROM Region")$nextId
      dbExecute(conn, "INSERT INTO Region (rid, chr, start, end, type, name) VALUES (?, ?, ?, ?, ?, ?)",
                params = list(newRid, input$targetRegionChr, input$targetRegionStart,
                              input$targetRegionEnd, input$targetRegionType, input$targetRegionName))
      rid <- newRid
    } else {
      rid <- as.integer(input$targetRegion)
    }
    
    # Handle Factor
    if (input$targetFactor == "other") {
      req(input$targetFactorName)
      newFid <- dbGetQuery(conn, "SELECT IFNULL(MAX(fid), 0) + 1 AS nextId FROM Factor")$nextId
      dbExecute(conn, "INSERT INTO Factor (fid, name) VALUES (?, ?)",
                params = list(newFid, input$targetFactorName))
      fid <- newFid
    } else {
      fid <- as.integer(input$targetFactor)
    }
    
    # Handle Cell
    if (input$targetCell == "other") {
      req(input$targetCellName, input$targetCellTissue)
      newCid <- dbGetQuery(conn, "SELECT IFNULL(MAX(cid), 0) + 1 AS nextId FROM Cell")$nextId
      dbExecute(conn, "INSERT INTO Cell (cid, name, tissue) VALUES (?, ?, ?)",
                params = list(newCid, input$targetCellName, input$targetCellTissue))
      cid <- newCid
    } else {
      cid <- as.integer(input$targetCell)
    }
    
    sid <- getOrCreateSource(input$targetSourceDoi)
    
    tryCatch({
      dbExecute(conn, "INSERT INTO Target (rid, fid, cid, sid, evidence, activity) VALUES (?, ?, ?, ?, ?, ?)",
                params = list(rid, fid, cid, sid, input$targetEvidence, input$targetActivity))
      output$submitStatus <- renderText(
        paste0("Target record added successfully"))
      refreshTriggers$Target(refreshTriggers$Target() + 1)
    }, error = function(e) {
      output$submitStatus <- renderText(paste("Failed to insert Target:", e$message))
    })
  })
  
  
  # Disconnect DB on session end
  session$onSessionEnded(function() {
    dbDisconnect(conn)
  })
}

shinyApp(ui, server)
