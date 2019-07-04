# Description ------------------------------------

# Opens viewer for plot files of peak windows from a csv table given a directory
# of MS1 tables

# User Input ------------------------------------------------------------------

# Path to folder with individual peak table directories
# for Windows operating system use \\ or / instead of \ to separate directories
# example:
# peakTableParentDir <- 'M:\\folder1\\foldwerWithFile'
imageDir <- "M:\\folder1\\foldwerWithFile"

# path to folder with tracking table
trackFolderPath <- "M:\\folder2\\foldwerWithFile"

# tracking table file name
trackingTableCSV <- "filename.csv"

# End user input ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# Library (packages)--------------------------------------------------------
library(shiny) # used for generating html apps

library(data.table)

# load tracking table ----------------------------------------------------

pathTrackingTable <- paste0(trackFolderPath,"/",trackingTableCSV)

dt.tracking <- fread(file = pathTrackingTable)

# get number of image rows in table
numTotalImages <- nrow(dt.tracking)

# initialize image number --------------------------------------------------

# sort by rowNum
setorder(dt.tracking, rowNum)

# subset to unviewed
dt.unviewed <- dt.tracking[viewed == 0]

# take first rowNum or set to last if all have been viewed
if (nrow(dt.unviewed) > 0) {
  startRow <- dt.unviewed$rowNum[1]
} else {
  startRow <- numTotalImages
}

# remove "unviewed"
rm(dt.unviewed)

numImageFile <- startRow

numImageFilePrevSave <- startRow

dt.tracking$viewed[startRow] <- 1

# define the User interface -------------------------------------------------

ui <- fluidPage(
  verticalLayout(
    # Use imageOutput to place the image on the page
    imageOutput(outputId = "myImage"),
    # breaklines so buttons aren't on top of plot
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    fluidRow(textOutput("textBad"),
             textOutput("textWindFlag"),
             textOutput("textBorderline")
             ),
    fluidRow(
      # HTML for hotkeys
      # find hotkeys at: http://keycode.info/
      tags$script(HTML("$(function(){
      $(document).keyup(function(e) {
      if (e.which == 32) {
        $('#buttonBad').click()
      }
      if (e.which == 71) {
        $('#buttonGood').click()
      }
      if (e.which == 37) {
        $('#buttonPrev').click()
      }
      if (e.which == 39) {
        $('#buttonNext').click()
      }
      if (e.which == 83) {
        $('#buttonSave').click()
      }
      if (e.which == 87) {
        $('#buttonWindFlag').click()
      }
      if (e.which == 66) {
        $('#buttonBorderline').click()
      }
      });
      })")),
      column(width = 1,
             actionButton(inputId = "buttonSave", label = "Save (s)"),
             textOutput("textSave")
             ),
      column(width = 2, offset = 1,
             actionButton(inputId = "buttonBad", label = "Bad (Space)"),
             actionButton(inputId = "buttonGood", label = "Good (g)")
             ),
      column(width = 3, offset = 1,
             actionButton(inputId = "buttonPrev", label = "Previous (<-)"),
             actionButton(inputId = "buttonNext", label = "Next (->)"),
             textOutput("textFileNum")
             ),
      column(width = 3, offset = 1,
             actionButton(inputId = "buttonWindFlag",
                          label = "(un)Flag Window Edges (w)"),
             actionButton(inputId = "buttonBorderline",
                          label = "(un)Flag Borderline (b)")
             ),
      # formatting for Bad text and Window flag
      tags$head(tags$style("#textBad{color: red;
                                 font-size: 40px;
                           font-style: bold;
                           }"
      )),
      tags$head(tags$style("#textWindFlag{color: red;
                                 font-size: 20px;
                                     font-style: bold;
                                     }"
      )),
      tags$head(tags$style("#textBorderline{color: red;
                                 font-size: 20px;
                                     font-style: bold;
                                     }"
      ))# end formatting
    ) # end well Panel
  ) # end vertical layout
) # end fluid page




# Define the server code ----------------------------------------------------


server <- function(input, output, session) {

  # use <<- to change global objects

  # create reactive state variables -----------------------------------------

  state <- reactiveValues()

  state$numImageFile <- startRow

  observe({

    # track the total button clicks
    state$sumButtonActivityCount <- input$buttonNext +
      input$buttonPrev + input$buttonBad

    # update image label
    state$peakGood <- dt.tracking$good1bad0[state$numImageFile]

    # update window bounds
    state$windowBad <- dt.tracking$windowBad[state$numImageFile]

    # update borderline variable
    state$borderline <- dt.tracking$borderline[state$numImageFile]

  })

  # track image information -------------------------------------------------------

  # update values when previous and next buttons are pressed
  observeEvent(input$buttonNext, {

    state$numImageFile <- min(numTotalImages, state$numImageFile + 1)
    # update viewed column
    dt.tracking$viewed[state$numImageFile] <<- 1
  })
  observeEvent(input$buttonPrev, {

    # advance image number
    state$numImageFile <- max(0, state$numImageFile - 1)

    # update viewed column
    dt.tracking$viewed[state$numImageFile] <<- 1
  })

  # text for tracking image number
  output$textFileNum <- renderText({
    # text output
    paste0("image number: ", state$numImageFile, " of ", numTotalImages)
  })


  # Display image -----------------------------------------------------------

  # A dynamically-sized plot
  output$myImage <- renderImage({

    # get path for current image file
    imagePath <- dt.tracking$path[state$numImageFile]

    # Read myImage's width and height. These are reactive values, so this
    # expression will re-run whenever they change.
    width  <- session$clientData$output_myImage_width
    # height <- session$clientData$output_myImage_height

    # Return a list containing the filename
    list(src = imagePath,
         contentType = 'image/jpg',
         width = width,
         alt = "No More Images")
  }, deleteFile = FALSE) # end render image

  # Indicate Image Bad -------------------------------------------------------

  observeEvent(input$buttonBad, {

    dt.tracking$good1bad0[state$numImageFile] <<- 0

    state$peakGood <- dt.tracking$good1bad0[state$numImageFile]

  })
  observeEvent(input$buttonGood, {

    dt.tracking$good1bad0[state$numImageFile] <<- 1

    state$peakGood <- dt.tracking$good1bad0[state$numImageFile]

  })

  output$textBad <- renderText({

    if (state$peakGood == 0) {
      "BAD"
    } else {""}

  })

  # flag window edges ---------------------------------------------------------------

  observeEvent(input$buttonWindFlag, {

    if (state$windowBad == 0) {

      dt.tracking$windowBad[state$numImageFile] <<- 1

    } else {

      dt.tracking$windowBad[state$numImageFile] <<- 0

    }

    state$windowBad <- dt.tracking$windowBad[state$numImageFile]

  })

  output$textWindFlag <- renderText({

    if (state$windowBad == 1) {
      "Window Edges Flagged"
    } else {""}

  })


  # flag as borderline ---------------------------------------------------------------



  observeEvent(input$buttonBorderline, {

    if (state$borderline == 0) {

      dt.tracking$borderline[state$numImageFile] <<- 1

    } else {

      dt.tracking$borderline[state$numImageFile] <<- 0

    }

    state$borderline <- dt.tracking$borderline[state$numImageFile]

  })

  output$textBorderline <- renderText({

    if (state$borderline == 1) {
      "Peak Flagged as Borderline Acceptable"
    } else {""}

  })

  # Saving table --------------------------------------------------------------------

  # initialize save tracker
  state$saveActivityVal <- 0

  observeEvent(input$buttonSave, {

    # write the tracking table
    fwrite(dt.tracking,pathTrackingTable)

    # update state$saveValPrev
    state$saveActivityVal <- state$sumButtonActivityCount

    # print to console
    print(paste0("save",Sys.time(), " ", state$saveActivityVal))

  })

  # text for saving
  output$textSave <- renderText({

    if (state$saveActivityVal == state$sumButtonActivityCount) {
      "Progress Saved"
    } else {""} # blank if not saved

  }) # end save text element

} # end server function

# Return a Shiny app object ---------------------------------------------------
shinyApp(ui = ui, server = server)

