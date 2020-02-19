#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(tidyverse)
library(stringr)
library(lubridate)
library(fs)
library(insol)
library(zoo)
library(forecast)
library(tsibble)
library(astrolibR)
library(httr)

get_new_data <- function(id=stringr::str_to_title(target), from_date='2010-01-01', to_date=current_date) {
    
    col_specs <- cols(
        .default = col_character(),
        JD = col_double(),
        mag = col_double(),
        uncert = col_double(),
        band = col_character(),
        by = col_character(),
        comCode = col_character(),
        val = col_character(),
        starName = col_character(),
        mtype = col_character(),
        obsID = col_skip(),
        fainterThan = col_skip(),
        obsType = col_skip(),
        software = col_skip(),
        obsName = col_skip(),
        obsCountry = col_skip()
    )
    
    fromjd <- JD(as.POSIXct(from_date))
    tojd <- JD(as.POSIXct(to_date))
    api_call <- paste0('https://www.aavso.org/vsx/index.php?view=api.delim&ident=%22',
                       id,
                       '%22&fromjd=',
                       fromjd,
                       '&tojd=',
                       tojd,
                       '&delimiter=@@@')
    
    response <- GET(api_call)
    
    new_data <- gsub("@@@", "\t", content(response))
    new_data <- gsub("\r", "", new_data)
    new_data %>%
        read_tsv(col_types = col_specs) %>%
        rename(jd = JD,
               observer = by,
               error = uncert) %>%
        filter(!is.na(mag)) %>%
        distinct(jd, observer, band, mag, error, .keep_all=TRUE)
}

add_dates <- function(magnitudes) {
    magnitudes %>%
        mutate(time = {if("jd" %in% names(.)) jd %% 1 else 0.5},
               day = {if("jd" %in% names(.)) round(jd) else {if("day" %in% names(.)) day else NA}},
               date = as_date(insol::JD(day, inverse = TRUE))) %>%
        separate(date, into=c("year","month","dom"), remove=FALSE) %>%
        mutate(year = as.integer(year),
               month = as.integer(month),
               dom = as.integer(dom),
               doy = (insol::daydoy(year, month, dom) + 175) %% 365)
}

### average multiple simultaneous observations by observer
clean_data <- function(raw_data) {
    raw_data %>%
        select(jd, mag, observer, band) %>%
        group_by(jd, observer, band) %>%
        summarize(delta = max(mag) - min(mag), mag = mean(mag), n = n()) %>%
        distinct() %>%
        ungroup() %>%
        select(jd, observer, band, mag, n, delta) %>%
        arrange(jd)
}

### calculate average daily magnitude where day spans [0.5, 0.5)
daily_magnitude <- function(magnitudes) {
    magnitudes %>%
        mutate(day = round(jd)) %>%
        select(day, mag, band, observer) %>%
        group_by(band, observer, day) %>%
        summarize(daily_observer_mag = mean(mag)) %>%
        group_by(band, day) %>%
        mutate(mean_daily_mag = mean(daily_observer_mag, na.rm=TRUE),
               sd_daily_mag = sd(daily_observer_mag, na.rm=TRUE),
               n_daily_mag = n(),
               delta_daily_mag = (daily_observer_mag - mean_daily_mag)) %>%
        ungroup() %>%
        arrange(day)
}

## get date of last observation in a set of observations
last_date <- function(ds) {
    ds %>%
        arrange(desc(jd)) %>%
        filter(row_number()==1) %>%
        add_dates() %>%
        pull(date)
} 

date_range <- function(ds) {
    ds %>%
        summarize(start_date=as_date(insol::JD(min(jd + 0.5), inverse=TRUE)),
                  end_date=as_date(insol::JD(max(jd + 0.5), inverse=TRUE)))
}

data_dir <- "data"

target <- "betelgeuse"
current_date <- today()
target_data_file <- fs::path(data_dir, paste0("aavso_", target), ext = "rds")

if(file_access(target_data_file, "read")) raw <- readRDS(target_data_file)
raw <- get_new_data(from_date=last_date(raw)) %>%
    bind_rows(raw) %>%
    distinct()
saveRDS(raw, file = target_data_file)
cleaned <- clean_data(raw)
latest_date <- last_date(cleaned)
dates_in_dataset <- date_range(cleaned) %>% c(.)

daily_mag <- cleaned %>%
    daily_magnitude() %>%
    distinct() %>%
    rename(mag = mean_daily_mag, n = n_daily_mag) %>%
    add_dates() %>%
    select(day, date, band, mag, n)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
#    titlePanel("AAVSO light curves"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            dateRangeInput(inputId = "daterange",
                           label = "Date range:",
                           start  = dates_in_dataset$end_date %m-% months(3),
                           end    = dates_in_dataset$end_date,
                           min    = dates_in_dataset$start_date,
                           max    = dates_in_dataset$end_date,
                           format = "mm/dd/yy",
                           separator = " - "
            ),
#            selectInput(inputId = "dateset",
#                        label = "Date quick set",
#                        choices = c("Past 3 months", "Past 6 months", "Past year"),
#                        selected = "Past 3 months",
#                        width = "220px"
#            ),
            numericRangeInput(inputId = "y_axis_range",
                              label = "y-axis limits:",
                              value = c(0, 4)
            ),
            selectInput(inputId = "bandselect",
                        label = "Band:",
                        choices = levels(as.factor(cleaned$band)),
                        selected = "Vis.",
                        width = "220px"
            ),
            radioButtons(inputId = "transformation",
                         label = "Data transformation:",
                         c("All observations" = "all",
                           "Daily average" = "average")
#                           "Daily average, >1 observation" = "avg_mult",
#                           "Running average" = "run")
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
            h3(textOutput("title")),
            plotOutput("magPlot"),
            textOutput("todayCaption"),
            textOutput("lastDateCaption"),
            textOutput("lastMeasurementCaption"),
            textOutput("lastAverageCaption")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    mag <- reactive({
        cleaned %>%
            filter(band==input$bandselect) %>%
            select(-band) %>%
            add_dates()
        })
    
    mag_ts <- reactive({
        mag() %>%
            as_tsibble(key = observer, index = jd, regular = FALSE)
    })
    
    avg_mag <- reactive({
        daily_mag %>%
            filter(band==input$bandselect) %>%
            select(-band) %>%
            distinct()
        })
    
    daily_mag_ts <- reactive({
        avg_mag() %>%
        as_tsibble(key=NULL, index=day) %>%
        fill_gaps() %>%
        mutate(in_run = if_else(is.na(date), FALSE, TRUE),
               last_of_run = lag(in_run, default = FALSE),
               next_of_run = lead(in_run, default = FALSE),
               run_break = in_run & !next_of_run)
    })
    
    output$title <- renderText({
        paste0("AAVSO Light curves for ", stringr::str_to_title(target))
    })
    
    output$todayCaption <- renderText({ 
        paste0("Today is ", current_date)
    })
    
    output$lastDateCaption <- renderText({ 
        paste0("Last date with measurement: ", dates_in_dataset$end_date)
    })
    
    output$lastMeasurementCaption <- renderText({ 
        paste0("Latest measurement = ",
               paste0(mag_ts() %>% top_n(1, day) %>% pull(mag), collapse = ','))
    })
    
    output$lastAverageCaption <- renderText({ 
        paste0("Latest average = ",
               paste0(avg_mag() %>% top_n(1, day) %>% pull(mag), collapse = ','))
    })
    
    output$magPlot <- renderPlot({
        
        plotdata <- switch(input$transformation,
                           all = mag_ts(),
                           average = daily_mag_ts(),
                           avg_mult = daily_mag_ts(),
                           run = daily_mag_ts(),
                           mag_ts())
        
        plottitle <- switch(input$transformation,
                            all = "all observations",
                            average = "daily average",
                            avg_mult = "daily average (>1 obs. per day)",
                            run = "running average",
                            mag_ts)
        
        plotdata %>%
            filter(date > input$daterange[1],
                   mag > input$y_axis_range[1] & mag < input$y_axis_range[2]) %>%
            ggplot(aes(x=date, y=mag)) + geom_point() + 
            scale_y_reverse() +
            ggtitle(paste0(input$bandselect, " band, ", plottitle))
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
