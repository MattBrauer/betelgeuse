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

