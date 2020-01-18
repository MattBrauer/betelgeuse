library(tidyverse)
library(lubridate)
library(insol)
library(zoo)
library(forecast)
library(tsibble)

files_betelgeuse <- c("data/aavso_betelgeuse_2.txt", "data/aavso_betelgeuse_1.txt", "data/aavso_betelgeuse.txt")
files_AC_Her <- "data/aavso_AC_Her.txt"

col_specs <- cols(
  .default = col_character(),
  JD = col_double(),
  Magnitude = col_double(),
  Band = col_character(),
  `Observer Code` = col_character(),
  `Comment Code(s)` = col_character(),
  `Validation Flag` = col_character(),
  `Star Name` = col_character(),
  `Measurement Method` = col_character()
)

load_files <- function(files) {
  do.call(rbind, lapply(as.list(files),
                        function(file) {
                          read_csv(file, col_types = col_specs)
                        })) %>%
    distinct() %>%
    rename(jd = JD,
           observer = `Observer Code`,
           mag = Magnitude,
           band = Band) %>%
    arrange(jd) %>%
    #    mutate(next_obs = lead(observer),
    #           next_mag = lead(mag)) %>%
    select(jd, mag, observer, band)
}

daily_magnitude <- function(magnitudes) {
  magnitudes %>%
    mutate(day = floor(jd)) %>%
    select(day, mag, band, observer) %>%
    group_by(band, observer, day) %>%
    summarize(daily_observer_mag = mean(mag)) %>%
    group_by(band, day) %>%
    mutate(date = as_date(insol::JD(day, inverse = TRUE)),
           mean_daily_mag = mean(daily_observer_mag, na.rm=TRUE),
           sd_daily_mag = sd(daily_observer_mag, na.rm=TRUE),
           n_daily_mag = n(),
           delta_daily_mag = (daily_observer_mag - mean_daily_mag)) %>%
    group_by(band, observer) %>%
    arrange(desc(day))
}

raw <- load_files(files_betelgeuse)

daily_vis <- raw %>%
  filter(band=="Vis.") %>%
  daily_magnitude() %>%
  group_by(observer) %>%
  select(-band)

multiday <- daily_vis %>%
  filter(n_daily_mag > 1) %>%
  pull(day)

raw %>%
  mutate(day = floor(jd), date = as_date(insol::JD(day, inverse = TRUE))) %>%
  filter(band=="Vis.", day %in% multiday, date > '1965-01-01') %>%
  arrange(jd) %>%
  mutate(rolling_avg = rollapply(mag, 10, mean, fill=NA, partial=TRUE, align="right")) %>%
  ggplot(aes(x=date, y=rolling_avg)) + geom_point(size=0.1)

         


raw %>% mutate(outlier = if_else(observer %in% unreliable, TRUE, FALSE))

unreliable <- c("CTOA","KCD","KJMB","PGD","RPT","RTH","WJD")

prolific <- raw %>%
  filter(band=="Vis.") %>%
  group_by(observer) %>%
  summarize(n = n()) %>%
  filter(n > 50) %>%
  pull(observer)


daily %>%
  filter(day %in% multiday) %>%
#  filter(observer %in% prolific) %>%
#  filter(date > '2005-01-01') %>%
  ac_her %>%
  arrange(day) %>%
  mutate("rolling_avg" = rollapply(daily_mag, 10, mean, fill=NA),
         "rolling_count" = rollapply(daily_mag, 10, length, fill=NA)) %>%
  ggplot(mapping=aes(x=day, y=ma2, group=`Observer Code`)) +
  geom_line(aes(color = `Observer Code`)) +
  scale_color_discrete(guide = 'none')

  mutate(count_class = cut(n_daily_mag, breaks=c(0, 1, 2, 5, Inf), labels=c("1","2","5","many"))) %>%
  arrange(date) %>%
  ggplot(aes(x=date, y=daily_observer_mag, color=count_class)) + geom_point(aes(color=count_class))

daily %>%
  #filter(day > 2443600 & day < 2447120) %>%
  filter(day %in% multiday) %>%
  mutate(obs=if_else(observer %in% c("CRK", "RB", "ROP"), observer, "other")) %>%
  ggplot(aes(x=date, y=daily_observer_mag, color=obs)) + geom_point(aes(color=obs))


daily %>%
  filter(!observer %in% prolific) %>%
  arrange(date) %>%
  ggplot(aes(x=date, y=daily_observer_mag)) + geom_point()


betelgeuse %>%
  group_by(observer) %>%
  summarize(delta = mean(delta_daily_mag, na.rm=TRUE),
            nobs = n(),
            start = min(date),
            end = max(date)) %>%
  mutate("outlier" = if_else(abs(delta_daily_mag/sd_daily_mag) > 1, TRUE, FALSE)) %>%
#  arrange(desc(delta)) %>%
#  group_by(outlier) %>%
#  summarize(n = n()) %>%
  ggplot(aes(delta)) +
  geom_histogram(binwidth = 0.1)


betelgeuse %>% 
  select(date, observer, mean_daily_mag, delta_daily_mag) %>%
  pivot_wider(names_from = `Observer Code`, values_from = daily_observer_mag)






ac_her <- load_files(files_AC_Her)




ac_her %>%
  arrange(`Observer Code`, day) %>%
  mutate("rolling_avg" = rollapply(daily_mag, 10, mean, fill=NA),
         "rolling_count" = rollapply(daily_mag, 10, length, fill=NA)) %>%
  ggplot(mapping=aes(x=day, y=ma2, group=`Observer Code`)) +
  geom_line(aes(color = `Observer Code`)) +
  scale_color_discrete(guide = 'none')

ac_her %>%
  group_by(day) %>%
  summarize(daily_mag = mean(daily_observer_mag)) %>%
  arrange(day) %>%
  mutate("rolling_avg" = rollapply(daily_mag, 10, mean, fill=NA)) %>%
  ggplot(mapping=aes(x=day, y=rolling_avg)) +
  geom_line()

betelgeuse %>%
  group_by(date) %>%
  summarize(daily_mag = mean(daily_observer_mag)) %>%
  arrange(date) %>%
  mutate("rolling_avg" = rollapply(daily_mag, 10, mean, fill=NA)) %>%
  filter(date > '1965-01-01' & date < '2024-01-01') %>%
  ggplot(mapping=aes(x=date, y=rolling_avg)) +
  geom_point(alpha=0.5)

betelgeuse %>%
  group_by(day) %>%
  summarize(daily_mag = mean(daily_observer_mag)) %>%
  ggplot(mapping=aes(x=day, y=daily_mag)) +
  geom_line()


