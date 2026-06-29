# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Big section ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## ---- |-Subsetion: More details ----
library <- function(...) suppressPackageStartupMessages(base::library(...,quietly=TRUE))

# libraries
start_time <- Sys.time()

library(magrittr) # better pipes
library(tidyverse) 
library(scales)
library(readxl)
library(purrr)
library(tidyr)
library(readr)
#library(rethinking)
library(patchwork)
library(viridis)
library(tidylog) # loaded only to wrap dplyr verbs with row-count logging; detached again further down
library(summarytools)
library(here)
# library(dagitty)
library(ISOweek)
library(fst)
#library( "EcdcColors" )
#library(cmdstanr)
library(crayon)
library(lubridate)
# library(Hmisc)
# library(data.table)
# library(dtplyr)

# libraries not to include 
# library(tsibble)

# remasking
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
date <- lubridate::date
intersect <- base::intersect
setdiff <- base::setdiff
union <- base::union
#rstudent <- rethinking::rstudent
expand <- tidyr::expand
map <- purrr::map
discard <- purrr::discard
col_factor <- readr::col_factor
area <- patchwork::area
view <- summarytools::view
#compare <- rethinking::compare
# keeping only select functions from tidy_log 
. %>% dfSummary %>% view() -> viewsummary
filter_log <- tidylog::filter
left_join_log <- tidylog::left_join
fill_log <- tidylog::fill
detach(package:tidylog, unload = T)
# simplify calling functions
g = glimpse

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Super basic functions ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clc = function() cat("\014")
mstand <- function(v){
  mout <- v/sum(v)
  return(mout)
}
my_comma <- scales::label_comma(accuracy = 1, big.mark = ",", decimal.mark = ".")


get_in_hms <- function(time1, time2) {
  format(as.POSIXct(as.numeric(difftime(time1, time2, units = 'secs')), 
                    origin = '1970-01-01', tz = 'UTC'), '%H:%M:%S')  
}

# return the last element(s) of a vector
mlast <- function(v , n=1 ){
  mout = v[length(v)]
  if (n>1) {
    last_element = length(v)
    used_n = n
    if (n > length(v)) {
      used_n = min(length(v),n)
      cat(red("mlast uses a smaller n"))
    }
    first_element = last_element - used_n + 1
    mout = v[first_element:last_element]
  }
  return(mout)
}
# return the first element of a vector
mfirst <- function(v ){
  mout = v[1]
  return(mout)
}
# 
odds <- function(p){
  (p/(1-p)) -> mout
  return(mout)
}
inv_odds <- function( odds ){
  # body
  mp <- odds / (1+odds)
  return(mp)
} # try it: inv_odds( (0.4/0.6) )
odds_log <- function( p ) {
  log(p/(1-p)) -> mout
  return(mout)
}
logit <- odds_log
#
inv_logit = function (x) 
{
  p <- 1/(1 + exp(-x))
  p <- ifelse(x == Inf, 1, p)
  return(p)
}
#
replace_inf = function(vector_with_inf,replacement=NA){
  vector_with_inf[is.infinite(vector_with_inf)] = replacement
  return(vector_with_inf)
}

zero_plus_eps = function(vector_with_zeros,eps=1/(100*10^6) ){
  (vector_with_zeros==0) %>% sum() -> n_zeros
  #print(paste("zero_plus_eps() is adding a mass of:", n_zeros*eps))
  vector_with_zeros[vector_with_zeros==0] = eps
  return(vector_with_zeros)
}
if (F) {
  ##### adventures in re-scaling ...
  x_data = seq(from=0.5,to=99.5,length=50)
  pop_i = 100
  x_odds = odds(x_data/pop_i)
  plot(x_data,x_odds)
  x_logodds = log(x_odds)
  plot(x_data,x_logodds)
  x_sqrtodds = sqrt(x_odds)
  plot(x_data,x_sqrtodds)
  # backwards: from sqrtodds to data scale
  pop_i = 100
  x_sqrtodds = seq(from=-0.1,to=0.1,length=50)
  x_odds = x_sqrtodds^2
  x_data = (inv_odds(x_odds))*pop_i
  plot(x_sqrtodds,x_data)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Other Options ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# rstan_options(auto_write = TRUE) # stan options
n_chains <-  2
options(mc.cores = parallel::detectCores())
# Reset R's most annoying default options
options(stringsAsFactors = FALSE, 
        scipen = 999, 
        dplyr.summarise.inform = FALSE,
        tibble.print_min=4)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Settings for plotting ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ECDC Font
# if ("Tahoma" %in% extrafont::fonts()) {
#   FONT <- "Tahoma"
#   suppressMessages(extrafont::loadfonts(device = "win"))
# } else if (Sys.info()["sysname"] == "Windows") {
#   suppressMessages(extrafont::font_import(pattern = 'tahoma', prompt = FALSE))
#   suppressMessages(extrafont::loadfonts(device = "win"))
#   FONT <- "Tahoma"
# } else {
#   FONT <- NULL
# }

cleancat <- function(astring, width=80) {
  # Reserves a line of 80 (default) characters
  # and uses it for serial updates
  # require("stringr")
  astring <- paste("\r", astring, sep="")
  cat(str_pad(astring, 80, "right"))
}

# Theme (got this from WNV repo)
# .plottheme <- ggplot2::theme(axis.text = ggplot2::element_text(size = 8, family = FONT),
#                              axis.title = ggplot2::element_text(size = 9, family = FONT),
#                              axis.line = ggplot2::element_line(colour = "black"),
#                              axis.line.x = ggplot2::element_blank(),
#                              # --- Setting the background
#                              panel.grid.major = ggplot2::element_blank(),
#                              panel.grid.minor = ggplot2::element_blank(),
#                              panel.background = ggplot2::element_blank(),
#                              # --- Setting the legend
#                              legend.position = "right",
#                              legend.title = ggplot2::element_blank(),
#                              legend.text = ggplot2::element_text(size = 8, family = FONT),
#                              legend.key.width = ggplot2::unit(0.8, "cm"),
#                              legend.key.size = ggplot2::unit(0.4, "cm"))

# amazing: overriding function defaults
geom_interval <- function(...) ggdist::geom_interval(...,alpha=0.4)
geom_lineribbon <- function(...) ggdist::geom_lineribbon(...,alpha=0.4)
geom_ribbon <- function(...) ggplot2::geom_ribbon(...,alpha=0.4)
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Dark2")
mean_qi <- function(...) ggdist::mean_qi(...,.width=0.80)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Medium-complex functions ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
countries <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", 
               "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", 
               "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein", 
               "Lithuania", "Luxembourg", "Malta", "Netherlands", "Norway", 
               "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", 
               "Sweden",
               # non EU/EEA
               "Switzerland","England","Northern Ireland","Scotland","EU/EEA")
countries_short <- c("AT", "BE", "BG", "HR", "CY", "CZ", 
                     "DK", "EE", "FI", "FR", "DE", "GR", 
                     "HU", "IS", "IE", "IT", "LV", "LI", 
                     "LT", "LU", "MT", "NL", "NO", 
                     "PL", "PT", "RO", "SK", "SI", "ES", 
                     "SE",
                     # non EU/EEA
                     "CH","GB-ENG","GB-NIR","GB-SCT","EU/EEA")


# EL 
#EU_short("Greece") <- "EL"
# EU_short("Greece","EL")
EU_short <- function(name_long,greece="GR" # or "EL
){
  name_short = name_long
  for (i in 1:length(name_long)) {
    name_short[i] <- countries_short[which(countries%in%name_long[i])]
    if (name_long[i]=="Greece"&greece=="GR") name_short[i]<-"GR"
    if (name_long[i]=="Greece"&greece!="GR") name_short[i]<-"EL"
  }
  
  return(name_short)
}
#name_short=c("DE","PL","DE","GR"); EU_long(name_short,"EL")
EU_long <- function(name_short,greece="GR" # or "EL
){
  name_long = name_short
  for (i in 1:length(name_long)) {
    if (name_short[i]=="EL"&greece=="EL") name_short[i]<-"GR"
    name_long[i] <- countries[which(countries_short%in%name_short[i])]
    
  } 
  return(name_long)
}

# Custom ggplot-based density function (similar to the rethinking::dens() )
dens <- function(x, title = "", fill_color = "steelblue", line_color = "black", alpha = 0.5) {
  # Convert to data frame for ggplot
  df <- data.frame(Value = x)
  
  # Create density plot
  p <- ggplot(df, aes(x = Value)) +
    geom_density(aes(y = after_stat(density) ), fill = fill_color, color = line_color, alpha = alpha) +
    labs(title = title, x = "Value", y = "Density") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  print(p)
}

ecdc_weektodate <- function( year_week ){
  if ( any( nchar( year_week )!=7 ) ){
    stop( "Input must be of the format yyyy-ww !")
  }
  
  date_out <- ISOweek2date( paste0( substr( year_week, 1, 4 ), "-W", substr( year_week, 6, 7 ), "-1"  ) ) 
  
  return( date_out )
}

ecdc_datetoweek <- function( date_in ){
  if ( any( class( date_in )!="Date" ) ){
    stop( "Input must be a date !")
  }
  iso_week <- date2ISOweek( date_in )
  
  year_week <- paste0( substr( iso_week, 1, 4 ), "-", substr( iso_week, 7, 8 ) )
  return( year_week )
}

#  less simple functions, more specific to project
quantile_df <- function(x, probs = c(0.25, 0.5, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

column_stats_ingroups = function( df , mycolumn,mygroup , ... ) {
  mycolumn = enquo(mycolumn)
  mygroup = enquo(mygroup)
  
  mysumm = df %>% ungroup() %>% 
    reframe( quantile_df( !!mycolumn , ... ), 
             .by = !!mygroup )
  return(mysumm)
}

df_to_list = function(df){
  mlist = list()
  for (i in 1:nrow(df) ){
    mlist[[ (df$id[i]) ]] = df$data[i][[1]]
  }
  return(mlist)
}

ggsave_as <- function(p,figname,height=10,width=16){
  ggsave(plot=p,filename=paste0(here(), "/figures/",figname,".pdf"),
         height=height,width=width,unit="cm")
  
}

ggsave_as_png <- function(p,figname,height=10,width=16){
  ggsave(plot=p,filename=paste0(here(), "/figures/",figname,".png"),
         height=height,width=width,unit="cm")
  
}

# adventures is crossings

if (F) {
  # 1: build a df for a time-varying indicators for 2 countries
  df_1 = crossing( 
    location=c("DK","AT"), # will be sorted by FIRST vector
    t=c(1,2)
  )
  # 2: repeating the above for scenarios
  # these 3 options are identical
  df_1 %>% crossing( nesting(scenario=c("A","B"),disease_severity=c("high","low")) )
  df_1 %>% crossing( tibble(scenario=c("A","B"),disease_severity=c("high","low")) )
  df_1 %>% full_join( crossing(location=c("DK","AT"),nesting(scenario=c("A","B"),disease_severity=c("high","low"))),
                      relationship = "many-to-many")
  # 3: use the nest functionality to hide the time-varying part
  df_1 %>% crossing( nesting(scenario=c("A","B"),disease_severity=c("high","low")) ) %>% 
    nest(.by=c("scenario","disease_severity"))
}

# very end: timing
end_time <- Sys.time()
pr=paste("> Setup script run:",round(end_time - start_time,2),"sec \n"); cat(green(pr))
