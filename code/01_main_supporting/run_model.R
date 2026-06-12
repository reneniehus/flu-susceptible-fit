run_model = function( params=NULL , data=NULL , models_in=NULL ){
  t1 <- Sys.time()
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Initiating output list ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  df_out = list(
    time_of_execution = now(),    # time-stamp
    duration = NULL,              # execution duration
    figs_prefit = NULL,           # figures of data prior to entering the fitting functions
    mout = NULL                   # for each country the model output
  )
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Running selected models ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  mout <- list()
  
  # ---- |-Run model for each country ----
  
  pr=paste("> Model run for all countries completed. \n"); cat(green(pr))
  
  #### output
  t2 <- Sys.time()
  df_out$duration = get_in_hms(t2, t1)
  
  return(df_out)
}

plot_all_season = function(all_season) {
  p = all_season %>% 
    filter(season%in%c(params$SIR_multiseason$seasons_include) ) %>% 
    unnest(respicompass_ili_plus) %>% 
    ggplot(aes(date,value,color=season)) + geom_line() + 
    facet_wrap(~country_short,scales="free_y") 
  return(p)
}