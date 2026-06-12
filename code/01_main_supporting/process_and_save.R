process_and_save = function(params=NULL, data=NULL, models_out=NULL,save_submission){
  pr=paste("Processing and saving ... \n"); cat(green(pr))
  
  ## ---- |-Put all countries together ----
  df_para = NULL
  df_submission = NULL
  df_data_summaries = NULL
  for (i in 1:length(models_out$mout)) {
    i_country = names(models_out$mout)[i]
    i_country_long = names(models_out$mout)[i] %>% EU_long()
    pr=paste(i_country_long,"\n"); cat(green(pr))
    # parameter estimates
    xpar = models_out$mout[[i]]$pars_df
    xpar = as_tibble(xpar, rownames = "para")
    xpar$country = i_country
    xpar$country_long = i_country_long
    df_para = rbind(df_para,xpar)
    # submission df
    xsubm = models_out$mout[[i]]$modelled_proj
    df_submission = rbind(df_submission,as_tibble(xsubm))
    # projected burden
    x1 = models_out$mout[[i]]$season_ili_mod_mean
    
    . %>% filter(pop_group=="total_vaxTotal") %>% 
      group_by(output_type_id,scenario_id) %>% summarise( value=sum(value) ) %>% ungroup()  -> mfu
    
    baseline_low = xsubm %>% filter(scenario_id=="G") %>% mfu() %>% summarise(val=quantile(value,0.1)) %>% as_tibble() %>% deframe()
    baseline_upp = xsubm %>% filter(scenario_id=="G") %>% mfu() %>% summarise(val=quantile(value,0.9)) %>% as_tibble() %>% deframe()
    pessimis_low = xsubm %>% filter(scenario_id%in%c("B","D","F") ) %>% mfu() %>% summarise(val=quantile(value,0.1)) %>% as_tibble() %>% deframe()
    pessimis_upp = xsubm %>% filter(scenario_id%in%c("B","D","F") ) %>% mfu() %>% summarise(val=quantile(value,0.9)) %>% as_tibble() %>% deframe()
    x_df = tibble(location=i_country_long,
                  baseline_low=baseline_low/x1,
                  baseline_upp=baseline_upp/x1,
                  pessimis_low=pessimis_low/x1,
                  pessimis_upp=pessimis_upp/x1)
    df_data_summaries = rbind(df_data_summaries,x_df)
    # 
  }
  ## ---- |-Sense checks ----
  # total fit + proj
  pdf(width = 12, height = 8, "code/03_report/fit_flip_tot.pdf")
  for (i in 1:length(models_out$mout)) {
    # figures
    models_out$mout[[i]]$plot_fit %>% print()
  }
  dev.off()
  # by-age fit + proj
  pdf(width = 12, height = 8, "code/03_report/fit_flip.pdf")
  for (i in 1:length(models_out$mout)) {
    # figures
    models_out$mout[[i]]$plot_fit_byage %>% print()
  }
  dev.off()
  
  # explore parameters
  ppar = list()
  df_para %>% filter(para=="prop_ili_mu") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) + 
    scale_x_log10() + labs(subtitle="prop_ili_mu") -> ppar$p1;ppar$p1
  df_para %>% filter(para=="SIR_ini_mu[1]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`))  + labs(subtitle="S_ini")-> ppar$p2; ppar$p2
  df_para %>% filter(para=="SIR_ini_mu[2]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) +
    scale_x_log10() + labs(subtitle="I_ini") -> ppar$p3; ppar$p3
  df_para %>% filter(para=="reciprocal_phi") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) +
    scale_x_log10() + labs(subtitle="reciprocal_phi") -> ppar$p4; ppar$p4
  # explore submissions: ordering
  df_submission %>% filter(pop_group=="total_vaxTotal") %>% group_by(scenario_id) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(cum_burden_log)
  df_submission %>% filter(pop_group=="total_vaxTotal") %>% group_by(scenario_id,location) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% ungroup() %>% arrange(location,cum_burden_log)
  #
  x = df_submission %>% filter(pop_group=="total_vaxTotal")  %>% group_by(scenario_id,location) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(location,cum_burden_log)
  
  # explore vax versus no-vax for individual countries
  # mcountry="RO"
  # x %>% filter(location==mcountry) %>% pull(cum_burden_log) -> mburd
  # names(mburd) = x %>% filter(location==mcountry) %>% pull(scenario_id)
  # mburd = exp(mburd)
  # (mburd["E"]-mburd["G"])/mburd["E"]
  
  ## ---- |-Save ----
  if (save_submission) {
    library(arrow)
    sub_cols = c("round_id","scenario_id","target","location","pop_group",
                 "horizon","target_end_date","output_type","output_type_id","value")
    df_submission %>% 
      select(any_of(sub_cols) ) %>% filter(output_type_id%in%c(1:300)) %>% 
      mutate(horizon=as.integer(horizon),output_type_id=as.character(output_type_id)) %>% 
      filter(scenario_id%in%c("A","B","C","D","E","F")) %>% 
      filter(!is.na(horizon)) -> x
    el_distinct = c("round_id"=1,"scenario_id"=6,"target"=1,"location"=18,
                    "pop_group"=15,"horizon"=43,"output_type_id"=300)
    cumprod(el_distinct)
    
    x %>% filter(value>60000,location=="MT") %>% as_tibble() %>% describe()
    
    x %>% write_parquet("../Big data/2024_2025_1_FLU-ECDC-flumod.parquet")
  }
  
  
  # pack the list
  rep_list = list(
    N_countries_fit=models_out$mout %>% length(),
    df_data_summaries=df_data_summaries,
    df_submission=df_submission,
    # 
    df_para=df_para,
    ppar=ppar
  )
  #save(rep_list,file="./output/rep_list.Rdata") # comment as this file is too large for github
  if (save_submission) save(data,params,rep_list,models_out,file="../Big data/RespiCompass_round1_models_out.Rdata") # ca 13 MBs
  
  return(rep_list)
}