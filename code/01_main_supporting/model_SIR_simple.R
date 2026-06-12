model_SIR_simple = function( params=NULL, dat=NULL, country_short_input, date_v_fit ){
  
  scenario_tag = "A"
  
  # ---- |-Fitting ----
  # prepare data for stan fit
  data_mock = dat %>% 
    filter(country_short == country_short_input, 
           target == params$SIR_simple$target, 
           agegroup == params$SIR_simple$agegroup) # 
  data_mock %<>% filter( date%in%date_v_fit )
  data_mock %>% ggplot(aes(date,value)) + geom_line()
  data_mock_fit = data_mock
  # Projection dates
  data_mock_project = data_mock
  data_mock_project$date = data_mock_project$date+365
  
  #mod2 <- cmdstan_model(stan_file='./stan/SIR_simple.stan') # This compiles the script
  stan_list = list(
    n_week_fit = nrow(data_mock_fit),
    severe_obs_fit = as.integer(data_mock_fit$value),
    n_week_project = nrow(data_mock_project),
    pop = 9e6,
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious
  )
  fit02=rstan::stan(
    file='./stan/SIR_simple.stan',
    chains=8 ,thin=8,iter=400,
    seed=12, cores = getOption("mc.cores", 1L),
    control=list(
      #adapt_delta=0.9,
      #max_treedepth=14
    ),
    data=stan_list
  ) # X mins
  
  # ---- |-Simulations ----
  #  inputs
  #- for each axis a list with available IDs
  
  ## transmission
  df_out = fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>%
    ungroup() %>%
    left_join(data_mock %>% select(date,country_short,agegroup,target) %>% mutate(t_vw = 1:n()), by="t_vw") %>%
    mutate(model = "SIR_simple",
           country_short = country_short,
           agegroup = agegroup,
           target = target,
           scenario_tag = scenario_tag,
           prediction_type = "sample"
    ) %>%
    rename(value = .value,
           sample_or_quantile = .draw) %>%
    select(-.chain, -.iteration, -.variable, -t_vw) 
  df_out %>% select(sample_or_quantile,date,value) %>% group_by(sample_or_quantile) %>% 
    nest() %>% rename(id=sample_or_quantile)  -> df
  list_transmission = df_to_list(df)
  list_transmission %>% names()
  ## vaccine uptake
  list_vaccine_id = list()
  list_vaccine_id[[1]] = list(
    vaccine_uptake = matrix(0*list_transmission[[1]]$value,ncol=1),
    VE_severe = 0.7
  )
  ## severity
  list_severe_id = list()
  list_severe_id[[1]] = list(
    
  )
  ## select needed IDs within each list & create mapping of the IDs across axes -> axis_ids_simulate
  axis_ids_simulate = tibble( round_id="2024-01",
                              transmission_id= c(1:100) ,
                              vaccine_id = 1,
                              severe_id = 1) %>% 
    ungroup() %>% mutate(sim_sample_id=1:n())
  axis_ids_simulate$sim = list(0)
  
  ## loop through axis_ids_simulate 
  for (sim_i in 1:nrow(axis_ids_simulate) ) {
    ## ---- |-Get axis IDs ----
    transmission_id = axis_ids_simulate$transmission_id[sim_i]
    vaccine_id = axis_ids_simulate$vaccine_id[sim_i]
    severe_id = axis_ids_simulate$severe_id[sim_i]
    
    ## ---- |-Prepare the details for each axis ----
    # prepare transmission
    transmission_df = list_transmission[[ transmission_id ]]
    incident_infections = transmission_df %>% select(value) %>% as.matrix() # format: [t,a]
    # prepare vaccination ( using list_vaccine_id )
    vaccine_uptake = list_vaccine_id[[vaccine_id]]$vaccine_uptake
    VE_severe = list_vaccine_id[[vaccine_id]]$VE_severe
    vax_sev = vaccine_severity_nowane(vaccine_uptake,VE_severe)
    # prepare natural severity ( using list_severe_id )
    severity_options = list_severe_id[[severe_id]]
    nat_sev = natural_severity(incident_infections,severity_options)
    severity_baseline = 0.5
    
    # run: transmission-to-severity 
    sev_fact = severity_factor(
      incident_infections,
      severity_baseline,
      severity_factor_vaccines=vax_sev$severity_factor_vaccines,
      severity_factor_natural=nat_sev$severity_factor_natural,
      severity_options
    )
    # run: combine all targets -> mysim
    mysim = combine_all_targets_SIR_simple(date_v=transmission_df$date,
                                           incident_infections,
                                           vaccine_uptake,
                                           incident_severe=sev_fact$incident_severe)
    axis_ids_simulate$sim[sim_i] = nest(mysim)[[1,1]]
    
    
    # Columns in the resulting df_out (as per: https://docs.google.com/document/d/13adcxpPdlDvJM5eiFSkMzlWMTcwsx6lVjY25JA26iS4/edit):
    # model_id
    # round_id ["2024_2025_1_FLU1"]
    # scenario_id ["A","B"], target [allowed targets], location ["DE","FR"] 
    # pop_group ["0-12","13-65"], horizon [week integer], target_end_date [Date string ('YYYY-MM-DD')]
    # output_type ["sample"], output_type_id [string: "1","2","3",...], value [float limited to 2 decimals]
    
    # making up a fictive result as placeholder
    mdf = tibble(
      model_id = "ECDC_lefluflu",
      round_id = params$scenario_round_id,
      scenario_id = scenario_tag,
      target = "inc infection",
      location = country_short_input,
      pop_group = params$SIR_simple$agegroup,
      horizon = c(1,2,3,4),
      target_end_date = today()+(c(1,2,3,4)-1)*7,
      output_type = "sample",
      output_type_id = 1,
      value = c(12,25,12,12)
    )
    
    return(mdf)
  }
  
  
  if (F){ # support debugging
    # create table of parameters
    fit02 %>% gather_draws(SIR_ini[state],
                           prop_severe,
                           pop_infect
    ) %>% 
      mean_qi() -> xp; xp
    
    axis_ids_simulate %>% unnest(cols=sim) %>% 
      ggplot(aes(date,inc_death,group=sim_sample_id)) + geom_line()
    
    round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
    (xp[xp$.variable=="SIR_ini"&xp$state==2,".value"]) %>% logit() %>% round(1) -> mI_ini
    (xp[xp$.variable=="SIR_ini"&xp$state==1,".value"]) %>% round(2) -> mS_ini
    (xp[xp$.variable=="pop_infect",".value"]) %>% round(2) -> mProp_inf
    fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>% 
      mean_qi() %>% left_join(data_mock_fit %>% mutate(t_vw = 1:n()),by="t_vw") %>%
      ggplot(aes(x=t_vw)) + 
      geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
      geom_line(aes(y=.value)) +
      geom_point(aes(y=value),col="black") +
      labs(subtitle = paste("Austria 2022/2023: fit |",
                            "prob_severe:", mprob_severe,"\n",
                            "| S_ini:",mS_ini,
                            "| I_ini:",mI_ini,
                            "| prop inf:",mProp_inf)) -> p_cf0; p_cf0
    
    
  }
  
  return(axis_ids_simulate)
}
