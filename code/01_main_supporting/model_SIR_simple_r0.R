model_SIR_simple_r0 = function( params=NULL, all_season=NULL, target_input=NULL,pop_country, country_short_input, date_v_fit,season){
  
  # ---- |-Filtering ----
  all_season %>% 
    filter(country_short == country_short_input) %>% 
    select(-typing_sentinel,-typing_nonsentinel,-typing_combined) %>% 
    unnest(inc_iliari) %>% 
    filter(target==params$SIR_simple$target) -> all_season_filtered
  # filter:age groups
  all_season_filtered %>% 
    filter(agegroup==params$SIR_simple$agegroup) -> all_season_filtered
  # filter:time
  all_season_filtered %>% 
    filter( date%in%date_v_fit ) %>% 
    select(country_short,season,date,value) %>% 
    mutate(n=1:n()) -> all_season_fit
  
  if (target_input=="ili_typing_sentinel") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_sentinel)) %>% filter(indicator=="positivity") %>% 
      filter( date%in%date_v_fit ) %>% 
      select(date,value_typing=value)
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing/100)) %>% select(-value_typing)
  }
  if (target_input=="ili_typing_all") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_combined)) %>% filter(indicator=="positivity") %>% 
      filter( date%in%date_v_fit ) %>% 
      select(date,value_typing=value_add_narm) %>% 
      mutate(value_typing=ifelse(is.nan(value_typing),NA,value_typing ) )
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing)) %>% select(-value_typing)
  }
  
  # prepare data for stan fit
  all_season_fit %>% ggplot(aes(date,value)) + geom_line()
  data_mock_fit = all_season_fit
  # Projection dates
  data_mock_project = all_season_fit
  data_mock_project$date = data_mock_project$date+365
  
  #mod2 <- cmdstan_model(stan_file='./stan/SIR_simple.stan') # This compiles the script
  stan_list = list(
    n_week_fit = nrow(data_mock_fit),
    severe_obs_fit = as.integer(replace_na(data_mock_fit$value,0)),
    severe_obs_notna = as.integer(!is.na(data_mock_fit$value)),
    n_week_project = nrow(data_mock_project),
    pop = pop_country,
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious
  )
  fit02=rstan::stan(
    file='./stan/SIR_simple_nas.stan',
    #chains=1 ,thin=1,iter=300,
    chains=6 ,thin=6,iter=1500,
    seed=12, cores = getOption("mc.cores", 1L),
    control=list(
      #adapt_delta=0.9,
      #max_treedepth=14
    ),
    data=stan_list
  ) # X mins
  est_Rnull = precis(fit02,pars = c("Rnull_eff"))
  
  mout = tibble(
    country_short = country_short_input,
    season = season,
    Rnull = est_Rnull$result[1],
    Rnull_Rhat = est_Rnull$result[6]
  )
  
  
  if (F){ # plotting to support debugging
    # create table of parameters
    fit02 %>% gather_draws(SIR_ini[state],
                           prop_severe,
                           pop_infect
    ) %>% 
      mean_qi() -> xp; xp
    
    
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
  return(mout)
}