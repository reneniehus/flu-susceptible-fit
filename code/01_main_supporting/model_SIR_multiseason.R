fit_with_stan = function(params,stan_list,mod_path,all_season_fit_wide,country_short_input,stan_mod_file) {
  # run the model fit
  fname = paste0("../Big data/fit",country_short_input,".Rdata")
  duration_fit = "NA:NA:NA"
  if (params$load_earlyfit) {
    load(file = fname)
  } else {
    start_fit <- Sys.time()
    m <- stan_model(file=stan_mod_file) # 
    # rstan:vb settings
    # grad_samples: samples to determine the gradient ( 2 is slower than 5, )
    # tol_rel_obj: default=0.01, smaller means more strict with convergence
    quick_vb = params$rapid_stan_fit
    if (!quick_vb) rstan_vb <- function(...) rstan::vb(...,grad_samples=5, tol_rel_obj = 0.005,output_samples = 300,iter=50000)
    if (!quick_vb) rstan_vb <- function(...) rstan::vb(...,grad_samples=5, tol_rel_obj = 0.010,output_samples = 300,iter=40000)
    if (!quick_vb&country_short_input%in%c("FI","IS") ) rstan_vb <- function(...) rstan::vb(...,grad_samples=10, tol_rel_obj = 0.005,output_samples = 400,iter=90000)
    if ( quick_vb) rstan_vb <- function(...) rstan::vb(...,grad_samples=2, tol_rel_obj = 0.018,output_samples = 300,iter= 10000)
    fit00=rstan_vb(m,algorithm = "meanfield",seed=12,data=stan_list) 
    end_fit <- Sys.time(); duration_fit = get_in_hms(end_fit, start_fit)
    if (!quick_vb) save(fit00,file = fname)
  }
  # plot model fit against fitted data
  p1 = plot_fit(fit00,stan_list,country_short_input)
  p2 = plot_fit_byage(fit00,stan_list,country_short_input)
  # extract data summaries
  season_ili_mean    =sum(stan_list$ili_obs_fit[stan_list$season_id_week_fit%in%c(1,2,3),])/3 # observed burden # 35517.33
  season_ili_mod_mean = fit00 %>% gather_draws(delta_ili_percap_weekly_sum[n]) %>% 
    filter(.draw%in%c(1:500)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(n) %>% summarise(mean_value = mean(.value)) %>% ungroup() %>% 
    right_join( stan_list$all_season_fit,
                by = join_by(n)) %>% group_by(season) %>% summarise(cum_ili=sum(mean_value) %>% log()) %>% 
    pull(cum_ili) %>% mean() %>% exp()
  
  # extract fitted parameters
  df = NULL
  mp="prop_ili_mu";    x=summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary; df=rbind(df,x)
  mp="prop_ili";       x=summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary; df=rbind(df,x[1:stan_list$n_age_groups,])
  mp="ar";             x=summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary; df=rbind(df,x)
  mp="SIR_ini_mu";     x=summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary; df=rbind(df,x)
  mp="cum_ili_log";    x=summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary; df=rbind(df,x)
  mp="reciprocal_phi"; x=summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary; df=rbind(df,x)
  
  if (F) source("code/01_main_supporting/old_stan_fit_code.R") # code of previous fitting implementations
  
  # output list 
  mout=list()
  mout$stan_list = stan_list
  mout$fit = fit00
  mout$modelled_proj = extract_projections(params,fit00,n_iter=300,
                                           stan_list$df_scenarios,
                                           stan_list$df_agegroups,
                                           stan_list$all_season_project) %>% lazy_dt()
  mout$plot_fit = p1
  mout$plot_fit_byage = p2
  mout$pars_df = df
  #
  mout$season_ili_mean = season_ili_mean
  mout$season_ili_mod_mean = season_ili_mod_mean
  mout$duration_fit
  return(mout)
}

# fitting with sequential ABC 
fit_with_eabc = function(params,stan_list,mod_path) {
  
  # Define priors
  myPriors <- list('S01' = c("unif",0,0.5),
                   'S02' = c("unif",0,0.5),
                   'S03' = c("unif",0,0.5))
  myPriors <- list('S01' = c("unif",0,0.5))
  
  if (T){
    tic()
    x<-generate_ili_epi_test( c(1,0.3),stan_list) 
    toc()
  }
  
  # Wrap up model in function that outputs summary stats
  myModel <- function(par){
    
    stan_list_f = generate_ili_epi_test(par,stan_list)
    
    return( stan_list_f$ili_obs_fit$age_1[stan_list_f$ili_obs_notna$age_1==1] ) 
  }
  
  # Define targets 
  myTarget <- c( stan_list$ili_obs_fit$age_1[stan_list$ili_obs_notna$age_1==1] )
  
  # 
  dist_euc <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
  dist_euc(myTarget,myModel( c(1,0.1)))
  
  # Run ABC-SMC (this should be parallelised, see package help)
  library(EasyABC)
  
  rval <- ABC_sequential(method = "Beaumont", 
                         model = myModel, 
                         prior = myPriors, 
                         nb_simul = 5, 
                         summary_stat_target = myTarget,
                         #n_cluster=8,
                         tolerance_tab = c(188521,188520),
                         use_seed=TRUE,
                         progress_bar=T
  )
  
  # Plot posteriors
  hist(rval$param[,1])
  hist(rval$param[,2])
  plot(rval$param[,1], rval$param[,2])
  
}

# computing the data frame from all_season_country
wrangle_fit_df = function(params,data,all_season_country,country_short_input,target_input){
  
  all_season = all_season_country
  
  if (target_input=="erviss_ili_plus") {
    all_season %>% 
      filter(country_short == country_short_input,
             season%in%params$SIR_multiseason$seasons_include) -> x
    
    # to compare the different ili_plus versions
    x %>% select(country_short,season,ili_plus_agesplit_sum,ili_plus_erviss_sent_agesplit_sum,ili_plus_erviss_nonsent_agesplit_sum)
    
    sent = NA
    if (country_short_input %in% params$ili_plus_sentinel   ) sent = T;
    if (country_short_input %in% params$ili_plus_nonsentinel) sent = F;
    if ( sent) {x %>% unnest(erviss_ili_plus_sentinel   ) -> y; y %>% ggplot(aes(date,value)) +geom_line()+labs(subtitle='erviss_sent')}
    if (!sent) {x %>% unnest(erviss_ili_plus_nonsentinel) -> y; y %>% ggplot(aes(date,value)) +geom_line()+labs(subtitle='erviss_nonsent')}
    if (is.na(sent)) warning("Country has unclear sentinal/nonsentinel ili-plus indicator")
    
    y %>% select(country_short,date,season,agegroup,value)-> all_season_fit
    
    # rescale countries where ILI have diff denominator in ERVISS
    if ( country_short_input%in%params$ili_diff_denom_erviss ) all_season_fit$value=1000*all_season_fit$value 
    
    # impute summer low-activity
    all_season_fit %>% mutate(
      summer_low_day = as.integer(date%in%params$summer_low_dates),
      value=ifelse( date%in%params$summer_low_dates & is.na(value), 0 , value )
    ) -> all_season_fit
    
  }
  
  if (target_input=="respicompass_ili_plus") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      unnest(respicompass_ili_plus) %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value)-> all_season_fit
    # impute summer low-activity
    all_season_fit %>% mutate(
      summer_low_day = as.integer(date%in%params$summer_low_dates),
      value=ifelse( date%in%params$summer_low_dates & is.na(value), 0 , value )
    ) -> all_season_fit
  }
  
  if (target_input=="ili") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      unnest(inc_iliari) %>% 
      filter(target=="ILIconsultationrate") %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value) -> all_season_fit
  }
  
  if (target_input=="own_ili_plus") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      unnest(inc_iliari) %>% 
      filter(target=="ILIconsultationrate") %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value) -> all_season_ili
    
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_combined)) %>% filter(indicator=="positivity") %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value_typing=value_add_narm) %>% 
      mutate(value_typing=ifelse(is.nan(value_typing),NA,value_typing ) )
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing)) %>% select(-value_typing)
  }
  
  # take care of age groups
  all_season_fit_wide = all_season_fit %>% 
    pivot_wider(names_from = agegroup, values_from = value) %>% 
    mutate(n=1:n())
  
  return(all_season_fit_wide)
} 

# computing a list with all input required by the model
make_stan_list = function(params,data,all_season_fit_wide,country_short_input,vax_country,pop_country,contacts,age_collapse="all"){
  # helpers for the fit and project dataframes
  start_year = params$proj_start_year
  season     = paste0(start_year,"/",start_year+1)
  start_date = ymd(paste0(start_year,params$season_start_monthday))
  end_date   = ymd(paste0(start_year+1,params$season_end_monthday))
  date_v = seq(from=start_date,to=end_date,by="day")
  date_v_wed = date_v[weekdays(date_v)=="Wednesday"][1:52]
  date_v_mon = date_v[weekdays(date_v)=="Monday"   ][1:52]
  
  # dataframe for projections
  all_season_project = tibble(country_short=country_short_input,
                              season=season,
                              date_mon=date_v_mon,
                              date_wed=date_v_wed,
                              value=NA) %>% 
    mutate(week_id=1:n()) %>% 
    # adding the RespiCompass indicators for projection horizons
    left_join(data$helpers_respicompass$iso_weeks,by=c("date_mon"="start_week_day"))
  
  # make daily version of the data frame - from some daily indicators that the model needs
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_fit_wide$date,season=all_season_fit_wide$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) %>% 
    group_by(season) %>%  mutate(h=1:n(),
                                 season_start=case_when(h==1~1,h==2~2,.default=0) ) %>% 
    ungroup() %>% select(-h) -> all_season_fit_daily
  
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_project$date_mon,season=all_season_project$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) %>% ungroup() %>% 
    group_by(season) %>%  mutate(h=1:n(),
                                 season_start=case_when(h==1~1,h==2~2,.default=0) ) %>% 
    ungroup() -> all_season_project_daily
 
  # plotting
  p1 = all_season_fit_wide %>% ggplot(aes(date,age_total)) + geom_line() + geom_rug()
  
  # helpers for stan list
  df_scenarios = params$scenarios
  pop_pyramid = data$demography_respicast$population_pyramid %>% filter(country==EU_long(country_short_input))
  pop_pyramid = pop_pyramid %>% select(age_group,population) %>% deframe()
  pop_age_group = pop_pyramid[params$four_age_groups]
  age_groups = params$SIR_multiseason$age_groups
  n_age_groups = length(age_groups)
  df_agegroups_ecdc = tibble(agegroup_id=c(0:n_age_groups),
                             age_group_ecdc=c("age_total", age_groups) )
  df_age_translate = tibble(age_group_ecdc=df_agegroups_ecdc$age_group_ecdc,
                            age_group_respicompass=c("total","0-4","5-14","15-64","65+"))
  df_agegroups = df_agegroups_ecdc %>% left_join(df_age_translate,by = join_by(age_group_ecdc))
  z_proj = rep(0,nrow(all_season_project_daily))
  z_fit  = rep(0,nrow(all_season_fit_daily))
  
  contact_matrix = contacts[[EU_long(country_short_input)]]
  if (country_short_input %in% c("NO","CZ")) {
    contact_matrix = contacts[["EU"]]
  }
  # first compute the total number of contacts made by a member of age group i with individuals across all age groups
  contact_activity_age = rowSums(contact_matrix)
  # then compute the activity factor per age group relative to the population weighted average number of contacts of all age groups
  a_factor = contact_activity_age/weighted.mean(contact_activity_age,w = pop_age_group)
  # normalise such that each row i contact_matrix[i,] gives the relative distribution of contact of any invidual of age group i
  for (a in 1:n_age_groups) {
    contact_matrix[a,] = contact_matrix[a,]/sum(contact_matrix[a,])
  }
  
  # ---- |-Stan list and fit----
  stan_list = list(
    ## EXTRA stuff good to carry forward
    all_season_fit=all_season_fit_wide,
    all_season_project=all_season_project,
    season_id_raw = fct_inorder(all_season_fit_daily$season) %>% levels() %>% enframe(),
    df_scenarios = df_scenarios,
    df_agegroups = df_agegroups,
    ili_fit_date = all_season_fit_wide$date,
    ili_proj_date = all_season_project$date_wed,
    ## data relevated for the fit
    n_season = n_distinct(all_season_fit_wide$season),
    n_season_cum_fit = params$n_season_cum_fit, # as per RespiCompass Round 1 
    n_week_fit = nrow(all_season_fit_wide),
    n_day_fit = nrow(all_season_fit_daily),
    n_age_groups = n_age_groups,
    #
    ili_obs_fit = all_season_fit_wide %>% 
      select( any_of(params$SIR_multiseason$age_groups) ) %>% 
      mutate_all(~ replace_na(.,0) ),
    ili_obs_notna = all_season_fit_wide %>% 
      select( any_of(params$SIR_multiseason$age_groups) ) %>% 
      mutate_all(~ !is.na(.) ) %>% mutate_all(~as.integer(.) ),
    ili_summer_low = all_season_fit_wide$summer_low_day,
    #
    season_start_fit = as.integer(all_season_fit_daily$season_start),
    season_id_day_fit = fct_inorder(all_season_fit_daily$season) %>% as.integer(),
    season_id_week_fit = fct_inorder(all_season_fit_wide$season) %>% as.integer(),
    #
    pop = sum(pop_age_group), 
    #
    pop_age_group=matrix(pop_age_group ,nrow=n_age_groups,ncol=1),
    contact_matrix=contact_matrix,
    a_factor=a_factor,
    delta_vax=tibble( A=z_fit,B=z_fit,C=z_fit,D=z_fit) %>% mnaming(age_groups),
    # data relevant for projected scenarios
    n_week_proj = nrow(all_season_project),
    n_day_proj= nrow(all_season_project)*7,
    season_start_proj = as.integer(all_season_project_daily$season_start),
    season_id_day_proj = fct_inorder(all_season_project_daily$season) %>% as.integer(),
    season_id_week_proj = fct_inorder(all_season_project_daily$season) %>% as.integer(),
    n_scenario = nrow(df_scenarios),
    axis_transmission = df_scenarios$axis_transmission,
    axis_vax = df_scenarios$axis_vax,
    delta_vax_real=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    delta_vax_opti=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    delta_vax_pess=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    delta_vax_null=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    # epi parameters
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious,
    ve_spread = params$ve_spread,
    ve_inf = params$ve_inf,
    ve_ili_cond_inf = params$ve_ili_cond_inf,
    # daily steps
    n_daily_time_steps = 1,
    # priors
    sigma_cum_ili = 2,
    prior_sigma_prop_ili = 2,
    prior_sigma_i = 5,
    prior_sigma_s = 2,
    sigma_prop_ili_age = 2
  )
  # transform ILI rates into incidence
  stan_list$ili_obs_incs = stan_list$ili_obs_fit
  for (a in 1:stan_list$n_age_groups) {
    stan_list$ili_obs_incs[,a] = (stan_list$ili_obs_fit[,a]*stan_list$pop_age_group[1,]/100000) %>% mutate_all(~as.integer(.) )
  }
  stan_list$age_total_incs = rowSums(stan_list$ili_obs_incs)
  # Add vaccination to Oct 1st to the oldest age group
  # for fitting
  my_date_v = all_season_fit_daily$date
  ind_vax = ( month(my_date_v)==10 & day(my_date_v)==1 )
  # use historical vaccine coverage
  hist_vax_ind = tibble(iso2_code=country_short_input, 
         season=all_season_fit_daily$season[ind_vax],
         target_group="65+y") %>% 
    left_join(data$vax$data_vax_history,join_by(iso2_code, season, target_group))
  hist_vax_ind$vaccine_coverage[is.na(hist_vax_ind$vaccine_coverage)] = (vax_country$higher_vax_coverage + vax_country$lower_vax_coverage)/2
  stan_list$delta_vax$age_65_99[ind_vax] = hist_vax_ind$vaccine_coverage
  # for projections
  my_date_v = all_season_project_daily$date
  ind_vax = which( my_date_v == paste0( year(min(my_date_v)),"-10-01" ) )
  stan_list$delta_vax_real$age_65_99[ind_vax] = (vax_country$higher_vax_coverage + vax_country$lower_vax_coverage)/2
  stan_list$delta_vax_opti$age_65_99[ind_vax] = vax_country$higher_vax_coverage
  stan_list$delta_vax_pess$age_65_99[ind_vax] = vax_country$lower_vax_coverage
  stan_list$delta_vax_null$age_65_99[ind_vax] = 0
  
  stan_list$daily_counter_fit = rep(1:stan_list$n_day_fit, each=stan_list$n_daily_time_steps)
  stan_list$daily_counter_proj = rep(1:stan_list$n_day_proj, each=stan_list$n_daily_time_steps)
  
  stan_list$daily_daystart_fit = rep(1:stan_list$n_daily_time_steps, each=stan_list$n_day_fit)
  stan_list$daily_daystart_proj = rep(1:stan_list$n_daily_time_steps, each=stan_list$n_day_proj)
  
  # summary targets
  stan_list$cum_ili_obs_log     = rowsum(x=stan_list$ili_obs_fit,group=stan_list$season_id_week_fit,na.rm = T) %>% rowSums() %>% zero_plus_eps(eps=1/10^6) %>% log()
  stan_list$cum_ili_obs_age_log = rowsum(x=stan_list$ili_obs_fit,group=stan_list$season_id_week_fit,na.rm = T) %>% zero_plus_eps(eps=1/10^6) %>% log()
  stan_list$cum_ili_obs_incs_age_log = rowsum(x=stan_list$ili_obs_incs,group=stan_list$season_id_week_fit,na.rm = T) %>% zero_plus_eps(eps=1/10^6) %>% log()
  stan_list$n_ili_obs_notna = rowsum(x=stan_list$ili_obs_notna,group=stan_list$season_id_week_fit,na.rm = T) %>% rowSums()
  stan_list$weight_obs_epi =  stan_list$ili_obs_fit*0 + params$weight_obs_epi #1/mean( stan_list$n_ili_obs_notna ) 
  # 
  stan_list$weight_cum_ili =  1.0 # 
  # support fit of AT
  # if (country_short_input %in% c("AT","IT")) {
  #   stan_list$weight_obs_epi[stan_list$ili_summer_low==0,] = 0 # don't get influenced by seasonal surveillance
  # }
  
  
  ###################################################################
  
  ### for debugging: make it 2 age groups: <65 and above
  if (age_collapse=="two") {
    stan_list$n_age_groups = 2
    stan_list$contact_matrix = matrix(data=c(1/2,1/2,1/2,1/2),nrow=2,ncol=2)
    stan_list$pop_age_group = matrix(data=c(sum(stan_list$pop_age_group[1:3,1]),stan_list$pop_age_group[4,1]),
                                     nrow=2,ncol=1)
    stan_list$ili_obs_fit =  stan_list$ili_obs_fit %>% transmute(
      age_1=replace_na(age_00_04+age_05_14+age_15_64,0) %>% as.integer(),
      age_2=replace_na(age_65_99,0) %>% as.integer()
    )
    stan_list$ili_obs_notna = stan_list$ili_obs_notna %>% transmute(
      age_1=as.integer((age_00_04+age_05_14+age_15_64)==3),
      age_2=as.integer((age_65_99)==1)
    )
    stan_list$delta_vax_real = stan_list$delta_vax_real %>% 
      transmute(age_1=rowSums(across(age_00_04:age_15_64)), age_2=age_65_99)
    stan_list$delta_vax_opti = stan_list$delta_vax_opti %>% 
      transmute(age_1=rowSums(across(age_00_04:age_15_64)), age_2=age_65_99)
    stan_list$delta_vax_pess = stan_list$delta_vax_pess %>% 
      transmute(age_1=rowSums(across(age_00_04:age_15_64)), age_2=age_65_99)
    stan_list$delta_vax_null = stan_list$delta_vax_null %>% 
      transmute(age_1=rowSums(across(age_00_04:age_15_64)), age_2=age_65_99)
    stan_list$delta_vax = stan_list$delta_vax %>% 
      transmute(age_1=rowSums(across(age_00_04:age_15_64)), age_2=age_65_99)
  }
  ### for debugging: make it 1 age group
  if (age_collapse=="one") {
    stan_list$n_age_groups = 1
    stan_list$contact_matrix = matrix(data=c(1),nrow=1,ncol=1)
    stan_list$pop_age_group = matrix(c(pop_country) ,nrow=1,ncol=1)
    stan_list$ili_obs_notna =  stan_list$ili_obs_notna %>% transmute(
      age_1=as.integer( (age_00_04+age_05_14+age_15_64+age_65_99) == 4)
    )
    stan_list$ili_obs_fit =  stan_list$ili_obs_fit %>% transmute(
      age_1= age_00_04+age_05_14+age_15_64+age_65_99
    )
    stan_list$delta_vax_real = stan_list$delta_vax_real %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax_opti = stan_list$delta_vax_opti %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax_pess = stan_list$delta_vax_pess %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax_null = stan_list$delta_vax_null %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax = stan_list$delta_vax %>% select(age_65_99) %>% rename(age_1=age_65_99)
  }
  return(stan_list)
}

# extract projections from the model in the format required by RespiCompass
extract_projections = function(params,fit00,n_iter,df_scenarios,df_agegroups,all_season_project){
  # extract projections
  modelled_proj = fit00 %>% 
    gather_draws(gen_ili_u_percap_obs_proj[scen_id,week_id,agegroup_id], # if you want to apply changes here, do look up the useful gather_draws {tidybayes} syntax
                 gen_ili_v_percap_obs_proj[scen_id,week_id,agegroup_id],
                 gen_ili_t_percap_obs_proj[scen_id,week_id,agegroup_id],
                 gen_ili_u_percap_obs_proj_sum[scen_id,week_id],
                 gen_ili_v_percap_obs_proj_sum[scen_id,week_id],
                 gen_ili_t_percap_obs_proj_sum[scen_id,week_id]) %>% 
    filter(.draw%in%c(1:n_iter)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% # remove unneeded columns and grouping
    mutate(vax_status=case_when(
      .variable=="gen_ili_u_percap_obs_proj"~"vaxNo",
      .variable=="gen_ili_v_percap_obs_proj"~"vaxYes",
      .variable=="gen_ili_t_percap_obs_proj"~"vaxTotal",
      .variable=="gen_ili_u_percap_obs_proj_sum"~"vaxNo",
      .variable=="gen_ili_v_percap_obs_proj_sum"~"vaxYes",
      .variable=="gen_ili_t_percap_obs_proj_sum"~"vaxTotal"
    )) %>% 
    mutate(agegroup_id=case_when(
      .variable=="gen_ili_u_percap_obs_proj_sum"~0,
      .variable=="gen_ili_v_percap_obs_proj_sum"~0,
      .variable=="gen_ili_t_percap_obs_proj_sum"~0,
      TRUE ~ agegroup_id
    )) %>% 
    left_join(df_scenarios,by = join_by(scen_id)) %>% # add scenario info
    left_join(df_agegroups,by = join_by(agegroup_id)) %>% # add agegroup info
    left_join(all_season_project, by=join_by(week_id)) %>% # add week info
    mutate(model_id=params$scenario_model,round_id=params$scenario_round_id, # add needed columns
           target="ili_plus",output_type="sample") %>% 
    unite(col="pop_group",age_group_respicompass,vax_status,sep="_") %>% 
    select( round_id, # select according to submission definition
            scenario_id,
            target=target,
            location=country_short,
            pop_group=pop_group,
            horizon=horizon,
            target_end_date=end_week_day,
            date_mon=date_mon, # keeping this in to have dates for all weeks beyond data$helpers_respicompass$iso_weeks
            output_type=output_type,
            output_type_id=.draw,
            value=.value) 
  # Required columns in df_for_submission (as per: https://github.com/european-modelling-hubs/RespiCompass/wiki/Submission-format):
  return(modelled_proj)
}

plot_fit = function(fit00,stan_list,country_short_input) {
  # ili rate 
  modelled_fit = fit00 %>% gather_draws(gen_ili_obs_percap_fit_sum[n]) %>% 
    filter(.draw%in%c(1:500)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(n) %>% summarise(mean_value = mean(.value),
                              low=quantile(.value,probs=0.1),
                              upp=quantile(.value,probs=0.9)) %>% ungroup() %>% 
    right_join( stan_list$all_season_fit,
                by = join_by(n)) 
  modelled_proj = fit00 %>% gather_draws(gen_ili_t_percap_obs_proj_sum[scen,week_id]) %>% 
    filter(.draw%in%c(1:500)) %>%
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(scen,week_id) %>% summarise(mean_value = mean(.value),
                                         low=quantile(.value,probs=0.1),
                                         upp=quantile(.value,probs=0.9)) %>% ungroup() %>% 
    right_join( stan_list$all_season_project ,
                by = join_by(week_id)) %>% mutate(date=date_mon)
  p1 = modelled_fit %>% ggplot() + geom_line(aes(date,age_total)) + 
    geom_line(data=. %>% filter(season=="2024/2025" ),aes(date,age_total),col="red") +
    geom_ribbon(aes(x=date,ymin=low,ymax=upp),fill="lightblue") +
    geom_ribbon(data=modelled_proj,aes(x=date,fill=as.factor(scen),ymin=low,ymax=upp)) + guides(fill="none")+
    labs(x="",y="",subtitle = paste0( EU_long(country_short_input)," (",country_short_input,") ILI rate") ); p1
  
  # ili incidence
  modelled_fit = fit00 %>% gather_draws(delta_ili_incs_weekly_sum[n]) %>% 
    filter(.draw%in%c(1:500)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(n) %>% summarise(mean_value = mean(.value),
                              low=quantile(.value,probs=0.1),
                              upp=quantile(.value,probs=0.9)) %>% ungroup() %>% 
    mutate( age_total=stan_list$age_total_incs, 
            date=stan_list$ili_fit_date,
            name=stan_list$season_id_week_fit) %>% 
    left_join(stan_list$season_id_raw %>% rename(season=value),by="name")
  
  modelled_proj = fit00 %>% gather_draws(gen_ili_t_obs_proj_sum[scen,week_id]) %>% 
    filter(.draw%in%c(1:500)) %>%
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(scen,week_id) %>% summarise(mean_value = mean(.value),
                                         low=quantile(.value,probs=0.1),
                                         upp=quantile(.value,probs=0.9)) %>% ungroup() %>% 
    right_join( stan_list$all_season_project ,
                by = join_by(week_id)) %>% mutate(date=date_mon)
  p2 = modelled_fit %>% ggplot() + geom_line(aes(date,age_total)) + 
    geom_line(data=. %>% filter(season=="2024/2025"),aes(date,age_total),col="red") +
    geom_ribbon(aes(x=date,ymin=low,ymax=upp),fill="lightblue") +
    geom_ribbon(data=modelled_proj,aes(x=date,fill=as.factor(scen),ymin=low,ymax=upp)) + guides(fill=guide_legend("Scenario")) +
    labs(x="",y="",subtitle = paste0( EU_long(country_short_input)," (",country_short_input,") ILI incidence") ); p2
  p=(p1+p2); p
  return(pp)
}

plot_fit_byage = function(fit00,stan_list,country_short_input){
  # age split version
  
  # ili rate 
  obs_data = stan_list$all_season_fit %>% select(age_00_04:age_65_99,date,season,n) %>% 
    pivot_longer(cols=1:stan_list$n_age_groups,names_to = "age") %>% 
    mutate(age=as.factor(age) %>% as.numeric())
  modelled_fit = fit00 %>% gather_draws(gen_ili_percap_obs_fit[n,age]) %>% 
    filter(.draw%in%c(1:100)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(n,age) %>% summarise(mean_value = mean(.value)) %>% ungroup() %>% 
    left_join(obs_data,by=c("n","age"))
  
  proj_df = stan_list$all_season_project %>% select(country_short,week_id,date=date_wed)
  modelled_proj = fit00 %>% gather_draws(gen_ili_t_percap_obs_proj[scen,week_id,age]) %>% 
    filter(.draw%in%c(1:100)) %>%
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(week_id,age,scen) %>% summarise(mean_value = mean(.value)) %>% ungroup() %>% 
    left_join( proj_df, by=join_by(week_id)  )
  
  p1 = modelled_fit %>% ggplot(aes(date,value)) + geom_line() + 
    geom_line(data=. %>% filter(season=="2024/2025"),aes(date,value),col="red",linewidth=1.4) +
    geom_line(aes(y=mean_value),col="lightblue",alpha=0.8) +
    geom_line(data=modelled_proj,aes(col=as.factor(scen),y=mean_value),alpha=0.5) + guides(col="none") +
    facet_wrap(~age,ncol=1,scales="free_y") +
    labs(x="",y="",subtitle = paste0( EU_long(country_short_input)," (",country_short_input,"), ILI rate") ); p1
  
  # ili incidence
  obs_data = cbind(stan_list$ili_obs_incs , stan_list$all_season_fit %>% select(date,season,n) ) %>% as_tibble() %>% 
    pivot_longer(cols=1:stan_list$n_age_groups,names_to = "age") %>% 
    mutate(age=as.factor(age) %>% as.numeric())
  modelled_fit = fit00 %>% gather_draws(gen_ili_obs_fit[n,age]) %>% 
    filter(.draw%in%c(1:100)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(n,age) %>% summarise(mean_value = mean(.value)) %>% ungroup() %>% 
    left_join(obs_data,by=c("n","age"))
  
  proj_df = stan_list$all_season_project %>% select(country_short,week_id,date=date_wed)
  modelled_proj = fit00 %>% gather_draws(gen_ili_t_obs_proj[scen,week_id,age]) %>% 
    filter(.draw%in%c(1:100)) %>%
    select(-.chain,-.iteration) %>% ungroup() %>% 
    group_by(week_id,age,scen) %>% summarise(mean_value = mean(.value)) %>% ungroup() %>% 
    left_join( proj_df, by=join_by(week_id)  )
  
  p2 = modelled_fit %>% ggplot(aes(date,value)) + geom_line() + 
    geom_line(data=. %>% filter(season=="2024/2025"),aes(date,value),col="red",linewidth=1.4) +
    geom_line(aes(y=mean_value),col="lightblue",alpha=0.8) +
    geom_line(data=modelled_proj,aes(col=as.factor(scen),y=mean_value),alpha=0.5) + guides(col=guide_legend("Scenario")) +
    facet_wrap(~age,ncol=1,scales="free_y") +
    labs(x="",y="",subtitle = paste0( EU_long(country_short_input)," (",country_short_input,"), ILI incidence") ); p2
  p = (p1 + p2); p
  
  return(p)
}
