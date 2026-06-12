// A model that fits to several seasons, and is structured by age and vaccine status

// Sentinel flu data is reported in rates, while the model produces incidences
// incs: denotes incidences in absolute numbers
// rate: denotes rates which is incidence/pop-size * 100000 

data {
  // data relevant for the fit 
  int n_age_groups;  // number of age groups
  int n_season;      // number of seasons
  int n_season_cum_fit; // number of seasons for the cumulative fit
  int n_week_fit;    // number of observable values, weekly
  int n_day_fit;     // number of obervatble values, daily
  real ili_obs_fit[n_week_fit, n_age_groups]; // observed ili
  int ili_obs_incs[n_week_fit, n_age_groups]; // observed ili as incidence
  real cum_ili_obs_log[ n_season ]; // observed cumulative ili rate (log-scale) by season
  real cum_ili_obs_age_log[ n_season,n_age_groups ]; // observed cumulative ili rate (log-scale) by season and age-group
  real cum_ili_obs_incs_age_log[ n_season,n_age_groups ]; // observed cumulative ili incidence (log-scale) by season and age-group
  int n_daily_time_steps; // number of daily steps
  array[n_week_fit,n_age_groups]int<lower=0,upper=1> ili_obs_notna; // indicating non-missing data with 1, otherwise 0
  array[n_day_fit] int<lower=0,upper=2> season_start_fit; // indicating first week of a season with 1, the second week with 2, otherwise 0
  array[n_day_fit] int<lower=1,upper=n_season> season_id_day_fit; // indicating which seasn each obervable day belongs to
  array[n_week_fit] int<lower=1,upper=n_season> season_id_week_fit; // indicating which seasn each obervable day belongs to
  int pop; // total population size
  array[n_age_groups,1] int pop_age_group; // population size per age group, required to be a matrix 
  matrix[n_age_groups, n_age_groups] contact_matrix; // contact matrix
  real a_factor[n_age_groups]; // age-specific modulation of beta, informed by contact matrix
  matrix[n_day_fit, n_age_groups] delta_vax; // daily fraction of newly vaccinated individuals per age group
  array[n_day_fit*n_daily_time_steps] int<lower=1,upper=n_day_fit> daily_counter_fit;
  array[n_day_fit*n_daily_time_steps] int<lower=1,upper=n_daily_time_steps> daily_daystart_fit;
  // data relevant for projected scenarios
  int n_week_proj; // number of projected weeks
  int n_day_proj; // number of projected days
  int n_scenario;// // number of projected scenarios
  int axis_transmission[n_scenario]; // indicator for the transmission scenario axis
  int axis_vax[n_scenario]; // indicator for the vaccine scenario axis
  matrix[n_day_proj, n_age_groups] delta_vax_real; // daily assumed vax uptake in projection period
  matrix[n_day_proj, n_age_groups] delta_vax_opti; // daily assumed vax uptake in projection period
  matrix[n_day_proj, n_age_groups] delta_vax_pess; // daily assumed vax uptake in projection period
  matrix[n_day_proj, n_age_groups] delta_vax_null; // daily assumed vax uptake in projection period
  array[n_day_proj] int<lower=0,upper=2> season_start_proj; // indicating first week of a season with 1, the second week with 2, otherwise 0
  array[n_day_proj*n_daily_time_steps] int<lower=1,upper=n_day_fit> daily_counter_proj;
  array[n_day_proj*n_daily_time_steps] int<lower=1,upper=n_daily_time_steps> daily_daystart_proj;
  array[n_day_proj] int<lower=1,upper=n_season> season_id_day_proj; // indicating which seasn each projected day belongs to
  // epi parameters
  real Rnull; // R0
  real rate_infectious; // infectious rate, such that beta = Rnull*rate_infectious
  // (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf)
  real ve_spread; // vaccine effectiveness on onward transmission/infectiousness
  real ve_inf; // vaccine effectiveness on susceptability
  real ve_ili_cond_inf; // vaccine effectiveness on severity, given infection
  //
  real<lower=0> sigma_cum_ili; // variability of observed cumulative ili (log-scale)
  real<lower=0> prior_sigma_prop_ili;
  real<lower=0> prior_sigma_i; 
  real<lower=0> prior_sigma_s;
  real<lower=0> sigma_prop_ili_age;
  //
  array[n_week_fit,n_age_groups] real weight_obs_epi; // give a weight to the observed epi likelihood
  real weight_cum_ili; // give a weight to the observed cumulative burden likelihood
}

transformed data {
  // epi parameters
  real beta = rate_infectious * Rnull;
  // time step manipulations
  int  n_multi_day_proj = n_day_proj * n_daily_time_steps; // time steps for projections = number of days for projections x sub-daily steps
  int  n_multi_day_fit = n_day_fit * n_daily_time_steps; // time steps for fit = number of days for fitting period x sub-daily steps
  real dt = 1.0/n_daily_time_steps;
}

parameters {
  // real mock_para; // can be actived when all other parameters are fitted
  // note: a simplex of length 3, has 2 free parameters
  // note: simplex[n] X[m,o] creates an m x o sized array of simplex, each of size n
  
  // initial state
  simplex[3] SIR_ini_mu; // overall season mean
  // real i_season[n_season];
  // real r_season[n_season];
  
  // ratio between observed and immunising
  real<lower=0> prop_ili_mu; // overall mean over season
  
  // real prop_ili_season[n_season];
  real prop_ili_age[n_age_groups];
  
  // dispersion parameters
  // real<lower=0> sigma_prop_ili_age;
  // real<lower=0> sigma_prop_ili_season;
  // real<lower=0> sigma_i;
  // variability of the incidence data
  real<lower=0, upper=1> reciprocal_phi; // overdipersion parameter for ili obs fit, var=mu+reciprocal_phi*mu^2

}

transformed parameters {
  //
  
  // 
  matrix<lower=0, upper=1>[n_season,n_age_groups] ar; // attack rate
  
  array[n_week_fit,n_age_groups] real<lower=0> delta_ili_incs_weekly; // ili/detectable incidence in absolute numbers, weekly aggregate
  array[n_week_fit,n_age_groups] real<lower=0> delta_ili_percap_weekly; // per 100000 of age group
  array[n_season,n_age_groups] real cum_ili_log ; // to store the sum of ili for each season
  real phi; // dispersion parameter of the observeation process, var=mu+reciprocal_phi*mu^2
  simplex[3] SIR_ini[n_season, n_age_groups]; // how to access: SIR_ini[season,age,compartment] // S I R initial values per season, 1 can be replaced by n_age_groups
  real<lower=0> prop_ili[n_season, n_age_groups]; // proportion of infections that are ili 
  phi = 1 / reciprocal_phi; // dispersion parameter: var=mu+reciprocal_phi*mu^2
  for (s in 1:n_season) {
    for (a in 1:n_age_groups) {
      cum_ili_log[ s,a ] = 0 ; // reset the counter
    }
  }
  
  // --------------------------------parameter hierarchical architecture
  for (s in 1:n_season) { // season effect
  for (a in 1:n_age_groups) { // age effect
  SIR_ini[s,a,1] = SIR_ini_mu[1] * 1 * 1; // S
  SIR_ini[s,a,2] = SIR_ini_mu[2] * 1 * 1; // I // * 2^(i_season[s]) *2^(i_age[a])
  SIR_ini[s,a,3] = SIR_ini_mu[3] * 1 * 1; // R // 2^(r_season[s]) *2^(r_age[a])
  prop_ili[s,a]  = prop_ili_mu   * 1 * 2^(prop_ili_age[a]); // 2^(prop_ili_season[s]) *2^(prop_ili_age[a])
  ar[s,a] = 0 ; // set attack rate to zero
  }
  }
  
  // --------------------------------SIR model architecture
  {
    // Declare local variables (only used within these curly brackets and then forgotten, cannot be constrained)
    // parameters that are altered in scenarios
    real beta_j = beta;
    matrix[n_day_fit, n_age_groups] delta_vax_j = delta_vax;
    // compartments
    matrix[n_day_fit,n_age_groups] S_u; // susceptible compartment, unvaccinated
    matrix[n_day_fit,n_age_groups] I_u; // infetious compartment, unvaccinated
    matrix[n_day_fit,n_age_groups] R_u; // recovered compartment, unvaccinated
    matrix[n_day_fit,n_age_groups] S_v; // susceptible compartment, vaccinated
    matrix[n_day_fit,n_age_groups] I_v; // infetious compartment, vaccinated
    matrix[n_day_fit,n_age_groups] R_v; // recovered compartment, vaccinated
    // increments
    real delta_S_u;
    real delta_I_u;
    real delta_R_u;
    real delta_S_v;
    real delta_I_v;
    real delta_R_v;
    real delta_infective_exposures_u;
    real delta_infective_exposures_v;
    array[n_day_fit,n_age_groups] real delta_ili_u; // unvaccinated ili/detectable incidence relative to population size
    array[n_day_fit,n_age_groups] real delta_ili_v; // vaccinated ili/detectable incidence relative to population size
    array[n_day_fit,n_age_groups] real delta_ili_u_abs; // unvaccinated ili/detectable incidence in absolute numbers
    array[n_day_fit,n_age_groups] real delta_ili_v_abs; // vaccinated ili/detectable incidence in absolute numbers
    // previous states
    real prev_S_u[n_age_groups];
    real prev_I_u[n_age_groups];
    real prev_R_u[n_age_groups];
    real prev_S_v[n_age_groups];
    real prev_I_v[n_age_groups];
    real prev_R_v[n_age_groups];
    real prev_delta_ili_u[n_age_groups];
    real prev_delta_ili_v[n_age_groups];
    real prev_delta_ili_u_abs[n_age_groups];
    real prev_delta_ili_v_abs[n_age_groups];
    // current states
    real curr_S_u[n_age_groups];
    real curr_I_u[n_age_groups];
    real curr_R_u[n_age_groups];
    real curr_S_v[n_age_groups];
    real curr_I_v[n_age_groups];
    real curr_R_v[n_age_groups];
    real curr_delta_ili_u[n_age_groups];
    real curr_delta_ili_v[n_age_groups];
    real curr_delta_ili_u_abs[n_age_groups];
    real curr_delta_ili_v_abs[n_age_groups];
    // counters and indicators
    int curr_day; // the current day in which we are looping 
    int n_multi_day = n_multi_day_fit;
    array[n_day_fit*n_daily_time_steps] int daily_counter = daily_counter_fit;
    array[n_day_fit*n_daily_time_steps] int daily_daystart= daily_daystart_fit;
    array[n_day_fit] int season_start = season_start_fit;
    array[n_day_fit] int season_id_day = season_id_day_fit;
    // end: local variables
    
    //////// START-CORE DISEASE DYNAMIC PROCESS (in model)
    // loop through all timesteps
    for (t in 1:n_multi_day){
      
      curr_day=daily_counter[t];
      
      // If it is first day of season AND first iteration of the day
      if ( season_start[curr_day]==1 && daily_daystart[t]==1 ){
        // initiate the compartments 
        for(a in 1:n_age_groups){
          // set the current states
          curr_S_u[a] = SIR_ini[ season_id_day[curr_day], a, 1] ; 
          curr_I_u[a] = SIR_ini[ season_id_day[curr_day], a, 2] ; 
          curr_R_u[a] = SIR_ini[ season_id_day[curr_day], a, 3] ; 
          curr_S_v[a] = 0;  // at start of season, no one is vaccinated
          curr_I_v[a] = 0;  // at start of season, no one is vaccinated
          curr_R_v[a] = 0;  // at start of season, no one is vaccinated
          curr_delta_ili_u[a] = 0; //
          curr_delta_ili_v[a] = 0; //
          curr_delta_ili_u_abs[a] = 0; //
          curr_delta_ili_v_abs[a] = 0; //
          // save the current states
          S_u[curr_day,a] = curr_S_u[a];
          I_u[curr_day,a] = curr_I_u[a];
          R_u[curr_day,a] = curr_R_u[a];
          S_v[curr_day,a] = curr_S_v[a];
          I_v[curr_day,a] = curr_I_v[a];
          R_v[curr_day,a] = curr_R_v[a];
          delta_ili_u[curr_day,a] = curr_delta_ili_u[a];
          delta_ili_v[curr_day,a] = curr_delta_ili_v[a];
          delta_ili_u_abs[curr_day,a] = curr_delta_ili_u_abs[a];
          delta_ili_v_abs[curr_day,a] = curr_delta_ili_v_abs[a];
        } // through age groups
        
      } else { // if not first day of the season, perform the normal step-wise processes
      for(a in 1:n_age_groups){
        // new infections
        delta_infective_exposures_u = dt * beta_j*a_factor[a] * prev_S_u[a] * sum(to_vector(contact_matrix[a, : ]) .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) );
        delta_infective_exposures_v = dt * beta_j*a_factor[a] * prev_S_v[a] * sum(to_vector(contact_matrix[a, : ]) .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) ) * (1 - ve_inf);
        // compute changes due to infections and curing
        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        delta_I_u = delta_infective_exposures_u - prev_I_u[a] * rate_infectious*dt;
        delta_I_v = delta_infective_exposures_v - prev_I_v[a] * rate_infectious*dt;
        delta_R_u = prev_I_u[a]*rate_infectious*dt;
        delta_R_v = prev_I_v[a]*rate_infectious*dt;
        // apply infection/curing changes
        curr_S_u[a] = prev_S_u[a] + delta_S_u;
        curr_S_v[a] = prev_S_v[a] + delta_S_v;
        curr_I_u[a] = prev_I_u[a] + delta_I_u;
        curr_I_v[a] = prev_I_v[a] + delta_I_v;
        curr_R_u[a] = prev_R_u[a] + delta_R_u;
        curr_R_v[a] = prev_R_v[a] + delta_R_v;
        // apply vaccination
        curr_S_u[a] = curr_S_u[a] - (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_S_u[a];
        curr_S_v[a] = curr_S_v[a] + (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_S_u[a];
        curr_R_u[a] = curr_R_u[a] - (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_R_u[a];
        curr_R_v[a] = curr_R_v[a] + (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_R_u[a];
        // collect the delta values over the steps within 1 day
        curr_delta_ili_u[a]     = curr_delta_ili_u[a]     + delta_infective_exposures_u * 1 * prop_ili[ season_id_day[daily_counter_fit[t]], a ];
        curr_delta_ili_v[a]     = curr_delta_ili_v[a]     + delta_infective_exposures_v * (1-ve_ili_cond_inf) * prop_ili[ season_id_day[daily_counter_fit[t]], a ];
        curr_delta_ili_u_abs[a] = curr_delta_ili_u_abs[a] + curr_delta_ili_u[a] * pop_age_group[a,1];
        curr_delta_ili_v_abs[a] = curr_delta_ili_v_abs[a] + curr_delta_ili_v[a] * pop_age_group[a,1];
        // save the current states (once per day)
        if ( daily_daystart_fit[t]==1 ) { 
          S_u[curr_day,a] = curr_S_u[a];
          S_v[curr_day,a] = curr_S_v[a];
          I_u[curr_day,a] = curr_I_u[a];
          I_v[curr_day,a] = curr_I_v[a];
          R_u[curr_day,a] = curr_R_u[a];
          R_v[curr_day,a] = curr_R_v[a];
          delta_ili_u[curr_day,a] = curr_delta_ili_u[a];
          delta_ili_v[curr_day,a] = curr_delta_ili_v[a];
          delta_ili_u_abs[curr_day,a] = curr_delta_ili_u_abs[a];
          delta_ili_v_abs[curr_day,a] = curr_delta_ili_v_abs[a];
          // reset delta values
          curr_delta_ili_u[a] = curr_delta_ili_u[a]*0;
          curr_delta_ili_v[a] = curr_delta_ili_v[a]*0;
          curr_delta_ili_u_abs[a] = curr_delta_ili_u_abs[a]*0;
          curr_delta_ili_v_abs[a] = curr_delta_ili_v_abs[a]*0;
          // backfill the first position of each season
          if ( season_start[curr_day] == 2 ) {
            delta_ili_u[curr_day-1,a] = delta_ili_u[curr_day,a];
            delta_ili_v[curr_day-1,a] = delta_ili_v[curr_day,a];
            delta_ili_u_abs[curr_day-1,a] = delta_ili_u_abs[curr_day,a];
            delta_ili_v_abs[curr_day-1,a] = delta_ili_v_abs[curr_day,a];
          } 
        }
      } // through age groups
      }
      // Update previous values with current values for next iteration
      prev_S_u = curr_S_u;
      prev_I_u = curr_I_u;
      prev_R_u = curr_R_u;
      prev_S_v = curr_S_v;
      prev_I_v = curr_I_v;
      prev_R_v = curr_R_v;
      prev_delta_ili_u = curr_delta_ili_u;
      prev_delta_ili_v = curr_delta_ili_v;
      prev_delta_ili_u_abs = curr_delta_ili_u_abs;
      prev_delta_ili_v_abs = curr_delta_ili_v_abs;
    } // end of multi-daily loop
    //////// END-CORE DISEASE DYNAMIC PROCESS
    
    // convert daily to weekly
    for (t in 1:n_week_fit) {
      for (a in 1:n_age_groups) {
        // define 2 local variables
        int day_start = (t-1)*7+1;
        int day_end = day_start+6;
        delta_ili_incs_weekly[t,a] = sum( delta_ili_u_abs[day_start:day_end,a] ) + sum( delta_ili_v_abs[day_start:day_end,a] );
        delta_ili_percap_weekly[t,a] = delta_ili_incs_weekly[t,a]/pop_age_group[a,1]*100000;
        // save summary stats
        if (ili_obs_notna[t,a]==1) cum_ili_log[ season_id_week_fit[t], a ] = cum_ili_log[ season_id_week_fit[t], a ] + delta_ili_incs_weekly[t,a] ;
      }
    }
    
  } // environment for local variables
  
  
  
  
  // convert to log scale
  for (s in 1:n_season_cum_fit ) {
    for (a in 1:n_age_groups) {
      cum_ili_log[s,a] = log(cum_ili_log[s,a]) ;
    }
  }
  
      
}

model {
  // --------------------------------likelihood part
  
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      if (ili_obs_notna[t,a]==1) target += weight_obs_epi[t,a]*neg_binomial_2_lpmf( ili_obs_incs[t,a] | delta_ili_incs_weekly[t,a]+1e-9, phi ) ; // TODO: remove this line and see if you get priors back
    }
  }
  
  for (s in 1:n_season_cum_fit ) {
    for (a in 1:n_age_groups) {
      target += weight_cum_ili*normal_lpdf( cum_ili_obs_incs_age_log[s,a] | cum_ili_log[s,a] , sigma_cum_ili ) ;
    }
  }
  
  // --------------------------------prior part
  // mock_para ~ normal(0,1);
  // target += normal_lpdf( i_season        | 0 , 0.0001 ) ;
  // target += normal_lpdf( r_season        | 0 , 0.0001 ) ;
  // target += normal_lpdf( prop_ili_season | 0 , 0.0001 ) ;
  
  // prior based on initial fit with flat priors (all countries by AT&IT are good fits)
  // target += normal_lpdf( log(prop_ili_mu) | -1.5 , prior_sigma_prop_ili );// check in R: rnorm(2000,logit(0.1), 3) %>% inv_logit() %>% dens()
  
  // prop_ili_season ~ normal( 0 , sigma_prop_ili_season);
  // prop_ili_age ~    normal( 0 , sigma_prop_ili_age);
  
  // I_ini determined the season timing and certainly be a very low value
  // logit(SIR_ini_mu[2]) ~ normal( logit(0.000003) , prior_sigma_i ); // check in R: rnorm(2000,logit(0.000002),0.4) %>% inv_logit() %>% dens()
  // logit(SIR_ini_mu[1]) ~ normal( logit(0.77) , prior_sigma_s ); // check in R: rnorm(2000,logit(0.85),0.2) %>% inv_logit() %>% dens()
  
  // logit(SIR_ini_mu[1,3]) ~ normal( logit(0.15) , 0.2 ); // check in R: rnorm(2000,logit(0.0015),0.4) %>% inv_logit() %>% dens()
  
  // logit(reciprocal_phi) ~ normal( logit(0.05) , 0.1 ); // check in R: rnorm(2000,logit(0.99),0.1) %>% inv_logit() %>% dens()
  
  // priors on dispersion parameters
  // sigma_prop_ili_age   ~ exponential(0.5); // parameter is the exponential RATE, so BIG numbers mean LOW mean
  // sigma_prop_ili_season~ exponential(0.5); // parameter is the exponential RATE, so BIG numbers mean LOW mean
  // sigma_i              ~ exponential(0.5); // parameter is the exponential RATE, so BIG numbers mean low mean
  
}

generated quantities {
  // --------------------------------declare generated variables
  // we give generated quantities the prefix "gen_"
  array[n_week_fit] real<lower=0> delta_ili_incs_weekly_sum;
  array[n_week_fit] real<lower=0> delta_ili_percap_weekly_sum;
  array[n_week_fit, n_age_groups] real<lower=0> gen_ili_obs_fit;
  array[n_week_fit, n_age_groups] real<lower=0> gen_ili_percap_obs_fit;
  array[n_week_fit] real<lower=0> gen_ili_obs_fit_sum;
  array[n_week_fit] real<lower=0> gen_ili_obs_percap_fit_sum;
  // note: stan does not have 3-dimensional matrices, thus opting for arrays or 2-dimensional matrixes
  // note: matrix[n,m] M[o] creates an array of length o, each element contraining an nxm matrix, M CONFUSINGLY has then dimension [o,n,m]
  array[n_scenario, n_week_proj, n_age_groups ] real<lower=0>  gen_ili_u_obs_proj; // unvaccinated
  array[n_scenario, n_week_proj, n_age_groups ] real<lower=0>  gen_ili_v_obs_proj; // vaccinated
  array[n_scenario, n_week_proj, n_age_groups ] real<lower=0>  gen_ili_t_obs_proj; // total
  array[n_scenario, n_week_proj, n_age_groups ] real<lower=0> gen_ili_u_percap_obs_proj; // per 100 000 of total
  array[n_scenario, n_week_proj, n_age_groups ] real<lower=0> gen_ili_v_percap_obs_proj; // per 100 000 of total
  array[n_scenario, n_week_proj, n_age_groups ] real<lower=0> gen_ili_t_percap_obs_proj; // per 100 000 of total
  array[n_scenario, n_week_proj  ] real<lower=0> gen_ili_u_obs_proj_sum;
  array[n_scenario, n_week_proj  ] real<lower=0> gen_ili_v_obs_proj_sum;
  array[n_scenario, n_week_proj  ] real<lower=0> gen_ili_t_obs_proj_sum;
  array[n_scenario, n_week_proj  ] real<lower=0> gen_ili_u_percap_obs_proj_sum; // per 100 000 of focus group 
  array[n_scenario, n_week_proj  ] real<lower=0> gen_ili_v_percap_obs_proj_sum; // per 100 000 of focus group
  array[n_scenario, n_week_proj  ] real<lower=0> gen_ili_t_percap_obs_proj_sum; // per 100 000 of focus group
  //
  array[n_scenario, n_week_proj, n_age_groups] real gen_delta_ili_u_abs_weekly; // unvaccinated
  array[n_scenario, n_week_proj, n_age_groups] real gen_delta_ili_v_abs_weekly; // vaccinated
  //
  array[n_scenario,n_week_proj,n_age_groups] real pop_a_u;
  array[n_scenario,n_week_proj,n_age_groups] real pop_a_v;
  array[n_scenario,n_week_proj] real pop_u;
  array[n_scenario,n_week_proj] real pop_v;
  //
  real Rnull_eff[n_season];
  real beta_noise;
  
  beta_noise = normal_rng( 1, 0.25 ); // to be applied as factor to log2( beta ), sd=1 interpreted as halfing / doubling of beta
  
  // --------------------------------simulate some quantities of interest
  // computation of Rnull_eff that is not quite correct due to age-structure
  for (season_i in 1:n_season) {
    Rnull_eff[season_i] = Rnull*(1-sum(SIR_ini[ season_i,,3 ]));
  }
  
  // --------------------------------simulate fitted observations
  // simulate fitted observations
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      gen_ili_obs_fit[t,a] = neg_binomial_2_rng( (delta_ili_incs_weekly[t,a]+1e-9) , phi );
      gen_ili_obs_fit[t,a] = fmin( gen_ili_obs_fit[t,a],pop_age_group[a,1] ); // ensure incidence does not exceed compartment
      gen_ili_percap_obs_fit[t,a] = gen_ili_obs_fit[t,a]/pop_age_group[a,1]*100000;
    }
    gen_ili_obs_fit_sum[t] = sum( gen_ili_obs_fit[t,] );
    gen_ili_obs_percap_fit_sum[t] =  gen_ili_obs_fit_sum[t]/pop*100000;
    delta_ili_incs_weekly_sum[t] = sum( delta_ili_incs_weekly[t, ] );
    delta_ili_percap_weekly_sum[t] = delta_ili_incs_weekly_sum[t]/pop*100000;
  }
  
  // --------------------------------simulate projected observations
  for (j in 1:n_scenario) {
    
    // settings for the scenarios
    // define 2 local variables
    real beta_j; // scenario-specific beta
    matrix[n_day_proj, n_age_groups] delta_vax_j; // scenario-specific vaccine uptake
    if ( axis_transmission[j]==0 ) beta_j = 1.0*beta; // more var: 1.0*2^( log2(beta)*beta_noise ) // status quo transmission
    if ( axis_transmission[j]==1 ) beta_j = 0.9*beta; // more var: 0.9*2^( log2(beta)*beta_noise ) // Ooptimistic transmission
    if ( axis_transmission[j]==2 ) beta_j = 1.1*beta; // more var: 1.1*2^( log2(beta)*beta_noise ) // optimisitc transmission
    if ( axis_vax[j]==0 ) delta_vax_j = delta_vax_real; // status quo vaccination
    if ( axis_vax[j]==1 ) delta_vax_j = delta_vax_opti; // optimistic vaccination
    if ( axis_vax[j]==2 ) delta_vax_j = delta_vax_pess; // pessimistic vaccination
    if ( axis_vax[j]==3 ) delta_vax_j = delta_vax_null; // zero vaccination
    // delta_vax_null < delta_vax_pess < delta_vax_real < delta_vax_opti
    // Declare local variables 
    // compartments
    matrix[n_day_proj,n_age_groups] S_u; // susceptible compartment, unvaccinated
    matrix[n_day_proj,n_age_groups] I_u; // infetious compartment, unvaccinated
    matrix[n_day_proj,n_age_groups] R_u; // recovered compartment, unvaccinated
    matrix[n_day_proj,n_age_groups] S_v; // susceptible compartment, vaccinated
    matrix[n_day_proj,n_age_groups] I_v; // infetious compartment, vaccinated
    matrix[n_day_proj,n_age_groups] R_v; // recovered compartment, vaccinated
    // increments
    real delta_S_u;
    real delta_I_u;
    real delta_R_u;
    real delta_S_v;
    real delta_I_v;
    real delta_R_v;
    real delta_infective_exposures_u; 
    real delta_infective_exposures_v;
    array[n_day_fit,n_age_groups] real delta_ili_u; // unvaccinated ili/detectable incidence relative to population size
    array[n_day_fit,n_age_groups] real delta_ili_v; // vaccinated ili/detectable incidence relative to population size
    array[n_day_fit,n_age_groups] real delta_ili_u_abs; // unvaccinated ili/detectable incidence in absolute numbers
    array[n_day_fit,n_age_groups] real delta_ili_v_abs; // vaccinated ili/detectable incidence in absolute numbers
    // previous states
    real prev_S_u[n_age_groups];
    real prev_I_u[n_age_groups];
    real prev_R_u[n_age_groups];
    real prev_S_v[n_age_groups];
    real prev_I_v[n_age_groups];
    real prev_R_v[n_age_groups];
    real prev_delta_ili_u[n_age_groups];
    real prev_delta_ili_v[n_age_groups];
    real prev_delta_ili_u_abs[n_age_groups];
    real prev_delta_ili_v_abs[n_age_groups];
    // current states
    real curr_S_u[n_age_groups];
    real curr_I_u[n_age_groups];
    real curr_R_u[n_age_groups];
    real curr_S_v[n_age_groups];
    real curr_I_v[n_age_groups];
    real curr_R_v[n_age_groups];
    real curr_delta_ili_u[n_age_groups];
    real curr_delta_ili_v[n_age_groups];
    real curr_delta_ili_u_abs[n_age_groups];
    real curr_delta_ili_v_abs[n_age_groups];
    // counters and indicators
    int curr_day; // the current day in which we are looping 
    int n_multi_day = n_multi_day_proj;
    array[n_day_proj*n_daily_time_steps] int daily_counter = daily_counter_proj;
    array[n_day_proj*n_daily_time_steps] int daily_daystart= daily_daystart_proj;
    array[n_day_proj] int season_start = season_start_proj;
    array[n_day_proj] int season_id_day = season_id_day_proj;
    
    //////// START-CORE DISEASE DYNAMIC PROCESS (in generated quantities)
    // loop through all timesteps
    for (t in 1:n_multi_day){
      
      curr_day=daily_counter[t];
      
      // If it is first day of season AND first iteration of the day
      if ( season_start[curr_day]==1 && daily_daystart[t]==1 ){
        // initiate the compartments 
        for(a in 1:n_age_groups){
          // set the current states
          curr_S_u[a] = SIR_ini[ season_id_day[curr_day], a, 1] ; 
          curr_I_u[a] = SIR_ini[ season_id_day[curr_day], a, 2] ; 
          curr_R_u[a] = SIR_ini[ season_id_day[curr_day], a, 3] ; 
          curr_S_v[a] = 0;  // at start of season, no one is vaccinated
          curr_I_v[a] = 0;  // at start of season, no one is vaccinated
          curr_R_v[a] = 0;  // at start of season, no one is vaccinated
          curr_delta_ili_u[a] = 0; //
          curr_delta_ili_v[a] = 0; //
          curr_delta_ili_u_abs[a] = 0; //
          curr_delta_ili_v_abs[a] = 0; //
          // save the current states
          S_u[curr_day,a] = curr_S_u[a];
          I_u[curr_day,a] = curr_I_u[a];
          R_u[curr_day,a] = curr_R_u[a];
          S_v[curr_day,a] = curr_S_v[a];
          I_v[curr_day,a] = curr_I_v[a];
          R_v[curr_day,a] = curr_R_v[a];
          delta_ili_u[curr_day,a] = curr_delta_ili_u[a];
          delta_ili_v[curr_day,a] = curr_delta_ili_v[a];
          delta_ili_u_abs[curr_day,a] = curr_delta_ili_u_abs[a];
          delta_ili_v_abs[curr_day,a] = curr_delta_ili_v_abs[a];
        } // through age groups
        
      } else { // if not first day of the season, perform the normal step-wise processes
      for(a in 1:n_age_groups){
        // new infections
        delta_infective_exposures_u = dt * beta_j*a_factor[a] * prev_S_u[a] * sum(to_vector(contact_matrix[a, : ]) .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) );
        delta_infective_exposures_v = dt * beta_j*a_factor[a] * prev_S_v[a] * sum(to_vector(contact_matrix[a, : ]) .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) ) * (1 - ve_inf);
        // compute changes due to infections and curing
        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        delta_I_u = delta_infective_exposures_u - prev_I_u[a] * rate_infectious*dt;
        delta_I_v = delta_infective_exposures_v - prev_I_v[a] * rate_infectious*dt;
        delta_R_u = prev_I_u[a]*rate_infectious*dt;
        delta_R_v = prev_I_v[a]*rate_infectious*dt;
        // apply infection/curing changes
        curr_S_u[a] = prev_S_u[a] + delta_S_u;
        curr_S_v[a] = prev_S_v[a] + delta_S_v;
        curr_I_u[a] = prev_I_u[a] + delta_I_u;
        curr_I_v[a] = prev_I_v[a] + delta_I_v;
        curr_R_u[a] = prev_R_u[a] + delta_R_u;
        curr_R_v[a] = prev_R_v[a] + delta_R_v;
        // apply vaccination
        curr_S_u[a] = curr_S_u[a] - (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_S_u[a];
        curr_S_v[a] = curr_S_v[a] + (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_S_u[a];
        curr_R_u[a] = curr_R_u[a] - (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_R_u[a];
        curr_R_v[a] = curr_R_v[a] + (delta_vax_j[daily_counter_fit[t-1],a]/n_daily_time_steps) * curr_R_u[a];
        // collect the delta values over the steps within 1 day
        curr_delta_ili_u[a]     = curr_delta_ili_u[a]     + delta_infective_exposures_u * 1 * prop_ili[ season_id_day[daily_counter_fit[t]], a ];
        curr_delta_ili_v[a]     = curr_delta_ili_v[a]     + delta_infective_exposures_v * (1-ve_ili_cond_inf) * prop_ili[ season_id_day[daily_counter_fit[t]], a ];
        curr_delta_ili_u_abs[a] = curr_delta_ili_u_abs[a] + curr_delta_ili_u[a] * pop_age_group[a,1];
        curr_delta_ili_v_abs[a] = curr_delta_ili_v_abs[a] + curr_delta_ili_v[a] * pop_age_group[a,1];
        // save the current states (once per day)
        if ( daily_daystart_fit[t]==1 ) { 
          S_u[curr_day,a] = curr_S_u[a];
          S_v[curr_day,a] = curr_S_v[a];
          I_u[curr_day,a] = curr_I_u[a];
          I_v[curr_day,a] = curr_I_v[a];
          R_u[curr_day,a] = curr_R_u[a];
          R_v[curr_day,a] = curr_R_v[a];
          delta_ili_u[curr_day,a] = curr_delta_ili_u[a];
          delta_ili_v[curr_day,a] = curr_delta_ili_v[a];
          delta_ili_u_abs[curr_day,a] = curr_delta_ili_u_abs[a];
          delta_ili_v_abs[curr_day,a] = curr_delta_ili_v_abs[a];
          // reset delta values
          curr_delta_ili_u[a] = curr_delta_ili_u[a]*0;
          curr_delta_ili_v[a] = curr_delta_ili_v[a]*0;
          curr_delta_ili_u_abs[a] = curr_delta_ili_u_abs[a]*0;
          curr_delta_ili_v_abs[a] = curr_delta_ili_v_abs[a]*0;
          // backfill the first position of each season
          if ( season_start[curr_day] == 2 ) {
            delta_ili_u[curr_day-1,a] = delta_ili_u[curr_day,a];
            delta_ili_v[curr_day-1,a] = delta_ili_v[curr_day,a];
            delta_ili_u_abs[curr_day-1,a] = delta_ili_u_abs[curr_day,a];
            delta_ili_v_abs[curr_day-1,a] = delta_ili_v_abs[curr_day,a];
          } 
        }
      } // through age groups
      }
      // Update previous values with current values for next iteration
      prev_S_u = curr_S_u;
      prev_I_u = curr_I_u;
      prev_R_u = curr_R_u;
      prev_S_v = curr_S_v;
      prev_I_v = curr_I_v;
      prev_R_v = curr_R_v;
      prev_delta_ili_u = curr_delta_ili_u;
      prev_delta_ili_v = curr_delta_ili_v;
      prev_delta_ili_u_abs = curr_delta_ili_u_abs;
      prev_delta_ili_v_abs = curr_delta_ili_v_abs;
    } // end of multi-daily loop
    //////// END-CORE DISEASE DYNAMIC PROCESS
    
    // convert projections daily into weekly
    for (i in 1:n_week_proj) {
      for (a in 1:n_age_groups) {
        // define 2 local variables
        int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
        int day_end = day_start+6;
        gen_delta_ili_u_abs_weekly[j,i,a] = sum( delta_ili_u_abs[day_start:day_end,a] );
        gen_delta_ili_v_abs_weekly[j,i,a] = sum( delta_ili_v_abs[day_start:day_end,a] );
        pop_a_u[j,i,a] = ( S_u[day_start,a]+I_u[day_start,a]+R_u[day_start,a] )*pop_age_group[a,1]+1e-9; // population size of the unvaccinated in age group a
        pop_a_v[j,i,a] = ( S_v[day_start,a]+I_v[day_start,a]+R_v[day_start,a] )*pop_age_group[a,1]+1e-9;// population size of the vaccinated in age group a
      }
      pop_u[j,i] = sum( pop_a_u[j,i, ] )+1e-9; // population size of the unvaccinated 
      pop_v[j,i] = sum( pop_a_v[j,i, ] )+1e-9; // popluation size of the vaccinated
    }
  }
  
  // simulate projections
  for (j in 1:n_scenario) {
    //for (j in 1:1) {
      for (t in 1:n_week_proj) {
        for (a in 1:n_age_groups) {
          gen_ili_u_obs_proj[j,t,a] = neg_binomial_2_rng( gen_delta_ili_u_abs_weekly[j,t,a]+1e-9 , phi ); // add small value to location parameter to avoid it being zero
          gen_ili_u_obs_proj[j,t,a] = fmin( gen_ili_u_obs_proj[j,t,a], pop_a_u[j,t,a] ); // ensure incidence does not exceed compartment
          gen_ili_v_obs_proj[j,t,a] = neg_binomial_2_rng( gen_delta_ili_v_abs_weekly[j,t,a]+1e-9 , phi ); // add small value to location parameter to avoid it being zero
          gen_ili_v_obs_proj[j,t,a] = fmin( gen_ili_v_obs_proj[j,t,a], pop_a_v[j,t,a] ); // ensure incidence does not exceed compartment
          gen_ili_t_obs_proj[j,t,a] = gen_ili_u_obs_proj[j,t,a] + gen_ili_v_obs_proj[j,t,a];
          // rate per 100,000 (as indicated with "percap")
          gen_ili_u_percap_obs_proj[j,t,a] = gen_ili_u_obs_proj[j,t,a]/    pop_a_u[j,t,a]*100000;
          gen_ili_v_percap_obs_proj[j,t,a] = gen_ili_v_obs_proj[j,t,a]/    pop_a_u[j,t,a]*100000;
          gen_ili_t_percap_obs_proj[j,t,a] = gen_ili_t_obs_proj[j,t,a]/pop_age_group[a,1]*100000;
        }
        // sums across age-groups
        gen_ili_u_obs_proj_sum[j,t]=  sum(gen_ili_u_obs_proj[j,t, ]);
        gen_ili_v_obs_proj_sum[j,t]=  sum(gen_ili_v_obs_proj[j,t, ]);
        gen_ili_t_obs_proj_sum[j,t]=  sum(gen_ili_t_obs_proj[j,t, ]);
        //  sums across age-groups, rate per 100,000 (as indicated with "percap")
        gen_ili_u_percap_obs_proj_sum[j,t]=  sum(gen_ili_u_obs_proj[j,t, ])/pop_u[j,t]*100000;
        gen_ili_v_percap_obs_proj_sum[j,t]=  sum(gen_ili_v_obs_proj[j,t, ])/pop_v[j,t]*100000;
        gen_ili_t_percap_obs_proj_sum[j,t]=  sum(gen_ili_t_obs_proj[j,t, ])/       pop*100000;
      }
  }
  
}

