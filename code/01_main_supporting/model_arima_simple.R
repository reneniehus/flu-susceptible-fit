arima_simple = function( params, data, country_short_input, scenario_tag ){
  
  # create the dataframe for fitting
  data_mock = data %>% 
    filter(country_short == country_short_input, 
           target == params$arima_simple$target, 
           agegroup == params$arima_simple$agegroup) # 
  # obtain fitting data frame
  data_mock_fit = data_mock
  
  # fit flu ARIMA model
  # make projections
  
  #TODO Lydia to continue 
  browser()
  # model = "last_year_burden",
  # country_short = country_short_input,
  # agegroup = params$SIR_simple$agegroup,
  # target = params$SIR_simple$target,
  # scenario_tag = scenario_tag,
  # prediction_type = "sample",
  # sample_or_quantile = sampl_i,
  # value = sample_value,
  # date=date_i+365
  #
  return(df_out)
} 