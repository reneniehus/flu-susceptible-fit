# flu_functions.R
#
# Shared helpers for the data layer. The canonical-table builders (make_data_timeseries_long /
# make_data_timeseries_wide / make_data_season_summary, plus the indicator label+unit helpers) live
# in code/01_main_supporting/gen_model_input.R, and the visual eyeballing() in
# code/01_main_supporting/eyeballing.R. What remains here is the season-time helpers consumed by
# gen_model_input() and the contact-matrix transform.

# ---- |-season helpers (shared, consumed by gen_model_input.R) ----
season_start_year_from_label <- function(season) {
  as.integer(stringr::str_sub(season, 1, 4))
}

add_season_time_columns <- function(df, params=NULL) {
  if (!"date" %in% names(df)) df$date <- as.Date(NA)
  df$date <- as.Date(df$date)
  if (!"season" %in% names(df)) df$season <- NA_character_
  season_start_year <- season_start_year_from_label(df$season)
  season_start_date <- lubridate::ymd(paste0(season_start_year, params$season_start_monthday))
  df %>% mutate(
    season_start_year = season_start_year,
    season_start_date = season_start_date,
    season_day = as.integer(date - season_start_date) + 1L,
    season_week = as.integer(floor((season_day - 1L) / 7L) + 1L),
    iso_week = ifelse(is.na(date), NA_character_, ISOweek::ISOweek(date))
  )
}

# ---- |-contact-matrix transform: 16-group reciprocity-corrected matrices -> 4 model age groups ----
transform_contracts = function(data,params) {

  #stop("Implement the 5th age group!")
  contacts_normalized_all = list()

  if (T){
    xlocations = data$helpers_respicompass$iso2_code
    for (country_i in xlocations$location_name){ # country_i = xlocations$location_name[1]

      # Load original contact matrix
      contacts_orig = data$contact[[country_i]]
      if (length(contacts_orig) == 1 ){ # If contact matrix for this country is not available, skip and go to the next
        next
      }

      # Get age-group sizes
      x_pop = data$demography_respicast$population_pyramid_fine %>%
        filter(country==country_i)
      x_pop_vec = x_pop$population
      x_pop_vec = c(x_pop_vec[1:16], sum(x_pop_vec[17:21]))

      # Add the 80+ age group, assuming it has same per person nr of contacts as 75-79y age group
      x_contacts = cbind(contacts_orig, contacts_orig[,16])
      x_contacts = rbind(x_contacts, x_contacts[16,])
      contacts_orig = x_contacts

      # Fix the contact matrix non-symmetry issue by taking the mean value of the two (taking population size into account)
      # computing m′ij as in [https://cran.r-project.org/web/packages/socialmixr/vignettes/socialmixr.html]
      # mij is the mean number of contacts made by members of age group i with members of age group j
      # thus ROWS are needed for the transmission process -> contacts[i,]
      contacts = NA*contacts_orig
      for (ii in 1:nrow(contacts_orig)){
        for (jj in 1:nrow(contacts_orig)){
          contacts[ii,jj] = (contacts_orig[ii,jj]*x_pop_vec[ii] + contacts_orig[jj,ii]*x_pop_vec[jj]) / (2*x_pop_vec[ii])
        }
      }

      # Get total number of contacts per age group; aka, each element is total number of contacts between age group i and j
      # Specifically, multiply contacts[i,j] by the population size of age group i (which is x_pop_vec[i])
      contacts_total = as.matrix(contacts) * x_pop_vec
      # Need to use the transpose in the pop matrix above such that columns of the population matrix have the same element
      # This is because value contact[i,j] represents number of contacts of person in age group i with persons in age group j
      if (!isSymmetric(contacts_total, check.attributes = FALSE)){
        stop("The contacts_total matrix is not symmetric!")
      }
      ####
      #total_nr_contacts_per_person = sum( contacts_total[row(contacts_total)>=col(contacts_total)] ) / sum(x_pop_vec)
      #contacts = contacts_total / (t(matrix(rep( x_pop_vec, 17), nrow = 17)) * total_nr_contacts_per_person)
      ####

      # Change from 16 age groups to 4
      contacts_total_new = matrix(NA,4,4)
      #
      contacts_total_new[1,1] = contacts_total[1,1]
      #
      tmp_matrix = contacts_total[2:3,2:3]
      contacts_total_new[2,2] = sum(tmp_matrix[row(tmp_matrix)>=col(tmp_matrix)]) # This part ensures that we add diagonal + only one off-diagonal part (e.g., upper but not lower), ensuring we don't count things twice
      contacts_total_new[1,2] = sum(contacts_total[1,2:3])
      contacts_total_new[2,1] = sum(contacts_total[2:3,1])
      #
      tmp_matrix = contacts_total[4:13,4:13]
      contacts_total_new[3,3] = sum(tmp_matrix[row(tmp_matrix)>=col(tmp_matrix)])
      contacts_total_new[1,3] = sum(contacts_total[1,4:13])
      contacts_total_new[3,1] = sum(contacts_total[4:13,1])
      contacts_total_new[2,3] = sum(contacts_total[2:3,4:13])
      contacts_total_new[3,2] = sum(contacts_total[4:13,2:3])
      #
      tmp_matrix = contacts_total[14:17,14:17]
      contacts_total_new[4,4] = sum(tmp_matrix[row(tmp_matrix)>=col(tmp_matrix)])
      contacts_total_new[1,4] = sum(contacts_total[1,14:17])
      contacts_total_new[4,1] = sum(contacts_total[14:17,1])
      contacts_total_new[2,4] = sum(contacts_total[2:3,14:17])
      contacts_total_new[4,2] = sum(contacts_total[14:17,2:3])
      contacts_total_new[3,4] = sum(contacts_total[4:13,14:17])
      contacts_total_new[4,3] = sum(contacts_total[14:17,4:13])

      # Go from total number of contacts back to mij: is the mean number of contacts made by members of age group i with members of age group j
      x_new_pop = data$demography_respicast$population_pyramid %>% filter(country == country_i) %>%
        select(age_group,population) %>% deframe()
      x_new_pop = x_new_pop[c("0-4","5-14","15-64","65+")] # ensure right ordering
      contacts_recovered = contacts_total_new / x_new_pop

      ## obtain matrix such at average number of contacts per person equalt to one
      # Get total mean number of contacts per person
      total_nr_contacts_per_person = sum( contacts_total_new[ row(contacts_total_new)>=col(contacts_total_new) ] ) / sum(x_pop_vec)
      # Get a new contact matrix with only 4 age groups, such that average number of contacts per person equals to one
      x_pop_matrix = t(matrix(rep(x_new_pop,4), nrow=4))
      # The new contact matrix where elements are per person contacts between age group i and j such that the population-weighted average number of contacts is 1
      contacts_normalized = contacts_total_new / (x_pop_matrix * total_nr_contacts_per_person)

      # NOTE (scaling choice -- which matrix to carry forward). We store contacts_recovered: the
      # reciprocity-corrected matrix in mean-contacts-per-person units, NOT rescaled. Left as-is per
      # decision, but for a future revision the more principled choice -- given that this project FIXES
      # R0 from the literature rather than fitting it -- is to carry a matrix normalised by its dominant
      # eigenvalue (spectral radius):
      #   * In age-structured models R0 is the dominant eigenvalue of the next-generation matrix (NGM),
      #     in which the contact matrix enters multiplicatively (Diekmann & Heesterbeek 2000; Diekmann,
      #     Heesterbeek & Roberts 2010, J R Soc Interface 7:873). Dividing the reciprocity-corrected
      #     matrix (Wallinga, Teunis & Kretzschmar 2006, Am J Epidemiol 164:936; socialmixr, Funk 2018)
      #     by its leading eigenvalue gives spectral radius 1, so the transmission rate beta calibrates
      #     DIRECTLY to the target R0, independent of the matrix's absolute scale.
      #   * That absolute scale is the least transferable part of synthetic matrices (Prem, Cook & Jit
      #     2017, PLoS Comput Biol 13:e1005697; Mossong et al. 2008, PLoS Med 5:e74): it varies with the
      #     survey and the reconstruction. An un-normalised matrix ties R0 to that arbitrary scale,
      #     conflating pathogen transmissibility with contact-data scaling; only the RELATIVE age-mixing
      #     structure is estimated robustly, and that is exactly what eigenvalue-normalisation preserves.
      #   * contacts_normalized (computed just above) instead rescales by the population-weighted MEAN
      #     contacts per person -- a convenient proxy, but it equals the spectral-radius normalisation
      #     only when the mean tracks the leading eigenvalue. For fixing R0 the eigenvalue normalisation
      #     is the exact one; do this when the age/vaccination-structured model is revisited.
      contacts_normalized_all[[country_i]] = contacts_recovered
    }
  }

  # create EU average by summing all locations
  contacts_collect = contacts_normalized_all[[1]]*0
  collect_counter = 0
  for (country_i in names(contacts_normalized_all)) {
    contacts_collect = contacts_collect + contacts_normalized_all[[country_i]]
    collect_counter = collect_counter + 1
  }
  EU_contacts = contacts_collect / collect_counter

  # fix assymmetry
  x_new_pop = data$demography_respicast$population_pyramid %>%
    group_by(age_group) %>% summarise(sum=sum(population)) %>% deframe()
  x_new_pop = x_new_pop[c("0-4","5-14","15-64","65+")]
  contacts_orig = EU_contacts
  contacts_new = NA*contacts_orig
  for (ii in 1:nrow(contacts_orig)){
    for (jj in 1:nrow(contacts_orig)){
      contacts_new[ii,jj] = (contacts_orig[ii,jj]*x_new_pop[ii] + contacts_orig[jj,ii]*x_new_pop[jj]) / (2*x_new_pop[ii])
    }
  }

  contacts_normalized_all[["EU"]] = contacts_new
  return(contacts_normalized_all)
}
