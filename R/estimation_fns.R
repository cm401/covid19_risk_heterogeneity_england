



synth_pop_poisson_reg_new <- function( # data specific inputs
  strata            = c('sex','RGN21NM','IMD_decile','ethnicity_simple','age_band_max','vaccine_v2'),
  case_definition   = 'isoweek2_death', 
  min_isoweek       = 17,       # don't go before May 1st
  max_isoweek       = 112,
  fixed_start_week  = 0,      # don't censor anything prior to this, this effectively becomes start of analysis
  variant_period    = "",
  tier_period       = "",
  IMD_region_data   = IMD_region_data,
  npi_tiers         = npi_tiers,
  region_id         = 'RGN21NM',
  con               = con,   # explicitly pass duckdb connection in
  # model specific inputs
  interaction_pairs = list( c(case_definition, 'RGN21NM') ),
  theta             = NULL,  # this is only if we fit a neg binomial model
  return_prediction = FALSE,
  variant_parameteric = FALSE,
  restriction_parameteric = FALSE,
  lambda            = NULL, # lambda = 0 is no regularization, NULL is default 
  db_table          = 'synth_pop_iso',
  keep_cross_validation_predictions = FALSE
){
  # Create data set
  agg_data_3 <- synth_pop_create_data(strata            = strata,
                                      case_definition   = case_definition, 
                                      min_isoweek       = min_isoweek,       # don't go before May 1st
                                      max_isoweek       = max_isoweek,
                                      fixed_start_week  = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                      variant_period    = variant_period,
                                      tier_period       = tier_period,
                                      IMD_region_data   = IMD_region_data,
                                      npi_tiers         = npi_tiers,
                                      region_id         = region_id,
                                      con               = con,
                                      db_table          = db_table ) 
  
  variable_in <- 'vacc_status' # 'age_v2'  #'IMD_national_quintile' #'ethnicity_simple'  #'vacc_status'
  reference   <- 'not_vaccinated' # 'under_40'  #'IMD1' #'white'  #'not_vaccinated'
  tt <- agg_data_3 %>% group_by(!!sym(variable_in),!!sym(case_definition)) %>% 
    summarise(person_risk_days=sum(person_risk_days),n=sum(n)) %>%
    group_by(!!sym(variable_in)) %>%
    summarise(person_at_risk_days = sum(person_risk_days),n=sum(n))
  reorder  <- tt[[variable_in]]
  reorder  <- c(reference,reorder[-which(reorder==reference)])
  tt <- tt %>% arrange(factor(!!sym(variable_in), levels = reorder))
  table <- epitools::ratetable(tt$n,tt$person_at_risk_days)
  rownames( table ) <- tt[[variable_in]]
  data_irr <- epitools::rateratio(table)
  
  data_irr_tbl      <- data_irr$data %>% as_tibble()
  data_irr_tbl$case_definition <- case_definition
  data_irr_tbl$vacc <- rownames(data_irr$data)
  data_irr_tbl      <- data_irr_tbl %>% mutate(Total     = data_irr_tbl[data_irr_tbl$vacc=='Total',]$Count,
                                               pct_cases = Count/Total * 100) %>%
    filter( ( str_detect(vacc,'unknown') | pct_cases < 0.01 ) | Count < 10 )

  agg_data_3 <- agg_data_3 %>% filter( !( vacc_status %in% data_irr_tbl$vacc ))
  
  if(variant_parameteric){
    agg_data_3 <- agg_data_3 |> 
      left_join(variant %>% mutate(isoweek2 = as.double(isoweek2)),
                by=join_by(!!sym(case_definition)==isoweek2,!!sym(region_id))) |>
      mutate(variant_num = case_when(variant=='WT' ~ 0,
                                     variant=='Alpha' ~ 1,
                                     variant=='Delta' ~ 2,
                                     variant=='Omicron' ~ 3))
  }
  
  if(restriction_parameteric){
    npi_tiers_region     <- get_npi_tiers_region(npi_tiers = npi_tiers, region_id = region_id)
    
    agg_data_3 <- agg_data_3 |> 
      left_join(npi_tiers_region,
                by=join_by(!!sym(case_definition)==isoweek2,!!sym(region_id))) |>
      mutate(restriction_num = case_when(tier_restriction_region=='none' ~ 1,
                                         tier_restriction_region=='level_1' ~ 2,
                                         tier_restriction_region=='level_2' ~ 3,
                                         tier_restriction_region=='level_3' ~ 4,
                                         tier_restriction_region=='national_lockdown'~5,
      ))
  }
  
  if(variant_parameteric){
    agg_data_3                  <- agg_data_3 %>% 
      mutate(offset = case_when(person_risk_days==0 ~ 0,
                                TRUE ~ log(person_risk_days))) %>%
      dplyr::select(c(strata,'vacc_status',!!sym(case_definition),'n',offset,variant_num)) 
    
    h2o_data <- agg_data_3 %>% as.h2o()
    
    covariates <- c(strata,'vacc_status',case_definition,'variant_num')
  } else if(restriction_parameteric) {
    agg_data_3                  <- agg_data_3 %>% 
      mutate(offset = case_when(person_risk_days==0 ~ 0,
                                TRUE ~ log(person_risk_days))) %>%
      dplyr::select(c(strata,'vacc_status',!!sym(case_definition),'n',offset,restriction_num)) 
    
    h2o_data <- agg_data_3 %>% as.h2o()
    
    covariates <- c(strata,'vacc_status',case_definition,'restriction_num')
  } else {
    agg_data_3                  <- agg_data_3 %>% 
      mutate(offset = case_when(person_risk_days==0 ~ 0,
                                TRUE ~ log(person_risk_days))) %>%
      dplyr::select(c(strata,'vacc_status',!!sym(case_definition),'n',offset)) 
    
    h2o_data <- agg_data_3 %>% as.h2o()
    
    covariates <- c(strata,'vacc_status',case_definition)
  }
  
  for(item in c(strata,'vacc_status',case_definition))
  {
    h2o_data[[item]] <- as.factor(h2o_data[[item]])
  }
  
  IMD_options <- c( 'IMD', 'Income','Employment','Education','Health','Crime','Housing','Environment')
  imd_check <- rep(FALSE,length(strata))
  for(comp in IMD_options) imd_check <- imd_check + str_detect(strata,comp)
  
  if( !is.null(lambda) )
    if(lambda==0)
    {
      IMD_var <- strata[as.logical(imd_check)]
      AGE_var <- strata[str_detect(strata,"age")]
      
      if(AGE_var=='age_v2') age_ref <- "50_59" else age_ref <- 64
      
      h2o_data[["vacc_status"]]      <- h2o.relevel(x = h2o_data[["vacc_status"]], y = "not_vaccinated")
      h2o_data[["ethnicity_simple"]] <- h2o.relevel(x = h2o_data[["ethnicity_simple"]], y = "white")
      h2o_data[["sex"]]              <- h2o.relevel(x = h2o_data[["sex"]], y = "F")
      
      if(region_id=='RGN21NM') {
        h2o_data[["RGN21NM"]]          <- h2o.relevel(x = h2o_data[["RGN21NM"]], y = "London")
      } else if(region_id == 'ITL221NM') {
        h2o_data[["ITL221NM"]]          <- h2o.relevel(x = h2o_data[["ITL221NM"]], y = "Inner London - West")
      } else if(region_id == 'ITL321NM') {
        h2o_data[["ITL321NM"]]          <- h2o.relevel(x = h2o_data[["ITL321NM"]], y = 'Kensington & Chelsea and Hammersmith & Fulham')
      }
      
      h2o_data[[IMD_var]]            <- h2o.relevel(x = h2o_data[[IMD_var]], y = paste0(str_split(IMD_var,'_')[[1]][1],"1"))
      h2o_data[[AGE_var]]            <- h2o.relevel(x = h2o_data[[AGE_var]], y = age_ref )
    }
  
  #h2o_data$offset           <- log(h2o_data$n_censored)
  
  if(is.null(theta))         #if we define theta run negative binomial model
  {
    fit <- h2o.glm(x = covariates,
                   y = 'n',
                   interaction_pairs = interaction_pairs,
                   offset_column = 'offset',
                   intercept = TRUE,
                   training_frame = h2o_data,
                   family = "poisson",
                   link = "log",
                   nfolds = 10, 
                   seed = 1234,
                   keep_cross_validation_predictions= keep_cross_validation_predictions,
                   remove_collinear_columns = TRUE,
                   lambda = lambda,  # 0 for no regularization
                   compute_p_values = TRUE
    )
  } else {
    fit <- h2o.glm(x = covariates,
                   y = 'n',
                   interaction_pairs = interaction_pairs,
                   offset_column = 'offset',
                   intercept = TRUE,
                   training_frame = h2o_data,
                   family = "negativebinomial",
                   link = "log",
                   nfolds = 10, 
                   seed = 1234,
                   keep_cross_validation_predictions= keep_cross_validation_predictions,
                   remove_collinear_columns = TRUE,
                   theta = 1/theta,
                   lambda = lambda, # 0 for no regularization
                   compute_p_values = TRUE
    )
  }

  if(return_prediction)
  {
    predict = h2o.predict(fit,h2o_data)
    return(list(fit=fit,prediction=predict,data=agg_data_3))
  } else {
    return(fit)
  }
}

