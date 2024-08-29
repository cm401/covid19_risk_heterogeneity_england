#### Helper functions to process data


# additional information for dataset --------------------------------------

# Add dominant variant by region
get_variant_data <- function(duckdb_con,
                             db_table = db_table 
)
{
  variant <- tbl(duckdb_con, db_table ) |> 
    group_by(isoweek2,!!sym(region_id)) |> 
    summarise(mean_sgtf = mean(s_negative_E1,na.rm = TRUE)) |> 
    mutate(variant = case_when(isoweek2 < 65 & is.na(mean_sgtf) ~ "WT", #early on we don't have this data
                               isoweek2 < 65 & mean_sgtf < 0.5 ~ "WT",
                               isoweek2 < 90 & mean_sgtf > 0.5 ~ "Alpha",
                               isoweek2 > 65 & mean_sgtf < 0.5 & isoweek2 < 104 ~ "Delta",
                               isoweek2 > 90 & mean_sgtf > 0.5 ~ "Omicron", 
                               isoweek2 > 104                  ~ "Omicron",
                               TRUE                            ~ "Unknown" )) |> 
    collect()
  
  return(variant)
}

# NPI tiers by ltla
get_npi_tiers_ltla <- function(npi_tiers=npi_tiers) 
{
  npi_tiers_ltla <- npi_tiers %>% dplyr::select(ltla_code,isoweek2,tier_restriction) %>%
    rename(tier_restriction_ltla=tier_restriction)
  
  return(npi_tiers_ltla)
}

# NPI tiers by region
get_npi_tiers_region <- function(npi_tiers=npi_tiers) 
{
  npi_tiers_region <- npi_tiers %>% dplyr::select(RGN21NM,isoweek2,tier_restriction) %>%
    group_by(isoweek2,RGN21NM,tier_restriction) %>%
    summarise(count = n()) %>% 
    group_by(isoweek2,RGN21NM) %>%
    slice(which.max(count)) %>%
    filter(!is.na(RGN21NM)) %>%
    dplyr::select(-count) %>%
    rename(tier_restriction_region=tier_restriction)
  
  return(npi_tiers_region)
}  

# create majority tier restriction by region aggregation 
get_region_restrictions <- function(npi_tiers = npi_tiers,
                                    region_id = region_id )
{
  region_restriction <- npi_tiers %>% 
    dplyr::select(ltla_code,isoweek2,tier_restriction,!!sym(region_id)) %>% 
    group_by_at(c('isoweek2',region_id,'tier_restriction')) %>% 
    summarise(n=n()) %>% 
    filter(!is.na(!!sym(region_id))) %>%
    group_by_at(c('isoweek2',region_id)) %>% 
    top_n(1, tier_restriction) %>%       # use majority tier restriction across region for that isoweek
    dplyr::select(-n) %>% ungroup()
  
  return(region_restriction)
}

# process dataset ---------------------------------------------------------

# censor by Pillar 2 PCR-confirmed cases
synth_pop_censor_pillar2pcr <- function(strata             = strata,
                                        fixed_start_week   = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                        region_id          = region_id,
                                        variant            = variant,
                                        region_restriction = region_restriction,
                                        npi_tiers_ltla     = npi_tiers_ltla,
                                        con                = con,
                                        db_table           = db_table 
)
{
  agg_data <- tbl(con, db_table ) |>
    left_join(IMD_region_data, by=join_by(ltla_code==LTLA19CD), copy = TRUE) |>
    filter(sex!='U' & ethnicity_simple != 'unknown' & !is.na(age_20200101) ) |>
    filter(fixed_start_week<isoweek2_pillar2pcr|is.na(isoweek2_pillar2pcr)) %>%
    mutate(backup_vacc = coalesce(vacc_first,vacc_second,vacc_third,vacc_fourth),
           backup_vacc = case_when((backup_vacc == "Pfizer" | backup_vacc == "Moderna") ~ 'mRNA',
                                   (backup_vacc == "AZ" | backup_vacc == "Janssen") ~ 'Adenovirus',
                                   TRUE ~ backup_vacc),
           vacc_first_dose = case_when((vacc_first == "Pfizer" | vacc_first == "Moderna") ~ 'mRNA',
                                       (vacc_first == "AZ" | vacc_first == "Janssen") ~ 'Adenovirus',
                                       !is.na(backup_vacc) ~ backup_vacc,
                                       TRUE ~ 'unknown'),
           vaccine_v2 = case_when((vacc_first == "Pfizer" | vacc_first == "Moderna") & (vacc_second == "Pfizer" | vacc_second == "Moderna") ~ "mRNA",
                                  is.na(vacc_first) & (vacc_second == "Pfizer" | vacc_second == "Moderna") ~ "mRNA",
                                  (vacc_first == "Pfizer" | vacc_first == "Moderna") & is.na(vacc_second) ~ "mRNA",
                                  (vacc_first == "AZ" | vacc_first == "Janssen") & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ "Adenovirus",
                                  (vacc_first == "AZ" | vacc_first == "Janssen") & is.na(vacc_second) ~ "Adenovirus",
                                  is.na(vacc_first) & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ "Adenovirus",
                                  (vacc_first == "AZ" | vacc_first == "Janssen") & (vacc_second == "Pfizer" | vacc_second == "Moderna" | vacc_second == 'Novavax' ) ~ 'Mixed_dose',
                                  (vacc_first == "Pfizer" | vacc_first == "Moderna" | vacc_first == 'Novavax') & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ 'Mixed_dose',
                                  (vacc_first == "Novavax" ) & (vacc_second == "Novavax" ) ~ "Novavax",
                                  !is.na(backup_vacc) ~ backup_vacc,
                                  TRUE ~ 'unknown'),
           vacc_status = case_when(is.na(date_first_iso) ~ 'not_vaccinated',
                                   date_first_iso >= isoweek2_pillar2pcr ~ 'not_vaccinated',
                                   date_first_iso < isoweek2_pillar2pcr & date_first_iso + 3 >= isoweek2_pillar2pcr & (date_second_iso >= isoweek2_pillar2pcr | is.na(date_second_iso) ) ~ paste0('first_dose_<21d_',vacc_first_dose),
                                   date_first_iso + 3 < isoweek2_pillar2pcr & (date_second_iso >= isoweek2_pillar2pcr | is.na(date_second_iso) ) ~ paste0('first_dose_21+d_',vacc_first_dose),
                                   date_second_iso < isoweek2_pillar2pcr & isoweek2_pillar2pcr - ( date_second_iso ) <= 2 ~ paste0('second_dose_<14d_',vaccine_v2),
                                   date_second_iso + 2 < isoweek2_pillar2pcr & isoweek2_pillar2pcr - ( date_second_iso ) < 10 ~ paste0('second_dose_14+d_',vaccine_v2),
                                   date_second_iso + 2 < isoweek2_pillar2pcr & isoweek2_pillar2pcr - ( date_second_iso ) <= 17 ~ paste0('wane_10_18_',vaccine_v2),
                                   date_second_iso + 2 < isoweek2_pillar2pcr & isoweek2_pillar2pcr - ( date_second_iso ) > 17 ~ paste0('wane_over_18_',vaccine_v2),
                                   TRUE ~ 'unknown')) |>
    left_join(variant,by=join_by(isoweek2_pillar2pcr==isoweek2,!!sym(region_id)), copy = TRUE) %>%
    left_join(npi_tiers_ltla,by=join_by(isoweek2_pillar2pcr==isoweek2,ltla_code), copy = TRUE) %>%
    left_join(region_restriction,by=join_by(isoweek2_pillar2pcr==isoweek2,!!sym(region_id)), copy = TRUE) %>%
    mutate(age_simple = case_when(age_20200101 < 19 ~ "Children",
                                  age_20200101 >= 19 & age_20200101 < 65 ~ "Adults",
                                  age_20200101 >= 65 ~ "old_age")) %>% 
    mutate(tier_restriction_ltla = case_when(tier_restriction_ltla %in% c("none", "national_lockdown", "level_3", "level_2", "level_1")~ tier_restriction_ltla,
                                             TRUE ~ 'none'),
           tier_restriction = case_when(tier_restriction %in% c("none", "national_lockdown", "level_3", "level_2", "level_1")~ tier_restriction,
                                        TRUE ~ 'none')) |>
    mutate(age_v2 = case_when(age_20200101 < 40 ~ "under_40",
                              age_20200101 >= 40 & age_20200101 < 50 ~ "40_49",
                              age_20200101 >= 50 & age_20200101 < 60 ~ "50_59",
                              age_20200101 >= 60 & age_20200101 < 70 ~ "60_69",
                              age_20200101 >= 70 & age_20200101 < 80 ~ "70_79",
                              age_20200101 >= 80 ~ "above_80")) %>% 
    filter(!is.na(isoweek2_pillar2pcr)) %>% # no infection event
    group_by_at(c(strata, 'vacc_status','isoweek2_pillar2pcr')) %>%  #always censor by pillar2pcr
    summarise(n=n()) |> 
    arrange(-desc('isoweek2_pillar2pcr'), .by_group = TRUE) |>
    mutate(cum_n = cumsum(n)) |> collect()
  
  return(agg_data)
}

# create event count
synth_pop_create_event_counts <- function(strata             = strata,
                                          case_definition    = case_definition,
                                          fixed_start_week   = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                          region_id          = region_id,
                                          variant            = variant,
                                          region_restriction = region_restriction,
                                          npi_tiers_ltla     = npi_tiers_ltla,
                                          con                = con,
                                          db_table           = db_table 
)
{
  agg_data_case_def <- tbl(con, db_table ) |>
    left_join(IMD_region_data, by=join_by(ltla_code==LTLA19CD), copy = TRUE) |>
    filter(sex!='U' & ethnicity_simple != 'unknown' & !is.na(age_20200101) ) |>
    filter(fixed_start_week<isoweek2_pillar2pcr|is.na(isoweek2_pillar2pcr)) %>%
    mutate(case_type   = case_definition,
           backup_vacc = coalesce(vacc_first,vacc_second,vacc_third,vacc_fourth),
           backup_vacc = case_when((backup_vacc == "Pfizer" | backup_vacc == "Moderna") ~ 'mRNA',
                                   (backup_vacc == "AZ" | backup_vacc == "Janssen") ~ 'Adenovirus',
                                   TRUE ~ backup_vacc),
           vacc_first_dose = case_when((vacc_first == "Pfizer" | vacc_first == "Moderna") ~ 'mRNA',
                                       (vacc_first == "AZ" | vacc_first == "Janssen") ~ 'Adenovirus',
                                       !is.na(backup_vacc) ~ backup_vacc,
                                       TRUE ~ 'unknown'),
           vaccine_v2 = case_when((vacc_first == "Pfizer" | vacc_first == "Moderna") & (vacc_second == "Pfizer" | vacc_second == "Moderna") ~ "mRNA",
                                  is.na(vacc_first) & (vacc_second == "Pfizer" | vacc_second == "Moderna") ~ "mRNA",
                                  (vacc_first == "Pfizer" | vacc_first == "Moderna") & is.na(vacc_second) ~ "mRNA",
                                  (vacc_first == "AZ" | vacc_first == "Janssen") & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ "Adenovirus",
                                  (vacc_first == "AZ" | vacc_first == "Janssen") & is.na(vacc_second) ~ "Adenovirus",
                                  is.na(vacc_first) & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ "Adenovirus",
                                  (vacc_first == "AZ" | vacc_first == "Janssen") & (vacc_second == "Pfizer" | vacc_second == "Moderna" | vacc_second == 'Novavax' ) ~ 'Mixed_dose',
                                  (vacc_first == "Pfizer" | vacc_first == "Moderna" | vacc_first == 'Novavax') & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ 'Mixed_dose',
                                  (vacc_first == "Novavax" ) & (vacc_second == "Novavax" ) ~ "Novavax",
                                  !is.na(backup_vacc) ~ backup_vacc,
                                  TRUE ~ 'unknown'),
           vacc_status = case_when(is.na(date_first_iso) ~ 'not_vaccinated',
                                   date_first_iso >= !!sym(case_definition) ~ 'not_vaccinated',
                                   date_first_iso < !!sym(case_definition) & date_first_iso + 3 >= !!sym(case_definition)  & (date_second_iso >= !!sym(case_definition) | is.na(date_second_iso) ) ~ paste0('first_dose_<21d_',vacc_first_dose),
                                   date_first_iso + 3 < !!sym(case_definition) & (date_second_iso >= !!sym(case_definition) | is.na(date_second_iso) ) ~ paste0('first_dose_21+d_',vacc_first_dose),
                                   date_second_iso < !!sym(case_definition) & !!sym(case_definition) - ( date_second_iso ) <= 2 ~ paste0('second_dose_<14d_',vaccine_v2),
                                   date_second_iso + 2 < !!sym(case_definition) & !!sym(case_definition) - ( date_second_iso ) < 10 ~ paste0('second_dose_14+d_',vaccine_v2),
                                   date_second_iso + 2 < !!sym(case_definition) & !!sym(case_definition) - ( date_second_iso ) <= 17 ~ paste0('wane_10_18_',vaccine_v2),
                                   date_second_iso + 2 < !!sym(case_definition) & !!sym(case_definition) - ( date_second_iso ) > 17 ~ paste0('wane_over_18_',vaccine_v2),
                                   TRUE ~ 'unknown')) |>
    mutate(is_first_episode = case_when(str_detect(case_type,'death') & abs(as.Date(dod) - as.Date(specimen_date_E1))>28 ~ FALSE,  
                                        str_detect(case_type,'hosp') & abs(as.Date(dateadmission_nhse_E1) - as.Date(specimen_date_E1))>28 ~ FALSE,
                                        TRUE ~ TRUE),
           is_pillar2pcr    = case_when( (pillar_E1=="Pillar 2" & str_detect(case_category_E1,"PCR")) ~ TRUE,
                                         TRUE ~ FALSE)) |>       # also ensure that these are pillar2 pcr confirmed cases                                                                                         # for deaths require dod to be within 28 days from specimen date 
    left_join(variant,by=join_by(!!sym(case_definition)==isoweek2,!!sym(region_id)), copy = TRUE) %>%
    left_join(npi_tiers_ltla,by=join_by(!!sym(case_definition)==isoweek2,ltla_code), copy = TRUE) %>%
    left_join(region_restriction,by=join_by(!!sym(case_definition)==isoweek2,!!sym(region_id)), copy = TRUE) %>%
    mutate(age_simple = case_when(age_20200101 < 19 ~ "Children",
                                  age_20200101 >= 19 & age_20200101 < 65 ~ "Adults",
                                  age_20200101 >= 65 ~ "old_age")) %>% 
    mutate(tier_restriction_ltla = case_when(tier_restriction_ltla %in% c("none", "national_lockdown", "level_3", "level_2", "level_1")~ tier_restriction_ltla,
                                             TRUE ~ 'none'),
           tier_restriction = case_when(tier_restriction %in% c("none", "national_lockdown", "level_3", "level_2", "level_1")~ tier_restriction,
                                        TRUE ~ 'none')) |>
    mutate(age_v2 = case_when(age_20200101 < 40 ~ "under_40",
                              age_20200101 >= 40 & age_20200101 < 50 ~ "40_49",
                              age_20200101 >= 50 & age_20200101 < 60 ~ "50_59",
                              age_20200101 >= 60 & age_20200101 < 70 ~ "60_69",
                              age_20200101 >= 70 & age_20200101 < 80 ~ "70_79",
                              age_20200101 >= 80 ~ "above_80")) %>% 
    filter(is_first_episode==TRUE ) %>%         # only deaths within 28days of E1 specimen (regardless of test)
    filter(!is.na(!!sym(case_definition)) ) |>  # no event took place
    group_by_at(c(strata,'vacc_status',case_definition)) %>%  
    summarise(n=n()) |> 
    arrange(-desc(case_definition), .by_group = TRUE) |>
    collect()
  
  return(agg_data_case_def)
}

# create overall dataset for regression
synth_pop_create_data <- function(strata            = strata,
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
{
  variant            <- get_variant_data(con,db_table)
  region_restriction <- get_region_restrictions(npi_tiers = npi_tiers, region_id = region_id)
  npi_tiers_ltla     <- get_npi_tiers_ltla(npi_tiers = npi_tiers)
  
  #get input tables
  input_data_tbl <- synth_pop_generate_input_tables( strata = strata, max_isoweek = 112, IMD_region_data = IMD_region_data, 
                                                     npi_tiers = npi_tiers, region_id = region_id, con = con, db_table = db_table )
  
  #input_data_tbl %>% group_by(isoweek,vacc_status) %>% summarise(n=sum(n_total)) %>% pivot_wider(names_from = vacc_status,values_from = n, values_fill = list(n = 0)) |> View()
  
  if('tier_restriction' %in% colnames(input_data_tbl))
  {
    input_data_tbl <- input_data_tbl %>%
      mutate(tier_restriction = replace_na(tier_restriction,'no_event')) 
  }
  
  # censor by infections only + take deaths associated with first infections only.
  agg_data_2 <- synth_pop_censor_pillar2pcr(strata             = strata,
                                            fixed_start_week   = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                            region_id          = region_id,
                                            variant            = variant,
                                            region_restriction = region_restriction,
                                            npi_tiers_ltla     = npi_tiers_ltla,
                                            con                = con,
                                            db_table           = db_table )
  
  # if we restrict to sub-time periods prior to computation of person-at-risk days
  if( !is.null(variant_period) & variant_period!="")
  {
    agg_data_2 <- agg_data_2 |> 
      left_join(variant,by=join_by(isoweek2_pillar2pcr==isoweek2,!!sym(region_id)), copy = TRUE) |>
      filter(variant==variant_period)
    
    min_isoweek <- max(min_isoweek, min(agg_data_2$isoweek2_pillar2pcr))
    max_isoweek <- min(max(agg_data_2$isoweek2_pillar2pcr),max_isoweek)
  }
  
  if( !is.null(tier_period) & tier_period!="")
  {
    agg_data_2 <- agg_data_2 |> 
      left_join(region_restriction,by=join_by(isoweek2_pillar2pcr==isoweek2,!!sym(region_id)), copy = TRUE)|> 
      filter(tier_restriction==tier_period)
    
    min_isoweek <- min(agg_data_2$isoweek2_pillar2pcr)
    max_isoweek <- max(agg_data_2$isoweek2_pillar2pcr)
  }
  
  if(is.na(min_isoweek)) min_isoweek <- 0
  if(is.na(max_isoweek)) max_isoweek <- 10000
  
  # Create person-days-at-risk
  agg_data_2 <- agg_data_2 %>% full_join(input_data_tbl %>% rename('isoweek2_pillar2pcr'='isoweek'), copy = TRUE ) %>% 
    filter(min_isoweek<=isoweek2_pillar2pcr & max_isoweek>=isoweek2_pillar2pcr) %>%            # CHECK THIS -- I THINK WE NEED TO BE CAREFUL HERE
    group_by_at(c(strata, 'vacc_status')) %>%
    dplyr::arrange(-desc(isoweek2_pillar2pcr), .by_group = TRUE) |>
    fill(cum_n) |>
    mutate(n          = replace_na(n,0),
           n_censored = replace_na(cum_n,0)) |>
    dplyr::mutate(person_risk_days = pmax(0,replace_na(7 * (replace_na(n_total, 0) - n_censored),0))) %>%    #we never allow this to go down
    dplyr::rename(n_pillar2pcr=n) %>%
    filter(!is.na(isoweek2_pillar2pcr) )         # these are all the individuals in the synthetic population with out any cases
  
  agg_data_case_def <- synth_pop_create_event_counts(strata             = strata,
                                                     case_definition    = case_definition,
                                                     fixed_start_week   = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                                     region_id          = region_id,
                                                     variant            = variant,
                                                     region_restriction = region_restriction,
                                                     npi_tiers_ltla     = npi_tiers_ltla,
                                                     con                = con,
                                                     db_table           = db_table )
  
  agg_data_2 <- agg_data_2 %>% 
    left_join(agg_data_case_def,by=c(strata,'vacc_status','isoweek2_pillar2pcr'= case_definition)) %>%
    rename(!!case_definition := 'isoweek2_pillar2pcr') %>%
    mutate(n = replace_na(n,0))    # if we can't match we must have had zero cases in that week
  #filter(!is.na(n))         
  
  agg_data_3 <-  agg_data_2 %>% 
    filter((!!sym(case_definition ))>min_isoweek & (!!sym(case_definition ))<max_isoweek ) %>%
    filter(n_censored>0) 
  
  return(agg_data_3)
}

# create input populations by week [save this for strata / db_tables so that we don't need to regenerate unnecessarily]
synth_pop_generate_input_tables <- function(strata            = c('sex','RGN21NM','IMD_decile','ethnicity_simple','age_band_max','vaccine_v2'),
                                            max_isoweek       = 112,
                                            IMD_region_data   = IMD_region_data,
                                            npi_tiers         = npi_tiers,
                                            region_id         = 'RGN21NM',
                                            con               = con,
                                            FORCE_REGENERATE  = FALSE,
                                            db_table          = db_table )
{
  variant            <- get_variant_data(con,db_table)
  region_restriction <- get_region_restrictions(npi_tiers = npi_tiers, region_id = region_id)
  npi_tiers_ltla     <- get_npi_tiers_ltla(npi_tiers = npi_tiers)
  
  #check if input table exists, otherwise create
  tables           <- dbListTables(con)
  input_table_name <- paste0("input_table_",paste0(strata,collapse = '.'))
  
  if(db_table != "synth_pop_iso")     # synth_pop_iso is already generated, if we use a different database with different variable contents append and regen
    input_table_name <- paste0(input_table_name,'.',db_table)
  
  if( !FORCE_REGENERATE & input_table_name %in% tables )
  {
    input_data_tbl <- tbl(con, input_table_name, check_from = FALSE) |> collect()
  } else {
    input_data_tbl <- tibble(NULL)
    
    for(input_week in 1:max_isoweek)
    {
      input_data_tmp <- tbl(con, db_table ) |> 
        left_join(IMD_region_data,by=join_by(ltla_code==LTLA19CD), copy = TRUE) |>
        filter(sex!='U' & ethnicity_simple != 'unknown' & !is.na(age_20200101) ) |>
        mutate(backup_vacc = coalesce(vacc_first,vacc_second,vacc_third,vacc_fourth),
               backup_vacc = case_when((backup_vacc == "Pfizer" | backup_vacc == "Moderna") ~ 'mRNA',
                                       (backup_vacc == "AZ" | backup_vacc == "Janssen") ~ 'Adenovirus',
                                       TRUE ~ backup_vacc),
               vacc_first_dose = case_when((vacc_first == "Pfizer" | vacc_first == "Moderna") ~ 'mRNA',
                                           (vacc_first == "AZ" | vacc_first == "Janssen") ~ 'Adenovirus',
                                           !is.na(backup_vacc) ~ backup_vacc,
                                           TRUE ~ 'unknown'),
               vaccine_v2 = case_when((vacc_first == "Pfizer" | vacc_first == "Moderna") & (vacc_second == "Pfizer" | vacc_second == "Moderna") ~ "mRNA",
                                      is.na(vacc_first) & (vacc_second == "Pfizer" | vacc_second == "Moderna") ~ "mRNA",
                                      (vacc_first == "Pfizer" | vacc_first == "Moderna") & is.na(vacc_second) ~ "mRNA",
                                      (vacc_first == "AZ" | vacc_first == "Janssen") & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ "Adenovirus",
                                      (vacc_first == "AZ" | vacc_first == "Janssen") & is.na(vacc_second) ~ "Adenovirus",
                                      is.na(vacc_first) & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ "Adenovirus",
                                      (vacc_first == "AZ" | vacc_first == "Janssen") & (vacc_second == "Pfizer" | vacc_second == "Moderna" | vacc_second == 'Novavax' ) ~ 'Mixed_dose',
                                      (vacc_first == "Pfizer" | vacc_first == "Moderna" | vacc_first == 'Novavax') & (vacc_second == "AZ" | vacc_second == "Janssen" ) ~ 'Mixed_dose',
                                      (vacc_first == "Novavax" ) & (vacc_second == "Novavax" ) ~ "Novavax",
                                      !is.na(backup_vacc) ~ backup_vacc,
                                      TRUE ~ 'unknown'),
               vacc_status = case_when(is.na(date_first_iso) ~ 'not_vaccinated',
                                       date_first_iso >= input_week ~ 'not_vaccinated',
                                       date_first_iso < input_week & date_first_iso + 3 >= input_week & (date_second_iso >= input_week | is.na(date_second_iso) ) ~ paste0('first_dose_<21d_',vacc_first_dose),
                                       date_first_iso + 3 < input_week & (date_second_iso >= input_week | is.na(date_second_iso) ) ~ paste0('first_dose_21+d_',vacc_first_dose),
                                       date_second_iso < input_week & input_week - ( date_second_iso ) <= 2 ~ paste0('second_dose_<14d_',vaccine_v2),
                                       date_second_iso + 2 < input_week & input_week - ( date_second_iso ) < 10 ~ paste0('second_dose_14+d_',vaccine_v2),
                                       date_second_iso + 2 < input_week & input_week - ( date_second_iso ) <= 17 ~ paste0('wane_10_18_',vaccine_v2),
                                       date_second_iso + 2 < input_week & input_week - ( date_second_iso ) > 17 ~ paste0('wane_over_18_',vaccine_v2),
                                       TRUE ~ 'unknown')) |>
        left_join(variant,by=join_by(isoweek2_pillar2pcr==isoweek2,!!sym(region_id)), copy = TRUE) %>%
        left_join(npi_tiers_ltla,by=join_by(isoweek2_pillar2pcr==isoweek2,ltla_code), copy = TRUE) %>%
        left_join(region_restriction,by=join_by(isoweek2_pillar2pcr==isoweek2,!!sym(region_id)), copy = TRUE) %>%
        mutate(age_simple = case_when(age_20200101 < 19 ~ "Children",
                                      age_20200101 >= 19 & age_20200101 < 65 ~ "Adults",
                                      age_20200101 >= 65 ~ "old_age")) %>% 
        mutate(tier_restriction_ltla = case_when(tier_restriction_ltla %in% c("none", "national_lockdown", "level_3", "level_2", "level_1")~ tier_restriction_ltla,
                                                 TRUE ~ 'none'),
               tier_restriction = case_when(tier_restriction %in% c("none", "national_lockdown", "level_3", "level_2", "level_1")~ tier_restriction,
                                            TRUE ~ 'none')
        ) |>
        mutate(age_v2 = case_when(age_20200101 < 40 ~ "under_40",
                                  age_20200101 >= 40 & age_20200101 < 50 ~ "40_49",
                                  age_20200101 >= 50 & age_20200101 < 60 ~ "50_59",
                                  age_20200101 >= 60 & age_20200101 < 70 ~ "60_69",
                                  age_20200101 >= 70 & age_20200101 < 80 ~ "70_79",
                                  age_20200101 >= 80 ~ "above_80")) %>% 
        group_by_at(c(strata,'vacc_status')) %>% summarise(n_total=n()) |>
        collect()
      
      input_data_tmp$isoweek <- input_week
      
      input_data_tbl <- bind_rows(input_data_tbl, input_data_tmp)
    }
    
    # save input table for future use
    dbWriteTable(con, input_table_name, input_data_tbl, overwrite = TRUE)
  }
  
  return(input_data_tbl)
}


# Helper functions for multiple regressions -------------------------------

# run model for different variant periods
estimate_model_for_variant_periods_new <- function(strata            = strata,
                                                   con               = con,
                                                   region_id         = region_id,
                                                   age_grouping      = age_grouping,
                                                   IMD_region_data   = IMD_region_data,
                                                   npi_tiers         = npi_tiers,
                                                   min_isoweek       = 17,
                                                   max_isoweek       = 112,
                                                   theta             = NULL,
                                                   lambda            = 0,
                                                   return_prediction = TRUE,
                                                   variant_periods   = c(""), #,"WT","Alpha","Delta", "Omicron"),
                                                   fixed_start_week  = 0,
                                                   db_table          = 'synth_pop_iso',
                                                   keep_cross_validation_predictions = FALSE
)
{
  model <- data <- prediction <- list()
  
  for(variant in variant_periods)
    for(cd in c('isoweek2_death','isoweek2_hosp','isoweek2_pillar2pcr') )   # ,'isoweek2_pillar2'
      model[[paste0(variant,'_',cd)]] <- list('variant'=variant,'case_definition'=cd,'start_end'=paste(min_isoweek,max_isoweek,sep = '_'))
  
  names <- rmse <- mae <- r2 <- AIC <- c()
  results <- fits <- fitted_model <- list()
  for(i in 1:length(model))
  {
    run       <- model[[i]]
    start_end <- as.double(str_split(run$start_end,'_')[[1]])
    print(paste0(unlist(run),collapse =", "))
    
    if(!exists( "use_variant", run)) run$use_variant = FALSE
    if(!exists( "use_tiers", run))   run$use_tiers   = FALSE
    if(!exists( "variant", run))     run$variant     = ""
    
    theta_in = NULL
    if(str_detect(run$case_definition,"pillar"))
    {
      if(run$variant %in% c('Alpha','Delta','Omicron')) theta_in <- NULL
      if(run$variant %in% c('WT')) theta_in = 1/3
    }
    
    #if(str_detect(run$case_definition,"hosp") & run$variant == "" & !is.null(theta)) theta_in = 2
    #if(str_detect(run$case_definition,"death") & run$variant == "" & !is.null(theta)) theta_in = 1
    
    if( sum(c("WT","Alpha","Delta","Omicron") %in% variant_periods) )     # need to deal with the case of switching between varient and tier periods
    {
      run$tier    = ""
    } else {
      run$tier    = run$variant
      run$variant = ""
    }
    
    try({
      fit <- synth_pop_poisson_reg_new( strata            = strata,
                                        case_definition   = run$case_definition,
                                        interaction_pairs = list(),#c(run$case_definition,region_id),
                                        #c(run$case_definition,age_grouping)),
                                        #c(region_id,age_grouping)),
                                        min_isoweek       = start_end[1],
                                        max_isoweek       = start_end[2],
                                        IMD_region_data   = IMD_region_data,
                                        npi_tiers         = npi_tiers,
                                        variant_period    = run$variant,
                                        tier_period       = run$tier,
                                        region_id         = region_id,
                                        con               = con,
                                        theta             = theta_in,
                                        lambda            = lambda,
                                        return_prediction = return_prediction,
                                        fixed_start_week  = fixed_start_week,
                                        restriction_parameteric   = TRUE,
                                        db_table          = db_table,
                                        keep_cross_validation_predictions = keep_cross_validation_predictions)
      results[names(model)[i]] <- h2o.performance(fit$fit, xval = TRUE)
      fits[names(model)[i]]    <- fit$fit #fit@model$coefficients_table %>% mutate(exp_coef = exp(coefficients)) 
      fitted_model[names(model)[i]] <- fit
      
      names[i]  <- names(model)[i]
      res       <- results[[names[i]]]@metrics
      rmse[i]   <- res$RMSE
      mae[i]    <- res$mae
      r2[i]     <- res$r2
      AIC[i]    <- res$AIC
      data[[i]]               <- fit$data
      predictions             <- as_tibble(fit$prediction)
      
      if(keep_cross_validation_predictions)
      {
        cv_prediction           <- combine_cv_predictions(h2o.cross_validation_predictions(fit$fit))
        data[[i]]$cv_prediction <- cv_prediction
      }
      
      data[[i]]$prediction    <- predictions$predict
      data[[i]]$prediction_se <- predictions$StdErr
      data[[i]]$name          <- names(model)[i]
    })
  }
  
  return(list(fits=fits,data=data, fitted_model=fitted_model))
}


# combine cross validation predictions
combine_cv_predictions <- function(cv_preds)
{
  cv_length       <- length(cv_preds)
  n               <- cv_preds[[1]]$predict |> dim()
  preds_matrix    <- matrix(ncol=cv_length,nrow=n[1])
  i               <- 1
  
  for(pred in cv_preds)
  {
    preds_matrix[,i] <- as.vector(pred$predict)
    i <- i+1
  }
  
  return(rowSums(preds_matrix))
}

# create overall summary table
create_summary_table <- function()
{
  IMD_region_data     <- create_IMD_region_data()
  npi_tiers           <- create_UK_restrictions_data()
  region_id           <- 'RGN21NM'          # ltla_code, RGN21NM, ITL221NM
  age_grouping        <- 'age_v2'       # age_band_max
  strata              <- c('sex',region_id,'IMD_national_quintile','ethnicity_simple',age_grouping) 
  max_isoweek         <- 112   #104           27Feb2022 is isoweek 8 of 2022
  min_isoweek         <- 17    #include from week 18
  fixed_start_week    <- 0
  db_table            <- 'synth_pop_multi'
  
  data_pillar2pcr        <- synth_pop_create_data(strata            = strata,
                                                  case_definition   = 'isoweek2_pillar2pcr', 
                                                  min_isoweek       = min_isoweek,       # don't go before May 1st
                                                  max_isoweek       = max_isoweek,
                                                  fixed_start_week  = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                                  variant_period    = '',
                                                  tier_period       = '',
                                                  IMD_region_data   = IMD_region_data,
                                                  npi_tiers         = npi_tiers,
                                                  region_id         = region_id,
                                                  con               = con,
                                                  db_table          = db_table ) 
  
  data_hospitalisation   <- synth_pop_create_data(strata            = strata,
                                                  case_definition   = 'isoweek2_hosp', 
                                                  min_isoweek       = min_isoweek,       # don't go before May 1st
                                                  max_isoweek       = max_isoweek,
                                                  fixed_start_week  = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                                  variant_period    = '',
                                                  tier_period       = '',
                                                  IMD_region_data   = IMD_region_data,
                                                  npi_tiers         = npi_tiers,
                                                  region_id         = region_id,
                                                  con               = con,
                                                  db_table          = db_table ) 
  
  data_death            <- synth_pop_create_data(strata            = strata,
                                                 case_definition   = 'isoweek2_death', 
                                                 min_isoweek       = min_isoweek,       # don't go before May 1st
                                                 max_isoweek       = max_isoweek,
                                                 fixed_start_week  = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                                 variant_period    = '',
                                                 tier_period       = '',
                                                 IMD_region_data   = IMD_region_data,
                                                 npi_tiers         = npi_tiers,
                                                 region_id         = region_id,
                                                 con               = con,
                                                 db_table          = db_table ) 
  
  
  categories <- c('sex','age_v2','RGN21NM','IMD_national_quintile','ethnicity_simple','vacc_status')
  data_summary <- as_tibble(NULL)
  for(cat in categories)
  {
    data_hosp_var <- data_hospitalisation %>% group_by_at(cat) %>% summarise(n_tot=sum(n_total[isoweek2_hosp==98]),#pick a week to get total number of people
                                                                             person_days_at_risk=sum(person_risk_days),
                                                                             n_hosp = sum(n))
    data_pillar2_var <- data_pillar2pcr %>% group_by_at(cat) %>% summarise(n_tot=sum(n_total[isoweek2_pillar2pcr==98]),
                                                                           person_days_at_risk=sum(person_risk_days),
                                                                           n_pillar2pcr = sum(n))
    data_death_var <- data_death %>% group_by_at(cat) %>% summarise(n_tot=sum(n_total[isoweek2_death==98]),
                                                                    person_days_at_risk=sum(person_risk_days),
                                                                    n_death = sum(n))
    data_var <- data_pillar2_var %>% left_join(data_hosp_var) %>% left_join(data_death_var) %>%
      rename(Category = !!cat)
    data_var$Variable <- cat
    
    if(cat=='age_v2')
      data_var$Category <- factor(data_var$Category, levels=c("under_40", "40_49", "50_59", "60_69","70_79","above_80"))
    
    if(cat=='ethnicity_simple')
    {
      data_var <- data_var %>% mutate(Category = case_when(Category=='white' ~ 'White',
                                                           Category=='south_asian' ~ 'South Asian',
                                                           Category=='other_asian' ~ 'Asian (other)',
                                                           Category=='black' ~ 'Black',
                                                           Category=='mixed_other' ~ 'Mixed/Other',
                                                           Category=='mixed' ~ 'Mixed',
                                                           Category=='other' ~ 'Other',
                                                           TRUE ~ Category)) %>%
        arrange(Category)
      
    }
    
    if(cat=='vacc_status')
    {
      data_var <- data_var %>% 
        mutate(Category = str_replace(str_to_title(str_replace_all(str_replace(Category, 'vacc_status.',''),'_',' ')),'Mrna','mRNA')) %>%
        mutate(Category = str_replace(Category,"14[+]D", ">14d")) %>%
        mutate(Category = str_replace(Category,"21[+]D", ">21d")) %>%
        mutate(Category = str_replace(Category,"10 18", "10-18w")) 
      order <- data_var$Category |> unique() |> sort()
      order <- (c('Not Vaccinated',order[order!='Not Vaccinated']))
      
      data_var$Category <- factor( data_var$Category, levels=order)
    }
    
    data_summary <- bind_rows(data_summary,data_var %>% arrange(Category))
  }
  
  data_summary <- data_summary %>% dplyr::select(Variable,Category,n_tot,person_days_at_risk,n_pillar2pcr,n_hosp,n_death) %>%
    rename(n=n_tot,'Cases (Pillar 2 PCR)'=n_pillar2pcr,Hospitalisations=n_hosp,Deaths=n_death, 'Person days at risk'=person_days_at_risk) %>%
    mutate(Variable = case_when(Variable=='sex'~'Sex',
                                Variable=='age_v2'~'Age',
                                Variable=='RGN21NM'~'Region',
                                Variable=='IMD_national_quintile'~'IMD Quintile',
                                Variable=='ethnicity_simple'~'Ethnicity',
                                Variable=='vacc_status'~'Vaccination status',
                                TRUE ~ Variable))
  
  return(data_summary)
}

# create data for map figure
create_summary_map <- function()
{
  IMD_region_data     <- create_IMD_region_data()
  npi_tiers           <- create_UK_restrictions_data()
  age_grouping        <- 'age_v2'       # age_band_max
  max_isoweek         <- 112   #104           27Feb2022 is isoweek 8 of 2022
  min_isoweek         <- 17    #include from week 18
  fixed_start_week    <- 0
  db_table            <- 'synth_pop_multi'
  
  region_id           <- 'ltla_code'
  strata              <- c('sex',region_id,'IMD_national_quintile','ethnicity_simple',age_grouping)
  
  data_pillar2pcr_ltla   <- synth_pop_create_data(strata            = strata,
                                                  case_definition   = 'isoweek2_pillar2pcr', 
                                                  min_isoweek       = min_isoweek,       # don't go before May 1st
                                                  max_isoweek       = max_isoweek,
                                                  fixed_start_week  = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                                  variant_period    = '',
                                                  tier_period       = '',
                                                  IMD_region_data   = IMD_region_data,
                                                  npi_tiers         = npi_tiers,
                                                  region_id         = 'ltla_code',
                                                  con               = con,
                                                  db_table          = db_table )
  data_hosp_ltla         <- synth_pop_create_data(strata            = strata,
                                                  case_definition   = 'isoweek2_hosp', 
                                                  min_isoweek       = min_isoweek,       # don't go before May 1st
                                                  max_isoweek       = max_isoweek,
                                                  fixed_start_week  = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                                  variant_period    = '',
                                                  tier_period       = '',
                                                  IMD_region_data   = IMD_region_data,
                                                  npi_tiers         = npi_tiers,
                                                  region_id         = 'ltla_code',
                                                  con               = con,
                                                  db_table          = db_table )
  data_death_ltla        <- synth_pop_create_data(strata            = strata,
                                                  case_definition   = 'isoweek2_death', 
                                                  min_isoweek       = min_isoweek,       # don't go before May 1st
                                                  max_isoweek       = max_isoweek,
                                                  fixed_start_week  = fixed_start_week,      # don't censor anything prior to this, this effectively becomes start of analysis
                                                  variant_period    = '',
                                                  tier_period       = '',
                                                  IMD_region_data   = IMD_region_data,
                                                  npi_tiers         = npi_tiers,
                                                  region_id         = 'ltla_code',
                                                  con               = con,
                                                  db_table          = db_table )
  
  
  ltla_cases <- data_pillar2pcr_ltla %>% group_by(ltla_code) %>% summarise(n_tot=round(sum(n_total)/(max_isoweek-min_isoweek)),
                                                                           person_days_at_risk=sum(person_risk_days),
                                                                           n_cases = sum(n),
                                                                           attack_rate_pillar2pcr = n_cases / n_tot *100) %>%
    filter(ltla_code != 'E06000053') 
  ltla_hosp  <- data_hosp_ltla %>% group_by(ltla_code) %>% summarise(n_tot=round(sum(n_total)/(max_isoweek-min_isoweek)),
                                                                     person_days_at_risk=sum(person_risk_days),
                                                                     n_hosp = sum(n),
                                                                     risk_hosp = n_hosp / n_tot *100 )
  ltla_death <- data_death_ltla %>% group_by(ltla_code) %>% summarise(n_tot=round(sum(n_total)/(max_isoweek-min_isoweek)),
                                                                      person_days_at_risk=sum(person_risk_days),
                                                                      n_death = sum(n),
                                                                      risk_death = n_death / n_tot * 100 )
  
  ltla_ethnicity <- data_pillar2pcr_ltla %>% group_by(ltla_code,ethnicity_simple) %>% 
    summarise(n_tot = sum(n_total)/(max_isoweek-min_isoweek)) %>%
    pivot_wider(names_from = ethnicity_simple,values_from = n_tot) %>%
    mutate(pct_non_white=1-white/(black+ mixed +other +other_asian+ south_asian + white) )
  
  baseline_ltla <- ltla_cases %>% left_join(ltla_hosp) %>% left_join(ltla_death)
  
  #all age
  ltla_vaccinated <- data_pillar2pcr_ltla %>% group_by(ltla_code) %>% 
    filter(isoweek2_pillar2pcr==111 & (vacc_status == 'not_vaccinated'|str_detect(vacc_status,'first_dose') )) %>%
    summarise(n_unvacc = sum(n_total)) %>%
    left_join(ltla_cases %>% dplyr::select(ltla_code,n_tot),by=c('ltla_code')) %>%
    mutate(pct_vaccinated = 1-n_unvacc/n_tot) %>%
    filter(ltla_code != 'E06000053') 
  
  england_ltla_sf   <- st_read('LAD_DEC_2021_GB_BFC.shp')
  england_region_sf <- st_read('RGN_DEC_2023_EN_BFE.shp')
  england_ltla_sf   <- england_ltla_sf %>% rename(LTLA19CD=LAD21CD,geo_name=LAD21NM) %>%
    right_join(IMD_region_data) %>% filter(!is.na(IMD_national_quintile)) %>%
    left_join(ltla_ethnicity,by=c("LTLA19CD"="ltla_code")) %>%
    left_join(ltla_vaccinated,by=c("LTLA19CD"="ltla_code")) %>%
    left_join(baseline_ltla,by=c("LTLA19CD"="ltla_code"))
  
  figure_1 <- create_map_figure(england_ltla_sf = england_ltla_sf, england_region_sf = england_region_sf)
  
  return(figure_1)
}

# create maps for England
create_map_figure <- function(england_ltla_sf,england_region_sf)
{
  # map1 - rank LTLA IMD quintiles
  england_plt1 <- england_ltla_sf %>% 
    ggplot() +
    geom_sf(data = england_ltla_sf, lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.integer(IMD_national_quintile)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf, lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", na.value = "grey80") +
    #coord_sf(expand = FALSE) +
    theme_void() +
    guides(fill = guide_colorbar(title = "IMD LTLA rank quintile")) +
    theme(legend.position = "bottom", legend.justification = "left")
  
  london_plt1 <- england_ltla_sf %>% filter(RGN21CD=='E12000007') %>%
    ggplot() +
    geom_sf(data = england_ltla_sf %>% filter(RGN21CD=='E12000007'), lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.integer(IMD_national_quintile)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf %>% filter(RGN23CD=='E12000007'), lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80") +
    theme_void() +
    guides(fill = guide_colorbar(title = "IMD LTLA rank quintile")) +
    theme(legend.position = "none", legend.justification = "left",
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ggtitle('London')
  
  plt1 <- england_plt1 + labs(tag = "A") + 
    inset_element(london_plt1, left = 0.025, bottom = 0.4, right = 0.325, top = 0.8, ignore_tag = TRUE)
  
  # map2 - percentage ethnicity non-white by LTLA
  min <- min(na.omit(england_ltla_sf$pct_non_white))
  max <- max(na.omit(england_ltla_sf$pct_non_white))
  
  england_plt2 <- england_ltla_sf %>% 
    ggplot() +
    geom_sf(data = england_ltla_sf, lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(pct_non_white)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf, lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(min-0.02,max+0.02)) +
    #coord_sf(expand = FALSE) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Percentage non-white population")) +
    theme(legend.position = "bottom", legend.justification = "left")
  
  london_plt2 <- england_ltla_sf %>% filter(RGN21CD=='E12000007') %>%
    ggplot() +
    geom_sf(data = england_ltla_sf %>% filter(RGN21CD=='E12000007'), lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(pct_non_white)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf %>% filter(RGN23CD=='E12000007'), lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(min-0.02,max+0.02)) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Percentage non-white population")) +
    theme(legend.position = "none", legend.justification = "left",
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ggtitle('London')
  
  plt2 <- england_plt2 + labs(tag = "B") + 
    inset_element(london_plt2, left = 0.025, bottom = 0.4, right = 0.325, top = 0.8, ignore_tag = TRUE)
  
  
  # map3 - Log Mean Population density by LTLA
  min <- min(na.omit(england_ltla_sf$pct_vaccinated))
  max <- max(na.omit(england_ltla_sf$pct_vaccinated))
  
  england_plt3 <- england_ltla_sf %>% 
    ggplot() +
    geom_sf(data = england_ltla_sf, lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(pct_vaccinated)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf, lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(min-0.02,max+0.02)) +
    #coord_sf(expand = FALSE) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Vaccination percentage \n(at least 2 doses)")) +
    theme(legend.position = "bottom", legend.justification = "left")
  
  london_plt3 <- england_ltla_sf %>% filter(RGN21CD=='E12000007') %>%
    ggplot() +
    geom_sf(data = england_ltla_sf %>% filter(RGN21CD=='E12000007'), lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(pct_vaccinated)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf %>% filter(RGN23CD=='E12000007'), lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(min-0.02,max+0.02)) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Vaccination percentage \n(at least 2 doses)")) +
    theme(legend.position = "none", legend.justification = "left",
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ggtitle('London')
  
  plt3 <- england_plt3 + labs(tag = "C") + 
    inset_element(london_plt3, left = 0.025, bottom = 0.4, right = 0.325, top = 0.8, ignore_tag = TRUE)
  
  # map4 - Attack rate infection (pillar2pcr) by LTLA 
  min <- min(na.omit(england_ltla_sf$attack_rate_pillar2pcr))
  max <- max(na.omit(england_ltla_sf$attack_rate_pillar2pcr))
  
  england_plt4 <- england_ltla_sf %>% 
    ggplot() +
    geom_sf(data = england_ltla_sf, lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(attack_rate_pillar2pcr)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf, lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(min-0.002,max+0.002)) +
    #coord_sf(expand = FALSE) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Attack rate Pillar 2 PCR positive cases")) +
    theme(legend.position = "bottom", legend.justification = "left")
  
  london_plt4 <- england_ltla_sf %>% filter(RGN21CD=='E12000007') %>%
    ggplot() +
    geom_sf(data = england_ltla_sf %>% filter(RGN21CD=='E12000007'), lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(attack_rate_pillar2pcr)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf %>% filter(RGN23CD=='E12000007'), lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(min-0.002,max+0.002)) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Attack rate Pillar 2 PCR positive cases")) +
    theme(legend.position = "none", legend.justification = "left",
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ggtitle('London')
  
  plt4 <- england_plt4 + labs(tag = "D") + 
    inset_element(london_plt4, left = 0.025, bottom = 0.4, right = 0.325, top = 0.8, ignore_tag = TRUE) 
  
  # map5 - relative risk of hospitalisation (pillar2pcr) by LTLA
  min <- min(na.omit(england_ltla_sf$risk_hosp))
  max <- max(na.omit(england_ltla_sf$risk_hosp))
  
  england_plt5 <- england_ltla_sf %>% 
    ggplot() +
    geom_sf(data = england_ltla_sf, lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(risk_hosp)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf, lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(0,max+0.0000002)) +
    #coord_sf(expand = FALSE) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Attack rate hospitalisation")) +
    theme(legend.position = "bottom", legend.justification = "left")
  
  london_plt5 <- england_ltla_sf %>% filter(RGN21CD=='E12000007') %>%
    ggplot() +
    geom_sf(data = england_ltla_sf %>% filter(RGN21CD=='E12000007'), lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(risk_hosp)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf %>% filter(RGN23CD=='E12000007'), lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(0,max+0.0000002)) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Attack rate hospitalisation")) +
    theme(legend.position = "none", legend.justification = "left",
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ggtitle('London')
  
  plt5 <- england_plt5 + labs(tag = "E") + 
    inset_element(london_plt5, left = 0.025, bottom = 0.4, right = 0.325, top = 0.8, ignore_tag = TRUE) 
  
  # map6 - relative risk of death (pillar2pcr) by LTLA
  min <- min(na.omit(england_ltla_sf$risk_death))
  max <- max(na.omit(england_ltla_sf$risk_death))
  
  england_plt6 <- england_ltla_sf %>% 
    ggplot() +
    geom_sf(data = england_ltla_sf, lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(risk_death)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf, lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(0,max+0.0000002)) +
    #coord_sf(expand = FALSE) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Attack rate death")) +
    theme(legend.position = "bottom", legend.justification = "left")
  
  london_plt6 <- england_ltla_sf %>% filter(RGN21CD=='E12000007') %>%
    ggplot() +
    geom_sf(data = england_ltla_sf %>% filter(RGN21CD=='E12000007'), lwd = 0.3, col = "grey60") +
    geom_sf(aes(fill = as.double(risk_death)), lwd = 0.3, col = "grey20") +
    geom_sf(data = england_region_sf %>% filter(RGN23CD=='E12000007'), lwd = 0.7, col = "black",  fill = NA) +
    scale_fill_hp(option = "NewtScamander", 
                  na.value = "grey80", 
                  limits = c(0,max+0.0000002)) +
    theme_void() +
    guides(fill = guide_colorbar(title = "Attack rate death")) +
    theme(legend.position = "none", legend.justification = "left",
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    ggtitle('London')
  
  plt6 <- england_plt6 + labs(tag = "F") + 
    inset_element(london_plt6, left = 0.025, bottom = 0.4, right = 0.325, top = 0.8, ignore_tag = TRUE) 
  
  design <- "
  12
  34
  56
"
  figure1_plot <- plt1 + plt2 + plt3 + plt4 + plt5 + plt6 + 
    plot_layout(design = design, tag_level = 'keep') + 
    plot_annotation(tag_levels = 'A')
  
  ggsave("figure_1.png", plot = figure1_plot, width = 12, height = 18)
  
  return(figure1_plot)
}

