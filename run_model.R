source('R/load_packages.R')
source('R/estimation_fns.R')
source('R/helper_fns.R')
source('R/data_fns.R')
source('R/plot_fns.R')

#### Data access ####
# While all data used in this analysis were anonymised, the individual-level nature 
# of the data used risks individuals being identified or being able to self-identify 
# if it is released publicly.  Requests for access to the underlying source data 
# should be directed to UKHSA.

if(file.exists("synth_pop_census2021.duckdb"))
{
  # load data / connect to duckdb table -- THE DATA IS ONLY AVAILABLE ON REQUEST AND NOT PROVIDED
  con <- dbConnect(duckdb(), dbdir = "synth_pop_census2021.duckdb", read_only = FALSE)
  
  # Table 1 - Table of synthetic population 
  data_summary <- create_summary_table( db_table = 'synth_pop_booster')     # in paper split into A and B panel, with B vaccination data
  write_csv(data_summary,'figures/data_summary.csv')
  
  # Figure 2 - map of England
  figure_1 <- create_summary_map(db_table = 'synth_pop_booster')
  
  # Figure 3 (ethnicity & IMD) + Figure 4 (vaccination)
  time_periods        <- c("","WT","Alpha","Delta","Omicron")
  region_id           <- 'RGN21NM'
  age_grouping        <- 'age_v2'
  strata              <- c('sex',region_id,'IMD_national_quintile','ethnicity_simple',age_grouping) 
  max_isoweek         <- 112
  min_isoweek         <- 17
  IMD_region_data     <- create_IMD_region_data()
  npi_tiers           <- create_UK_restrictions_data()
  
  model_fit    <- estimate_model_for_variant_periods_new(strata = strata,
                                                         con               = con,
                                                         region_id         = region_id,
                                                         age_grouping      = age_grouping,
                                                         IMD_region_data   = region_id,
                                                         npi_tiers         = npi_tiers,
                                                         min_isoweek       = min_isoweek,
                                                         max_isoweek       = max_isoweek,
                                                         theta             = 3,
                                                         variant_periods   = time_periods,  # only run for full period of vacc
                                                         fixed_start_week  = 0,
                                                         db_table = 'synth_pop_booster' )    # multi
  
  model_fit$fits |> saveRDS('data/main_model_fits.RDS')
  
  ethnicity <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr')) #+ theme(legend.position = 'none')
  IMD       <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='IMD') + theme(legend.position = 'none')
  sex       <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='sex') #+ theme(legend.position = 'bottom')
  age       <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='age_v2') #+ theme(legend.position = 'bottom')
  vacc_status <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='vacc_status', min_p_value = 0.05,restrict_vacc_plot = TRUE,plot_VE = TRUE) #+ theme(legend.position = 'bottom')
  res       <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='restriction', min_p_value = 0.05) #+ theme(legend.position = 'bottom')
  region    <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='RGN21NM', min_p_value = 0.05) #+ theme(legend.position = 'bottom')
  
  vacc_full    <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='vacc_status', min_p_value = 0.05,restrict_vacc_plot = FALSE,plot_VE = FALSE) #+ theme(legend.position = 'bottom')
  vacc_full_VE <- plot_IRR_results(model_fit$fits,include = c('death','hosp','pillar2pcr'),filter_by='vacc_status', min_p_value = 0.05,restrict_vacc_plot = FALSE,plot_VE = TRUE) #+ theme(legend.position = 'bottom')
  
  #IMD$data |> filter(variant=='all'& names=='IMD5') |> View()
  
  plot_variant_period_simple <- (IMD / ethnicity) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
  plot_variant_period_SI <- ((sex + age) / region  ) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
  plot_vacc_SI <- vacc_full / vacc_full_VE + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
  ggsave("figures/figure3.png", plot = plot_variant_period_simple, width = 12, height = 8)
  ggsave("figures/figure3_SI.png", plot = plot_variant_period_SI, width = 12, height = 8)
  ggsave("figures/figure_vaccination.png", plot = vacc_status, width = 12, height = 5)
  ggsave("figures/figureSI_vaccination.png", plot = plot_vacc_SI, width = 12, height = 12)
  ggsave("figures/figure_restriction.png", plot = res, width = 12, height = 2)
  
  # SI analysis -------------------------------------------------------------
  
  layout_SI_IMD <- "
12
12
34
34
56
56
56
77
77
77
77
77
77"
  
  IMD_subcomponents <- c('IMD_national_quintile','Income_national_quintile','Employment_national_quintile','Education_national_quintile','Health_national_quintile','Crime_national_quintile','Housing_national_quintile','Environment_national_quintile')
  subcomponent_fits <- list()
  
  for(IMD_subcomponent in IMD_subcomponents)
  {
    time_periods <- c("","WT","Alpha","Delta","Omicron")
    strata              <- c('sex',region_id,IMD_subcomponent,'ethnicity_simple',age_grouping) 
    max_isoweek         <- 112
    min_isoweek         <- 17
    
    model_fit_subcomp <- estimate_model_for_variant_periods_new(strata = strata,
                                                                con               = con,
                                                                region_id         = region_id,
                                                                age_grouping      = age_grouping,
                                                                IMD_region_data   = IMD_region_data,
                                                                npi_tiers         = npi_tiers,
                                                                min_isoweek       = min_isoweek,
                                                                max_isoweek       = max_isoweek,
                                                                theta             = 3,
                                                                variant_periods   = time_periods,  
                                                                fixed_start_week  = 0,
                                                                db_table          = 'synth_pop_booster')    
    
    subcomponent_fits[[IMD_subcomponent]] <- model_fit_subcomp
    
    ethnicity <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr')) + theme(legend.position = 'none')
    IMD       <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr'),filter_by=str_replace(IMD_subcomponent,"_national_quintile","")) + theme(legend.position = 'none')
    sex       <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr'),filter_by='sex') #+ theme(legend.position = 'bottom')
    age       <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr'),filter_by='age_v2') #+ theme(legend.position = 'bottom')
    vacc      <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr'),filter_by='vacc_status', min_p_value = 0.05)
    res       <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr'),filter_by='restriction', min_p_value = 0.05)   
    region    <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr'),filter_by='RGN21NM', min_p_value = 0.05)
    
    plot_variant_subcomp <- IMD + ethnicity + sex + age + region + res +  vacc + plot_layout(design=layout_SI_IMD,guides = "collect") + plot_annotation(tag_levels = 'A')
    ggsave(paste0("figures/figureSI_", str_to_lower(str_replace(IMD_subcomponent,"_national_quintile","")),".png"), plot = plot_variant_subcomp, width = 12, height = 14)
  }
  
  IMD_plots <- list()
  for(comp in IMD_subcomponents)
  {
    model_fit_subcomp <- subcomponent_fits[[comp]]  
    IMD_plots[[comp]]       <- plot_IRR_results(model_fit_subcomp$fits,include = c('death','hosp','pillar2pcr'),filter_by=str_replace(comp,"_national_quintile","")) + theme(legend.position = 'none')
  }
  
  IMD_subcomp_plot <- patchwork::wrap_plots(IMD_plots,ncol=1) + plot_annotation(tag_levels = 'A')
  ggsave(paste0("figures/figureSI_IMD_subcomponents",".png"), plot = IMD_subcomp_plot, width = 12, height = 18)
  
  # Table for vaccination IRRs + vaccine effectiveness 
  vacc_csv <- vaccination_table(model_fit$fits, min_p_value = 0.05 )
  write.table(vacc_csv,'figures/vaccination_effectiveness.csv', sep = ',',
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  ethnicity_csv <- ethnicity_table(subcomponent_fits[['IMD_national_quintile']]$fits, min_p_value = 1 )
  write.table(ethnicity_csv,'figures/ethnicity_IRR.csv', sep = ',',
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  IMD_csv <- IMD_table(subcomponent_fits[['IMD_national_quintile']]$fits, min_p_value = 1 )
  write.table(IMD_csv,'figures/IMD_IRR.csv', sep = ',',
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  # Model choice table
  model_selection_tbl <- as_tibble(NULL)
  files <- c('model_selection_NEW.RDS','model_selection_WT.RDS','model_selection_ALPHA.RDS',
             'model_selection_DELTA.RDS','model_selection_OMICRON.RDS')
  
  for(file in files)
  {
    tmp    <- readRDS(paste0('data/',file))
    period <- str_to_title(str_replace(str_split(file,'[.]')[[1]][1],'model_selection_',''))
    
    if(period=='New') period <- 'Full time period'
    tmp$time_period <- period
    
    model_selection_tbl <- bind_rows(model_selection_tbl,tmp)
  }
  
  model_selection_preprocess <- model_selection_tbl %>% filter(!is.na(index)) %>% 
    mutate(AIC = replace_na(AIC,0), 
           R2  = replace_na(R2,0)) %>%
    filter(R2>0) %>%
    dplyr::select(-c(vaccination,restriction,index,region)) %>%  # take out region as all RGN21NM
    mutate(interaction = case_when(interaction==1~'none',
                                   interaction==2~'(WxR)',
                                   interaction==3~'(WxA)',
                                   interaction==4~'(WxR, WxA)',
                                   interaction==5~'(AxR)',
                                   interaction==6~'(AxR, WxA)',
                                   interaction==7~'(WxR, AxR)',
                                   interaction==8~'(WxRxA)',
                                   interaction==9~'(AxV)',
                                   interaction==10~'(WxA,AxV)'),
           case_definition = case_when(case_definition=='isoweek2_death'~'Death',
                                       case_definition=='isoweek2_hosp' ~'Hospitalisation',
                                       case_definition=='isoweek2_pillar2pcr' ~ 'Pillar 2 PCR case',
                                       TRUE ~ case_definition),
           age = case_when(age=='age_v2'~'Age simplified',
                           age=='age_band_max' ~ 'Age (10yr bands)',
                           TRUE ~ age),
           deprivation = case_when(deprevation=="IMD_national_quintile" ~ 'IMD quintiles',
                                   deprevation=="IMD_national_decile" ~ 'IMD deciles',
                                   TRUE ~ deprevation),
           ethnicity = case_when(db_table=="synth_pop_iso" ~ '4 groups',
                                 db_table=="synth_pop_multi" ~ '5 groups',
                                 TRUE ~ db_table))%>%
    arrange(case_definition,desc(-AIC)) %>% 
    dplyr::select(case_definition,time_period,age,deprivation,ethnicity,interaction,RMSE,MAE,R2,AIC) %>% 
    rename(`Case definition`=case_definition,`Time period`=time_period,Age=age,Deprivation=deprivation,Ethnicity=ethnicity,Interactions=interaction) 
  
  model_selection_preprocess |> saveRDS('data/model_selection_preprocess.RDS')
  
  create_model_selection_gt(model_selection_preprocess,'Death')
  create_model_selection_gt(model_selection_preprocess,'Hospitalisation')
  create_model_selection_gt(model_selection_preprocess,'Pillar 2 PCR case')
  
  # NEJM rep ----------------------------------------------------------------
  
  time_periods <- c("Alpha","Delta")
  strata              <- c('sex',region_id,'IMD_national_quintile','ethnicity_simple',age_grouping) 
  
  model_fit_nejm <- estimate_model_for_variant_periods_new(strata = strata,
                                                           con               = con,
                                                           region_id         = region_id,
                                                           age_grouping      = age_grouping,
                                                           IMD_region_data   = IMD_region_data,
                                                           npi_tiers         = npi_tiers,
                                                           min_isoweek       = 65,
                                                           max_isoweek       = 73,
                                                           theta             = 3,
                                                           variant_periods   = time_periods,  
                                                           fixed_start_week  = 0,
                                                           db_table          = 'synth_pop_multi')  
  
  ethnicity <- plot_IRR_results(model_fit_nejm$fits,include = c('death','hosp','pillar2pcr'),min_p_value = 0.05) + theme(legend.position = 'none')
  IMD       <- plot_IRR_results(model_fit_nejm$fits,include = c('death','hosp','pillar2pcr'),min_p_value = 0.05,filter_by='IMD') + theme(legend.position = 'none')
  vacc      <- plot_IRR_results(model_fit_nejm$fits,include = c('death','hosp','pillar2pcr'),filter_by='vacc_status', min_p_value = 0.05)
  
  ggsave("figures/figureSI_nejm_vacc.png", plot = vacc, width = 12, height = 5)
  
  nejm_vacc_table <- vaccination_table(model_fit_nejm$fits, min_p_value = 0.05 ) %>%
    filter(!str_detect(Vaccine,'Mixed Dose')&!str_detect(Vaccine,'Wane'))
  
  write.table(nejm_vacc_table,'figures/nejm_vacc_table.csv', sep = ',',
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

