

plot_IRR_results <- function(fits, 
                             filter_by = 'ethnicity', 
                             include = c('death','hosp','pillar2','pillar2pcr'),
                             min_p_value = 1,
                             restriction_levels = FALSE,
                             restrict_vacc_plot = FALSE,
                             plot_VE = FALSE)
{
  number_of_fits <- length(fits)
  names_of_fits  <- names(fits)
  
  results_to_plot <- tibble(NULL)
  
  for(i in 1:number_of_fits)
  {
    details <- str_split(names_of_fits[i],'_')[[1]]
    
    if(details[1]=='level'|details[1]=='national')
    {
      details[1] = paste(details[1],details[2],sep="_")
      details    = details[-2]
    }
    
    if(details[1]=='') variant = 'all' else variant = details[1]
    case_definition = details[3]
    
    fit <- fits[[i]]@model$coefficients_table %>% mutate(exp_coef = exp(coefficients)) %>%
      filter(str_detect(names, filter_by)) %>% filter(p_value < min_p_value|is.na(p_value)) # need to include NA as this is for the 
    # reference values
    if(dim(fit)[1]==0) next
    
    # need to add all the reference groups back in!
    if(str_detect('ethnicity_simple',filter_by))
    {
      tmp <- as_tibble(list(names='ethnicity_simple.white',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    if(filter_by %in% c('IMD','Income','Employment','Education','Health','Crime','Housing','Environment','IMD_national_decile')) # str_detect('IMD',filter_by))
    {
      tmp <- as_tibble(list(names=paste(filter_by,'1'),coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    if(str_detect('sex',filter_by))
    {
      tmp <- as_tibble(list(names='sex.F',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    if(str_detect('vacc_status',filter_by))
    {
      tmp <- as_tibble(list(names='vacc_status.not_vaccinated',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    if(str_detect('age_v2',filter_by))
    {
      tmp <- as_tibble(list(names='age_v2.50_59',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    if(str_detect('RGN21NM',filter_by))
    {
      tmp <- as_tibble(list(names='RGN21NM.London',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    fit$variant <- variant
    fit$case_definition <- case_definition
    
    results_to_plot <- bind_rows(results_to_plot,fit)
  }
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  if(restriction_levels)
    results_to_plot$variant <- factor(results_to_plot$variant,levels=c("none","level_1","level_2","level_3","national_lockdown")) 
  
  
  # other groupings
  if(filter_by == 'age')
  {
    results_to_plot$names <- factor(results_to_plot$names, levels=c("age_band_max.5", "age_band_max.17", "age_band_max.24",
                                                                    "age_band_max.34", "age_band_max.44", "age_band_max.54", 
                                                                    "age_band_max.64", "age_band_max.74", "age_band_max.84", "age_band_max.130"))
  }
  
  if(filter_by == 'age_simple')
  {
    results_to_plot$names <- factor(results_to_plot$names, levels=c("age_simple.Children", "age_simple.Adults", "age_simple.old_age"))
  }
  
  if(filter_by == 'age_v2')
  {
    results_to_plot$names <- str_replace(results_to_plot$names,'age_v2.','')
    results_to_plot$names <- factor(results_to_plot$names, levels=c("under_40", "40_49", "50_59", "60_69",
                                                                    "70_79","above_80"))
  }
  
  if(filter_by == 'RGN21NM')
  {
    results_to_plot$names <- str_replace(results_to_plot$names,'RGN21NM.','')
  }
  
  if(filter_by == 'weeks_since_vacc')
  {
    results_to_plot$names <- factor(results_to_plot$names, levels=c("E1_weeks_since_vacc.unvaccinated", "E1_weeks_since_vacc.2-6 weeks", 
                                                                    "E1_weeks_since_vacc.6-10 weeks", "E1_weeks_since_vacc.10-18 weeks", 
                                                                    "E1_weeks_since_vacc.18+ weeks"))
  }
  
  if(filter_by == 'restriction')
  {
    results_to_plot$variant <- factor(results_to_plot$variant, levels=c("all",'WT','Alpha',
                                                                        'Delta','Omicron'))
  }
  
  if(filter_by == 'ethnicity')
  {
    results_to_plot$names <- str_replace(results_to_plot$names,'ethnicity_simple.','')
    
    results_to_plot <- results_to_plot %>% mutate(names = case_when(names=='white' ~ 'White',
                                                                    names=='south_asian' ~ 'South Asian',
                                                                    names=='other_asian' ~ 'Asian (other)',
                                                                    names=='black' ~ 'Black',
                                                                    names=='mixed_other' ~ 'Mixed/Other',
                                                                    names=='mixed' ~ 'Mixed',
                                                                    names=='other' ~ 'Other',
                                                                    TRUE ~ names))
    
    results_to_plot$names <- factor(results_to_plot$names, levels=c("White", 'South Asian', 'Asian (other)','Black',
                                                                    'Mixed/Other', 'Mixed', 'Other', "unknown", 
                                                                    "south_asian","other_asian",
                                                                    "asian", "black", "mixed_other"))
  }
  
  if(filter_by == 'IMD' & sum(str_detect(results_to_plot$names, "prop_LSOA") ))
  {
    seq_1 <- seq(0,length(unique(results_to_plot$names))-1)*10
    seq_2 <- seq(1,length(unique(results_to_plot$names)))*10
    group <- as.integer(substrRight(results_to_plot$names,1))
    results_to_plot$names <- paste0("IMD_prop_LSOA_measure.",group*10,"%-",(group+1)*10, "%")  
    results_to_plot$names <- factor(results_to_plot$names, levels=paste0("IMD_prop_LSOA_measure.",seq_1,"%-",seq_2, "%"))
  }
  
  if( (filter_by == 'IMD' | filter_by == 'IMD_national_decile') & length(unique(results_to_plot$names)) > 9 & !sum(str_detect(results_to_plot$names, "prop_LSOA") ))
  {
    results_to_plot$names <- str_replace(results_to_plot$names,'IMD_national_decile.','')
    results_to_plot$names <- factor(results_to_plot$names, levels=paste0("IMD",seq(1,length(unique(results_to_plot$names)))))
  }
  
  if(filter_by == 'IMD' & sum(str_detect(results_to_plot$names, 'quintile')) & !sum(str_detect(results_to_plot$names, "prop_LSOA") ))
  {
    results_to_plot$names <- str_replace(results_to_plot$names,'IMD_national_quintile.','')
  }
  
  if(filter_by %in% c('Income','Employment','Housing','Crime','Health','Environment','Education') & sum(str_detect(results_to_plot$names, 'quintile')))
  {
    results_to_plot$names <- str_replace(results_to_plot$names,paste0(filter_by,'_national_quintile.'),'')
  }
  
  if(filter_by == 'vaccine_v2')
  {
    results_to_plot$names <- str_replace(results_to_plot$names,'vaccine_v2.','')
  }
  
  if(filter_by == 'Moran' & length(unique(results_to_plot$names)) > 9 )
  {
    results_to_plot$names <- factor(results_to_plot$names, levels=paste0("Moran_I_national_decile.Moran_I_",seq(1,length(unique(results_to_plot$names)))))
  }
  
  if(filter_by == 'sex' )
  {
    results_to_plot <- results_to_plot %>% mutate(names = case_when(names=='sex.F'~'Female',
                                                                    names=='sex.M'~'Male'))
  }
  
  #if(filter_by %in% c('IMD','age','Moran','age_simple','age_v2','Income', 'Employment','Housing','Crime','Health', 'Environment','Education'))
  if(filter_by %in% c('age','Moran','age_simple','age_v2'))
  {
    plt <- results_to_plot %>% #group_by(names,variant) %>%
      filter(!is.na(exp_coef) & !is.na(names) & case_definition %in% include ) %>%
      mutate(case_definition = case_when(case_definition=='hosp'~'Hospitalisations',
                                         case_definition=='death'~'Deaths',
                                         case_definition=='pillar2pcr'~ 'Pillar 2 PCR Cases')) %>%
      mutate(lower_bound = exp(coefficients-1.96*std_error),
             upper_bound = exp(coefficients+1.96*std_error)) %>%
      ggplot(aes(x=names,y=exp_coef,color=variant,group=variant)) + geom_line() +
      geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound,fill=variant), alpha = 0.2) +
      facet_wrap(~case_definition,scales = "free_y") + ylab('IRR') + xlab('') + theme_bw() +
      #scale_color_hp_d(option='NewtScamander') +
      scale_color_lancet() +
      scale_fill_lancet() +
      #scale_color_viridis(discrete = TRUE, end=0.9) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
            plot.tag = element_text(face = 'bold')) 
  } else if (filter_by %in% c('vacc_status','restriction','RGN21NM')) {
    
    results_to_plot <- results_to_plot %>% 
      filter(!is.na(exp_coef) & !is.na(names) & case_definition %in% include ) %>%
      mutate(names = str_replace(str_to_title(str_replace_all(str_replace(names, 'vacc_status.',''),'_',' ')),'Mrna','mRNA')) %>%
      mutate(names = str_replace(names,"14[+]D", ">2w")) %>%
      mutate(names = str_replace(names,"21[+]D", ">3w")) %>%
      mutate(names = str_replace(names,"14d", "2w")) %>%
      mutate(names = str_replace(names,"10 18", "10-18w")) %>%
      mutate(names = str_replace(names,"Over 18", "over 18w")) %>%
      mutate(names = str_replace(names,"Wane", "Second Dose")) %>%
      mutate(names = str_replace(names,'Second Dose >14d', 'Second Dose 2-10w'))
    
    results_to_plot <- results_to_plot %>% mutate(names = str_replace(names,'Booster Dose','Booster'),
                                                  names = str_replace(names,'Booster','Booster Dose') ) #ensure we always call this a booster dose
                                                  
    order        <- c('Not Vaccinated')
    names_unique <- results_to_plot$names %>% unique()
    for(type in c('First','Second','Booster'))
    {
      names_tmp <- sort( names_unique[ str_starts(names_unique, type) ] )
      order     <- c(order, names_tmp)
    }
    #order <- results_to_plot$names |> unique() |> sort()
    #order <- rev(c('Not Vaccinated',order[order!='Not Vaccinated']))
    order <- rev(order)
    
    results_to_plot$names <- factor(results_to_plot$names, levels=order)
    
    if(filter_by=='vacc_status' & restrict_vacc_plot)
    {
      results_to_plot <- results_to_plot %>% filter(!str_detect(names,'First Dose') & 
                                                      !str_detect(names,'Novavax') &
                                                      !str_detect(names,'Mixed Dose'))
    }
    
    y_label <- 'IRR'
    if(plot_VE)
    {
      y_label <- 'Vaccine effectiveness'
      results_to_plot <- results_to_plot %>% filter(names!='Not Vaccinated') %>%
        mutate(case_definition = case_when(case_definition=='hosp'~'Hospitalisations',
                                           case_definition=='death'~'Deaths',
                                           case_definition=='pillar2pcr'~ 'Pillar 2 PCR Cases')) %>%
        mutate(exp_coef    = 1 - exp_coef,
               lower_bound = 1 - exp(coefficients+1.96*std_error),
               upper_bound = 1 - exp(coefficients-1.96*std_error))
    } else {
      results_to_plot <- results_to_plot %>% 
        mutate(case_definition = case_when(case_definition=='hosp'~'Hospitalisations',
                                           case_definition=='death'~'Deaths',
                                           case_definition=='pillar2pcr'~ 'Pillar 2 PCR Cases')) %>%
        mutate(lower_bound = exp(coefficients-1.96*std_error),
               upper_bound = exp(coefficients+1.96*std_error))
    }
    
    x_label <- ''
    if(filter_by=='restriction')
    {
      results_to_plot <- results_to_plot %>% mutate(names = '')
      x_label         <- 'Restriction Level'
    }
    
    plt <- results_to_plot %>% 
      ggplot(aes(x=names,y=exp_coef,color=variant,group=variant)) + 
      geom_point(position=position_dodge(width=0.4)) +
      #scale_color_hp_d(option='NewtScamander') +
      scale_color_lancet() +
      #scale_color_viridis(discrete = TRUE, end=0.9) +
      geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound),alpha=0.5,position=position_dodge(width=0.4),width=0.4) +
      coord_flip() +
      facet_wrap(~case_definition,scales = "free_x") + ylab(y_label) + xlab(x_label) + theme_bw() +
      theme(plot.tag = element_text(face = 'bold'))
    
    if(filter_by=='vacc_status' & plot_VE )
      plt <- plt + ylim(.2, 1)
    
    if(filter_by=='vacc_status' & !plot_VE )
      plt <- plt + ylim(0, 1.4)
    
  } else {
    plt <- results_to_plot %>% #group_by(names,variant) %>%
      filter(!is.na(exp_coef) & !is.na(names) & case_definition %in% include) %>%
      mutate(case_definition = case_when(case_definition=='hosp'~'Hospitalisations',
                                         case_definition=='death'~'Deaths',
                                         case_definition=='pillar2pcr'~ 'Pillar 2 PCR Cases')) %>%
      mutate(lower_bound = exp(coefficients-1.96*std_error),
             upper_bound = exp(coefficients+1.96*std_error)) %>%
      ggplot(aes(x=names,y=exp_coef,color=variant,group=variant)) + 
      geom_point(position=position_dodge(width=0.4)) +
      #scale_color_hp_d(option='NewtScamander') +
      scale_color_lancet() +
      #scale_color_viridis(discrete = TRUE, end=0.9) +
      geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound),alpha=0.5,position=position_dodge(width=0.4),width=0.4) +
      facet_wrap(~case_definition,scales = "free_y") + ylab('IRR') + xlab('') + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
            plot.tag = element_text(face = 'bold')) 
  }
  
  return(plt)
}




# table functions ---------------------------------------------------------


insert_blank_rows_variant <- function(dataframe, column) {
  plain_labels <- unique(dataframe[[column]])
  tlabels      <- as.character(unique(dataframe[[column]]))
  
  pre_pend <- letters[1:length(tlabels)]
  
  for(i in 1:length(tlabels))
    tlabels[i] <- paste(pre_pend[i],tlabels[i])
  
  dataframe <- dataframe %>% mutate(variant = case_when(variant=='Full time period' ~ "a Full time period",
                                                        variant=='WT' ~ 'b WT',
                                                        variant=='Alpha' ~ 'c Alpha',
                                                        variant=='Delta' ~ 'd Delta',
                                                        variant=='Omicron' ~ 'e Omicron'))
  
  df_split <- split(dataframe, dataframe[[column]])
  dataframe <- do.call(rbind, lapply(df_split, function(group) {
    rbind(group, rep(NA, ncol(group)))
  }))
  dataframe <- rbind(NA,dataframe)#add NAs for top row header
  dataframe <- dataframe[-nrow(dataframe),]#remove NAs at bottom
  dataframe[[column]] <- NULL#remove column
  inds                 <- which(!complete.cases(dataframe))
  dataframe[inds,1]    <- plain_labels
  dataframe[[1]][inds] <- paste0("\\bfseries{", dataframe[[1]][inds], "}")
  dataframe            <- dataframe %>% mutate_all(~ ifelse(is.na(.), "", .))
  
  return(dataframe)
}



vaccination_table <- function(fits,
                              min_p_value = 1 )
{
  filter_by = 'vacc_status'
  include = c('death','hosp','pillar2pcr')
  
  number_of_fits <- length(fits)
  names_of_fits  <- names(fits)
  
  results_to_plot <- tibble(NULL)
  
  for(i in 1:number_of_fits)
  {
    details <- str_split(names_of_fits[i],'_')[[1]]
    
    if(details[1]=='level'|details[1]=='national')
    {
      details[1] = paste(details[1],details[2],sep="_")
      details    = details[-2]
    }
    
    if(details[1]=='') variant = 'all' else variant = details[1]
    case_definition = details[3]
    
    fit <- fits[[i]]@model$coefficients_table %>% mutate(exp_coef = exp(coefficients)) %>%
      filter(str_detect(names, filter_by)) %>% filter(p_value < min_p_value|is.na(p_value)) # need to include NA as this is for the 
    # reference values
    if(dim(fit)[1]==0) next
    
    # need to add all the reference groups back in!
    if(str_detect('vacc_status',filter_by))
    {
      tmp <- as_tibble(list(names='vacc_status.not_vaccinated',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    fit$variant <- variant
    fit$case_definition <- case_definition
    
    results_to_plot <- bind_rows(results_to_plot,fit)
  }
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  results_to_plot <- results_to_plot %>% 
    filter(!is.na(exp_coef) & !is.na(names) & case_definition %in% include ) %>%
    mutate(names = str_replace(str_to_title(str_replace_all(str_replace(names, 'vacc_status.',''),'_',' ')),'Mrna','mRNA')) %>%
    mutate(names = str_replace(names,"14[+]D", ">2w")) %>%
    mutate(names = str_replace(names,"21[+]D", ">3w")) %>%
    mutate(names = str_replace(names,"14d", "2w")) %>%
    mutate(names = str_replace(names,"21d", "3w")) %>%
    mutate(names = str_replace(names,"10 18", "10-18w")) %>%
    mutate(names = str_replace(names,"Over 18", "over 18w")) %>%
    mutate(names = str_replace(names,"Wane", "Second Dose")) %>%
    mutate(names = str_replace(names,'Second Dose >14d', 'Second Dose 2-10w'))
  
  order <- results_to_plot$names |> unique() |> sort()
  order <- (c('Not Vaccinated',order[order!='Not Vaccinated']))
  
  results_to_plot$names   <- factor(results_to_plot$names, levels=order)
  results_to_plot$variant <- factor(results_to_plot$variant, levels=c('all','WT','Alpha','Delta','Omicron'))
  
  x_label <- ''
  if(filter_by=='restriction')
  {
    results_to_plot <- results_to_plot %>% mutate(names = '')
    x_label         <- 'Restriction Level'
  }
  
  tbl_vac <- results_to_plot %>% arrange(variant,names) %>%
    mutate(variant = case_when(variant=='all' ~ 'Full time period', TRUE ~ variant)) %>%
    mutate(names = str_replace(names,'<','$<$'),
           names = str_replace(names,'>','$>$')) %>%
    dplyr::select(names,coefficients,std_error,case_definition,variant) %>%
    mutate(case_definition = case_when(case_definition=='hosp'~'Hospitalisations',
                                       case_definition=='death'~'Deaths',
                                       case_definition=='pillar2pcr'~ 'Pillar 2 PCR Cases')) %>%
    mutate(coef        = exp(coefficients),
           lower_bound = exp(coefficients-1.96*std_error),
           upper_bound = exp(coefficients+1.96*std_error)) %>%
    dplyr::select(-c(coefficients,std_error)) %>%
    mutate(vacc_eff    = 1 - coef,
           vacc_eff_lb = 1 - upper_bound,
           vacc_eff_up = 1 - lower_bound ) %>%
    mutate(IRR = paste0(round(coef,3), " (",round(lower_bound,3),'-',round(upper_bound,3),')'),
           'Vaccine effectiveness' = paste0(round(vacc_eff,3), " (",round(vacc_eff_lb,3),'-',round(vacc_eff_up,3),')')) %>%
    dplyr::select(-c(coef,upper_bound, lower_bound,vacc_eff,vacc_eff_lb,vacc_eff_up)) %>%
    pivot_wider(names_from = case_definition, values_from = c(IRR,'Vaccine effectiveness'),values_fill = '..', names_glue = "{.value} ({case_definition})") %>%
    mutate(names=as.character(names)) %>% rename(Vaccine=names)
  
  tbl_vac_csv <- insert_blank_rows_variant(tbl_vac ,'variant')
  
  return(tbl_vac_csv)
}

ethnicity_table <- function(fits,
                            min_p_value = 1 )
{
  filter_by = 'ethnicity_simple'
  include = c('death','hosp','pillar2pcr')
  
  number_of_fits <- length(fits)
  names_of_fits  <- names(fits)
  
  results_to_plot <- tibble(NULL)
  
  for(i in 1:number_of_fits)
  {
    details <- str_split(names_of_fits[i],'_')[[1]]
    
    if(details[1]=='level'|details[1]=='national')
    {
      details[1] = paste(details[1],details[2],sep="_")
      details    = details[-2]
    }
    
    if(details[1]=='') variant = 'all' else variant = details[1]
    case_definition = details[3]
    
    fit <- fits[[i]]@model$coefficients_table %>% mutate(exp_coef = exp(coefficients)) %>%
      filter(str_detect(names, filter_by)) %>% filter(p_value < min_p_value|is.na(p_value)) # need to include NA as this is for the 
    # reference values
    if(dim(fit)[1]==0) next
    
    # need to add all the reference groups back in!
    if(str_detect('ethnicity_simple',filter_by))
    {
      tmp <- as_tibble(list(names='ethnicity_simple.white',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    fit$variant <- variant
    fit$case_definition <- case_definition
    
    results_to_plot <- bind_rows(results_to_plot,fit)
  }
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  results_to_plot$names <- str_replace(results_to_plot$names,'ethnicity_simple.','')
  
  results_to_plot <- results_to_plot %>% mutate(names = case_when(names=='white' ~ 'White',
                                                                  names=='south_asian' ~ 'South Asian',
                                                                  names=='other_asian' ~ 'Asian (other)',
                                                                  names=='black' ~ 'Black',
                                                                  names=='mixed_other' ~ 'Mixed/Other',
                                                                  names=='mixed' ~ 'Mixed',
                                                                  names=='other' ~ 'Other',
                                                                  TRUE ~ names))
  
  results_to_plot$names <- factor(results_to_plot$names, levels=c("White", 'South Asian', 'Asian (other)','Black',
                                                                  'Mixed/Other', 'Mixed', 'Other', "unknown", 
                                                                  "south_asian","other_asian",
                                                                  "asian", "black", "mixed_other"))
  
  results_to_plot <- results_to_plot %>% 
    filter(!is.na(exp_coef) & !is.na(names) & case_definition %in% include )
  
  x_label <- ''
  
  ethnicity_table <- results_to_plot %>% 
    mutate(variant = case_when(variant=='all' ~ 'Full time period', TRUE ~ variant)) %>%
    dplyr::select(names,coefficients,std_error,case_definition,variant) %>%
    mutate(case_definition = case_when(case_definition=='hosp'~'Hospitalisations',
                                       case_definition=='death'~'Deaths',
                                       case_definition=='pillar2pcr'~ 'Pillar 2 PCR Cases')) %>%
    mutate(coef        = exp(coefficients),
           lower_bound = exp(coefficients-1.96*std_error),
           upper_bound = exp(coefficients+1.96*std_error)) %>%
    mutate(IRR = paste0(round(coef,3), " (",round(lower_bound,3),'-',round(upper_bound,3),')')) %>%
    dplyr::select(-c(coefficients, std_error,coef,upper_bound, lower_bound)) %>%
    pivot_wider(names_from = case_definition, values_from = c(IRR),values_fill = '..', names_glue = "{.value} ({case_definition})") %>%
    mutate(names=as.character(names)) %>% rename(Ethnicity=names)
  
  ethnicity_table_csv <- insert_blank_rows_variant(ethnicity_table ,'variant')
  
  return(ethnicity_table_csv)
}

IMD_table <- function(fits,
                      min_p_value = 1 )
{
  filter_by = 'IMD'
  include = c('death','hosp','pillar2pcr')
  
  number_of_fits <- length(fits)
  names_of_fits  <- names(fits)
  
  results_to_plot <- tibble(NULL)
  
  for(i in 1:number_of_fits)
  {
    details <- str_split(names_of_fits[i],'_')[[1]]
    
    if(details[1]=='level'|details[1]=='national')
    {
      details[1] = paste(details[1],details[2],sep="_")
      details    = details[-2]
    }
    
    if(details[1]=='') variant = 'all' else variant = details[1]
    case_definition = details[3]
    
    fit <- fits[[i]]@model$coefficients_table %>% mutate(exp_coef = exp(coefficients)) %>%
      filter(str_detect(names, filter_by)) %>% filter(p_value < min_p_value|is.na(p_value)) # need to include NA as this is for the 
    # reference values
    if(dim(fit)[1]==0) next
    
    # need to add all the reference groups back in!
    if(str_detect('IMD',filter_by))
    {
      tmp <- as_tibble(list(names='IMD_national_quintile.IMD1',coefficients=0,std_error=0,z_value=10,p_value=0,standardized_coefficients=0,exp_coef=1))
      fit <- bind_rows(tmp,fit)  
    }
    
    fit$variant <- variant
    fit$case_definition <- case_definition
    
    results_to_plot <- bind_rows(results_to_plot,fit)
  }
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  results_to_plot$names <- str_replace(results_to_plot$names,'IMD_national_quintile.','')
  
  results_to_plot <- results_to_plot %>% 
    filter(!is.na(exp_coef) & !is.na(names) & case_definition %in% include )
  
  x_label <- ''
  
  IMD_table <- results_to_plot %>% 
    mutate(variant = case_when(variant=='all' ~ 'Full time period', TRUE ~ variant)) %>%
    dplyr::select(names,coefficients,std_error,case_definition,variant) %>%
    mutate(case_definition = case_when(case_definition=='hosp'~'Hospitalisations',
                                       case_definition=='death'~'Deaths',
                                       case_definition=='pillar2pcr'~ 'Pillar 2 PCR Cases')) %>%
    mutate(coef        = exp(coefficients),
           lower_bound = exp(coefficients-1.96*std_error),
           upper_bound = exp(coefficients+1.96*std_error)) %>%
    mutate(IRR = paste0(round(coef,3), " (",round(lower_bound,3),'-',round(upper_bound,3),')')) %>%
    dplyr::select(-c(coefficients, std_error,coef,upper_bound, lower_bound)) %>%
    pivot_wider(names_from = case_definition, values_from = c(IRR),values_fill = '..', names_glue = "{.value} ({case_definition})") %>%
    mutate(names=as.character(names)) %>% rename(Deprivation=names)
  
  IMD_table_csv <- insert_blank_rows_variant(IMD_table ,'variant')
  
  return(IMD_table_csv)
}
