library(readxl)
library(webchem)
library(stringr)
library(ggplot2)
library(httr)
library(jsonlite)
library(stringdist)
library(latex2exp) # plot formatting

CENSUS_API = "https://geo.fcc.gov/api/census/area?lat=%f&lon=%f&censusYear=2020&format=json"

make_simple = function(input_strs){
  
  if(length(input_strs) == 1){
    input_strs = input_strs[[1]]
    
    if(length(input_strs) == 1 && is.na(input_strs)){
      return(NA)
    }
  }
  
  # lower case
  input_strs = tolower(input_strs)
  # remove spaces
  input_strs = gsub(" ","",input_strs)
  
  # remove commas, dashes, parentheses, brackets, equals
  input_strs = gsub("\\=|\\:|\\-|\\,","",input_strs)
  input_strs = gsub("\\(|\\)","",input_strs)
  input_strs = gsub("\\[|\\]","",input_strs)
  
}
# ------- PFC Synonyms --------- #
get_syns = function(){
  # Make tables for each contaminant
  PUBCHEM_SYNS = rbind(pc_synonyms("11CL-PF3OUDS"),
                       pc_synonyms("9CL-PF3ONS"),
                       pc_synonyms("ADONA"),
                       pc_synonyms("HFPO-DA"),
                       pc_synonyms("NETFOSAA"),
                       pc_synonyms("NMEFOSAA"),
                       pc_synonyms("PFBS"),
                       pc_synonyms("PFDA"),
                       pc_synonyms("PFDOA"),
                       pc_synonyms("PFHPA"),
                       pc_synonyms("PFHXA"),
                       pc_synonyms("PFHXS"),
                       pc_synonyms("PFNA"),
                       pc_synonyms("PFOA"),
                       c("Total PFOA and PFOS"), # pc_synonyms("PFOA+PFOS"),
                       pc_synonyms("PFOS"),
                       pc_synonyms("PFTA"),
                       pc_synonyms("PFTRDA"),
                       pc_synonyms("PFUNA" ),
                       c("PFAS Sum") # pc_synonyms("Total PFAS")
  )
  colnames(PUBCHEM_SYNS) = "PUBCHEM_SYNS"
  SIMPLE_SYNS = lapply(PUBCHEM_SYNS, make_simple)
  
  syns_df = data.frame(cbind(PUBCHEM_SYNS=PUBCHEM_SYNS,
                             SIMPLE_SYNS=SIMPLE_SYNS),
                       row.names=sort(colnames(contaminant_df)))
  return(syns_df)
}


# ------- Drinking Water --------- #
get_drinking_water = function(){
  drinking_water_df = read_excel("data/drinking_water.xlsx")
  drinking_water_df$Concentration = as.numeric(drinking_water_df$Concentration)
  drinking_water_df = drinking_water_df[!is.na(drinking_water_df$Concentration),]
  
  drinking_water_df = drinking_water_df[drinking_water_df$Contaminant == 
                                          "Total PFAS",]
  drinking_water_df$County = str_replace(drinking_water_df$County, 
                                         ", CO", "")
  drinking_water_df$County = str_replace(drinking_water_df$County, 
                                         "City and County of ", "")
  drinking_water_df$County = str_replace(drinking_water_df$County, 
                                         " County", "")
  drinking_water_df$County = str_to_upper(drinking_water_df$County)
  
  return(drinking_water_df)
}


# ------- Industries of Interest --------- #
get_ioi = function(counties){
  ioi_df = read_excel("data/industries_of_interest.xlsx")
  ioi_df = ioi_df[ioi_df$Status == "Active", ]
  ioi_df = ioi_df[which(ioi_df$FAC_COUNTY %in% counties), ]
  return(ioi_df)
}


# ------- Discharge Sites --------- #
get_discharge = function(){
  discharge_df = read_excel("data/discharge_monitoring.xlsx")
  discharge_df$Year = as.numeric(discharge_df$Year)
  discharge_df = discharge_df[discharge_df$Year < 2021, ]
  
  discharge_facilities = unique(discharge_df$Facility)
  
  zips = substring(unique(discharge_df$`ZIP Code`),1,5)
  sub_discharge_df = data.frame(PERMIT = character(length(discharge_facilities)),
                                FACILITY = discharge_facilities,
                                ZIP = numeric(length(discharge_facilities)),
                                COUNTY = character(length(discharge_facilities)))
  
  for(i in 1:length(discharge_facilities)){
    facility = discharge_facilities[i]
    lat = as.numeric(unique(discharge_df[discharge_df$Facility == facility, 
                                         "Latitude"]))
    lon = as.numeric(unique(discharge_df[discharge_df$Facility == facility, 
                                         "Longitude"]))
    zip = as.numeric(substring(unique(discharge_df[discharge_df$Facility == facility, 
                                                   "ZIP Code"]), 1, 5))
    
    if(length(lat) > 1){
      message(sprintf("Error for locations for facility %d", i))
    }
    else{
      get_str = sprintf(CENSUS_API, lat, lon)
      res = GET(get_str)
      data = fromJSON(rawToChar(res$content))
      county = str_to_upper(str_replace(unique(data$results$county_name), 
                                        " County", ""))
    }
    
    permit = unique(discharge_df[discharge_df$Facility == facility, 
                                 "NPDES Permit Number"]$`NPDES Permit Number`)
    total = sum(discharge_df[discharge_df$Facility == facility, 
                             "Average Concentration (mg/L)"], na.rm=T)
    sub_discharge_df[sub_discharge_df$FACILITY == facility, 
                     "PERMIT"] = permit
    
    if(length(county) > 1){
      message(sprintf("Error for number of counties for facility %d", i))
      num_fill = length(county)
      sub_discharge_df[seq(i + num_fill, nrow(sub_discharge_df) + num_fill - 1),] = sub_discharge_df[seq(i + num_fill - 1, nrow(sub_discharge_df)),]
      
      sub_discharge_df[[i, "PERMIT"]] = permit
      sub_discharge_df[[i, "ZIP"]] = zip
      sub_discharge_df[[i, "COUNTY"]] = county[1]
      sub_discharge_df[[i, "TOTAL_AVG_CONC_mgL"]] = total
      for(j in 1:length(county) - 1){
        sub_discharge_df[[i + j, "PERMIT"]] = permit
        sub_discharge_df[[i + j, "FACILITY"]] = facility
        sub_discharge_df[[i + j, "ZIP"]] = zip
        sub_discharge_df[[i + j, "COUNTY"]] = county[j + 1]
        sub_discharge_df[[i + j, "TOTAL_AVG_CONC_mgL"]] = total
      }
      
    }
    else{
      sub_discharge_df[sub_discharge_df$FACILITY == facility, 
                       "ZIP"] = zip
      sub_discharge_df[sub_discharge_df$FACILITY == facility, 
                       "COUNTY"] = county
      sub_discharge_df[sub_discharge_df$FACILITY == facility, 
                       "TOTAL_AVG_CONC_mgL"] = total
    }
  }
  return(sub_discharge_df)
}
# ------- Federal Sites --------- #
get_fed_sites = function(){
  fed_sites_df = read_excel("data/federal_sites.xlsx")
  
  for(i in 1:nrow(fed_sites_df)){
    
    lat = as.numeric(unique(fed_sites_df[[i, "Latitude"]]))
    lon = as.numeric(unique(fed_sites_df[[i, "Longitude"]]))
    
    if(length(lat) > 1){
      message(sprintf("Error for locations for federal site %d", i))
    }
    else{
      get_str = sprintf(CENSUS_API, lat, lon)
      res = GET(get_str)
      data = fromJSON(rawToChar(res$content))
      county = str_to_upper(str_replace(unique(data$results$county_name), 
                                        " County", ""))
      fed_sites_df[[i, "County"]] = county
    }
  }
}


# ------- Land Area Statistics per County --------- #
get_area = function(){
  area_df = read.csv("data/State of Colorado Counties - 2020 Census - Data as of January 1, 2020.csv")
  area_df$BASENAME = str_to_upper(area_df$BASENAME)
  return(area_df)
}


# ------- Biosolids Permits --------- #
get_biosolids = function(){
  active_permits_df = read_excel("data/Active permits 3-1-23.xlsx")
  
  i_keep = which(rowSums(is.na(active_permits_df)) != ncol(active_permits_df))
  active_permits_df = active_permits_df[i_keep, ]
  
  i_permit = active_permits_df$GeneralPermitType == "COBMP-Biosolids user"
  biosolids_df = active_permits_df[i_permit,]
  
  i_na = which(rowSums(is.na(biosolids_df)) != ncol(biosolids_df))
  biosolids_df = biosolids_df[i_na, ]
  
  # Refactor date column
  biosolids_df$EffectiveDate = as.Date(biosolids_df$EffectiveDate, 
                                       format="%Y/%m/%d")
  
  # Filter by end date of sampling project
  end_date = as.Date("2020-04-11")
  biosolids_df = biosolids_df[biosolids_df$EffectiveDate < end_date, ]
  biosolids_df = biosolids_df[!is.na(biosolids_df$FacilityCounty), ]
  
  # Clean up the counties to match main dataframe
  biosolids_df = biosolids_df[which(biosolids_df$FacilityCounty %in% 
                                      c("Milwaukee", "West") == FALSE), ]
  biosolids_df$COUNTY = str_to_upper(biosolids_df$FacilityCounty)
  biosolids_df[biosolids_df$COUNTY == "ARCHLETA", "COUNTY"] = "ARCHULETA"
  biosolids_df[biosolids_df$COUNTY == "CONEJOE", "COUNTY"] = "CONEJOS"
  biosolids_df[biosolids_df$COUNTY == "LA PLATE", "COUNTY"] = "LA PLATA"
  biosolids_df[biosolids_df$COUNTY == "RIO GRAND", "COUNTY"] = "RIO GRANDE"
  biosolids_df$COUNTY = str_replace(biosolids_df$COUNTY, " COUNTY", "")
  biosolids_df$COUNTY = str_replace(biosolids_df$COUNTY, "`", "")
  biosolids_df$COUNTY = str_replace(biosolids_df$COUNTY, " AND ", " & ")
  
  # Add an ALT_COUNTY option to count later
  biosolids_df[c("COUNTY","ALT_COUNTY")] = str_split_fixed(biosolids_df$COUNTY, 
                                                           "\\ & ", 2)
  return(biosolids_df)
}


# ------- Final Covariate Data Frame --------- #
get_pfc = function(){
  
  drinking_water_df = get_drinking_water()
  num_samples = length(drinking_water_df$Contaminant)
  counties = unique(drinking_water_df$County)
  
  
  ioi_df = get_ioi(counties)
  area_df = get_area()
  sub_discharge_df = get_discharge()
  fed_sites_df = get_fed_sites()
  biosolids_df = get_biosolids()
  
  num_ioi = numeric(num_samples)
  areas = numeric(num_samples)
  
  fed_sites_indicator = logical(num_samples)
  num_fed_sites = numeric(num_samples)
  
  discharge_indicator = logical(num_samples)
  num_discharge = numeric(num_samples)
  discharge_conc = numeric(num_samples)
  
  biosolids_indicator = logical(num_samples)
  num_biosolids = numeric(num_samples)
  
  final_counties = character(num_samples)
  concentration = numeric(num_samples)
  
  for(i in 1:num_samples){
    concentration[i] = as.numeric(drinking_water_df[[i, "Concentration (ng/L)"]])
    county = drinking_water_df[[i, "County"]]
    final_counties[i] = str_replace(str_to_upper(county), " ", "_")
    num_ioi[i] = sum(ioi_df$FAC_COUNTY == county)
    areas[i] = area_df[area_df$BASENAME == county, "AREALAND"]/1e6
    
    discharge_indicator[i] = county %in% sub_discharge_df$COUNTY
    if(discharge_indicator[i]  == 1){
      num_discharge[i]  = sum(sub_discharge_df$COUNTY == county)
      discharge_conc[i] = mean(sub_discharge_df[sub_discharge_df$COUNTY == county, 
                                                "TOTAL_AVG_CONC_mgL"])
    }
    
    fed_sites_indicator[i] = county %in% fed_sites_df$County
    if(fed_sites_indicator[i]  == 1){
      num_fed_sites[i]  = sum(fed_sites_df$County == county)
    }
    
    biosolids_indicator[i] = county %in% rbind(biosolids_df$COUNTY, 
                                               biosolids_df$ALT_COUNTY)
    
    if(biosolids_indicator[i] == 1){
      num_biosolids[i] = sum(biosolids_df$COUNTY == county) + 
        sum(biosolids_df$ALT_COUNTY == county)
    }
  }
  total_pfas_df = data.frame(CONCENTRATION_ngL = concentration,
                        COUNTY = final_counties,
                        IOI_COUNT = num_ioi,
                        AREA_KM2 = areas,
                        NUM_IOI_PER_KM2 = num_ioi/areas,
                        BOOL_DISCHARGE = discharge_indicator,
                        NUM_DISCHARGE = num_discharge,
                        AVG_AVG_CONC_mgL = discharge_conc,
                        BOOL_FED_SITES = fed_sites_indicator,
                        NUM_FED_SITES = num_fed_sites,
                        BOOL_BIOSOLIDS = biosolids_indicator,
                        NUM_BIOSOLIDS = num_biosolids,
                        NUM_BIOSOLIDS_PER_KM2 = num_biosolids/areas)
  
  return(total_pfas_df)
}




# --------------- Compare Covariates to Data -------------------------- #
plot_pfcs_log = function(total_pfas_df){
  ggplot(total_pfas_df, aes(x=NUM_IOI_PER_KM2, y=log10(CONCENTRATION_ngL), 
                            shape = BOOL_FED_SITES,
                            color = BOOL_DISCHARGE)) + 
    geom_point(alpha = 0.7) + 
    scale_x_log10() + 
    scale_shape_discrete(name="Federal Sites") + 
    scale_color_discrete(name="Discharge Permitted Sites") + 
    labs(x = TeX("Number of Industries of Interest $\\textit{n/km}^2$"), 
         y = TeX("$\\log(Concentration) $$(\\textit{ng/l})$$")) + 
    theme(text = element_text(family = "Times New Roman"))
  ggsave("ioi_fed_discharge.png")
  
  ggplot(total_pfas_df, aes(x=NUM_BIOSOLIDS_PER_KM2, y=log10(CONCENTRATION_ngL), 
                            shape = BOOL_FED_SITES,
                            color = BOOL_DISCHARGE)) + 
    geom_point(alpha = 0.7) + 
    scale_x_log10() + 
    scale_shape_discrete(name="Federal Sites") + 
    scale_color_discrete(name="Discharge Permitted Sites") + 
    labs(x = TeX("Number of Biosolids Sites $\\textit{n/km}^2$"), 
         y = TeX("$\\log(Concentration) $$(\\textit{ng/l})$$")) + 
    theme(text = element_text(family = "Times New Roman"))
  
  ggsave("biosolids_fed_discharge.png")
  
  ggplot(total_pfas_df, aes(x=NUM_IOI_PER_KM2, y=log10(CONCENTRATION_ngL), 
                            color = COUNTY,
                            shape = BOOL_BIOSOLIDS)) + 
    geom_point(alpha = 0.7) + 
    scale_x_log10() + 
    labs(x = TeX("Number of Industries of Interest $\\textit{n/km}^2$"), 
         y = TeX("Concentration $$(\\textit{ng/l})$$")) + 
    theme(text = element_text(family = "Times New Roman"))
}

plot_pfcs = function(drinking_water_df){
  contaminant_names = unique(drinking_water_df$Contaminant)
  contaminant_counts = numeric(length(contaminant_names))
  colors = viridis::viridis_pal(option="D", 
                                begin=0.8, 
                                end=0.2, alpha=0.8)(length(contaminant_names))
  
  for(i in 1:length(contaminant_names)){
    
    vals = drinking_water_df[drinking_water_df$Contaminant == contaminant_names[i], 
                             "Concentration"]
    vals$y = median(vals$Concentration)
    title = paste(contaminant_names[i], "Density")
    p = ggplot(vals, aes(x=Concentration)) + 
      geom_density(fill=colors[i], alpha=0.6) +
      scale_x_log10() + 
      geom_vline(data=vals, aes(xintercept=y, color=colors[i]),
                 linetype="dashed", show.legend = FALSE) + 
      # geom_text(data=vals, aes(label=median(Concentration), 
      #                          x=median(Concentration) + 0.3,
      #                          y=max(density(Concentration)$y) + 0.01)) + 
      labs(x="Concentration (ppt)", y="Density", title=title) + 
      theme_bw() 
    
    
    file = paste0(contaminant_names[i], "_density.pdf")
    ggsave(file,p)
  }
}

get_acs = function(counties){
  
  unique_counties = unique(counties)
  covar_idx = c(2:4,6:18,67:72,75,80)
  acs_df = read.csv("data/ACSDP5Y2020.DP05-2023-04-30T233324.csv")
  acs_df$Label..Grouping. = str_trim(acs_df$Label..Grouping.)
  i_percent = which(!is.na(rowSums(str_locate(colnames(acs_df),
                                              "..Percent"))) & 
                      is.na(rowSums(str_locate(colnames(acs_df),
                                               "..Percent.Margin.of.Error"))) |
                      !is.na(rowSums(str_locate(colnames(acs_df), 
                                                "..Estimate"))))
  keep_names = colnames(acs_df)[i_percent]
  keep_names = c("Label..Grouping.", keep_names)
  acs_df = acs_df[covar_idx, keep_names]
  colnames(acs_df) = str_replace(colnames(acs_df), 
                                 ".County..Colorado", "")
  colnames(acs_df) = str_replace(colnames(acs_df), 
                                 "..Percent", "_Percent")
  colnames(acs_df) = str_replace(colnames(acs_df), 
                                 "..Estimate", "_Estimate")
  colnames(acs_df) = str_replace(colnames(acs_df), 
                                 "..Grouping.", "_Grouping")
  colnames(acs_df) = str_replace_all(colnames(acs_df), 
                                     "[.]", "_")
  colnames(acs_df) = str_to_upper(colnames(acs_df))
  
  pop_df = acs_df[which(acs_df$LABEL_GROUPING == "Total population"),
                  c("LABEL_GROUPING", paste0(unique_counties, "_ESTIMATE"))]
  colnames(pop_df) = str_replace(colnames(pop_df), "_ESTIMATE", "")
  
  per_df = acs_df[which(acs_df$LABEL_GROUPING != "Total population"), 
                  c("LABEL_GROUPING", paste0(unique_counties, "_PERCENT"))]
  colnames(per_df) = str_replace(colnames(per_df), "_PERCENT", "")
  
  
  final_df = rbind(pop_df, per_df)
  final_df = final_df[,c("LABEL_GROUPING",counties)]
  final_df = data.frame(t(final_df))
  names(final_df) = final_df["LABEL_GROUPING",]
  final_df = final_df[c(2:nrow(final_df)),]
  
  for (i in 1:ncol(final_df)){
    if(i == 1){
      final_df[, i] = as.numeric(sub(",","",final_df[,i]))
    }else{
      final_df[, i] = as.numeric(sub("%","",final_df[,i]))
    }
  }
  
  final_df[, "35 years and over"] = rowSums(final_df[,c("35 to 44 years",
                                                        "45 to 54 years",
                                                        "55 to 59 years",
                                                        "60 to 64 years",
                                                        "65 to 74 years",
                                                        "75 to 84 years",
                                                        "85 years and over")])
  return(final_df)
}


total_pfas_df = get_pfc()
acs_df = get_acs(total_pfas_df$COUNTY)
total_pfas_df = cbind(total_pfas_df, acs_df)

demos = c("Female",
          "35 years and over",
          "Black or African American",
          "American Indian and Alaska Native",
          "Asian",
          "Native Hawaiian and Other Pacific Islander",
          "Some other race",
          "Hispanic or Latino (of any race)")

pfas_vars = c("AREA_KM2", "IOI_COUNT", "NUM_IOI_PER_KM2", 
              "BOOL_DISCHARGE", "NUM_DISCHARGE", "AVG_AVG_CONC_mgL" ,                         
              "BOOL_FED_SITES", "NUM_FED_SITES",                             
              "BOOL_BIOSOLIDS", "NUM_BIOSOLIDS","NUM_BIOSOLIDS_PER_KM2")

counties = unique(total_pfas_df$COUNTY)

for (i in 1:length(counties)){
  cty = (total_pfas_df$COUNTY == counties[i]) + 0
  total_pfas_df[,counties[i]] = cty
}


ggplot(total_pfas_df, aes(x=`Black or African American`, y=log10(CONCENTRATION_ngL), 
                          shape = BOOL_FED_SITES,
                          color = BOOL_DISCHARGE)) + 
  geom_point(alpha = 0.7) + 
  # scale_x_log10() + 
  scale_shape_discrete(name="Federal Sites") + 
  scale_color_discrete(name="Discharge Permitted Sites") + 
  labs(x = TeX("Number of Industries of Interest $\\textit{n/km}^2$"), 
       y = TeX("$\\log(Concentration) $$(\\textit{ng/l})$$")) + 
  theme(text = element_text(family = "Times New Roman"))


ggplot(total_pfas_df) + 
  geom_density(aes(`Black or African American`))

ggplot(total_pfas_df) + 
  geom_density(aes(`American Indian and Alaska Native`))
                                
ggplot(covid_20210307) + 
  geom_density(aes(hs))

ggplot(total_pfas_df, aes(x=NUM_IOI_PER_KM2, y=log10(CONCENTRATION_ngL), 
                          shape = BOOL_FED_SITES,
                          color = BOOL_DISCHARGE)) + 
  geom_point(alpha = 0.7, size=3) + 
  scale_shape_discrete(name="Federal Sites") + 
  scale_color_discrete(name="Discharge Permitted Sites") + 
  labs(x = TeX("Number of Industries of Interest $\\textit{n/km}^2$"), 
       y = TeX("$\\log(Concentration) $$(\\textit{ng/l})$$")) + 
  theme(text = element_text(family = "Times New Roman"))
ggsave("ioi_fed_discharge.png")

ggplot(total_pfas_df, aes(x=NUM_BIOSOLIDS_PER_KM2, y=log10(CONCENTRATION_ngL), 
                          shape = BOOL_FED_SITES,
                          color = BOOL_DISCHARGE)) + 
  geom_point(alpha = 0.7, size=3) + 
  scale_shape_discrete(name="Federal Sites") + 
  scale_color_discrete(name="Discharge Permitted Sites") + 
  labs(x = TeX("Number of Biosolids Sites $\\textit{n/km}^2$"), 
       y = TeX("$\\log(Concentration) $$(\\textit{ng/l})$$")) + 
  theme(text = element_text(family = "Times New Roman"))

ggsave("biosolids_fed_discharge.png")


contaminant_names = unique(drinking_water_df$Contaminant)
contaminant_counts = numeric(length(contaminant_names))
colors = viridis::viridis_pal(option="D", 
                              begin=0.8, 
                              end=0.2, alpha=0.8)(10)


vals = drinking_water_df[drinking_water_df$Contaminant == contaminant_names[i], 
                         "Concentration"]
vals$y = median(vals$Concentration)
title = paste("Total PFAS", "Density")
ggplot(total_pfas_df, aes(x=CONCENTRATION_ngL)) + 
  geom_density(fill=colors[1], alpha=0.6) +
  # scale_x_log10() + 
  # geom_vline(data=vals, aes(xintercept=y, color=colors[i]),
  #            linetype="dashed", show.legend = FALSE) + 
  # geom_text(data=vals, aes(label=median(Concentration), 
  #                          x=median(Concentration) + 0.3,
  #                          y=max(density(Concentration)$y) + 0.01)) + 
  labs(x="Concentration (ppt)", y="Density", title=title) + 
  theme_bw() 

ggsave("total_pfas_density.png")
  
