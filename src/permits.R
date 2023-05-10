# Then for each facility in the discharge list, get the discharge flow rate 
# based on the permit number
discharge_permits = unique(discharge_df$`NPDES Permit Number`)
permit_dir = "data/discharge_permits/"
start_date = as.Date("2017-01-01")
end_date = as.Date("2020-04-11")

read_permits = function(index, discharge_permits, discharge_df, permit_dir="data/discharge_permits/", 
                        start_date, end_date){
  permit_file = paste0(permit_dir,paste0(discharge_permits[index],".csv"))
  permit_df = read.csv(permit_file, skip=3)
  
  permit_idx = discharge_df$`NPDES Permit Number` == discharge_permits[index]
  # permit_params = unique(discharge_df[permit_idx, "Parameter Description"]$`Parameter Description`)
  permit_codes = unique(discharge_df[permit_idx, "Parameter Code"]$`Parameter Code`)
  
  permit_df = permit_df[which(permit_df$Parameter.Code %in% permit_codes),]
  permit_df = permit_df[permit_df$Limit.Type == "30DA AVG",]
  
  permit_df$Monitoring.Period.Date = as.Date(permit_df$Monitoring.Period.Date, 
                                             format="%m/%d/%Y")
  
  permit_df = permit_df[permit_df$Monitoring.Period.Date > start_date, ]
  permit_df = permit_df[permit_df$Monitoring.Period.Date < end_date, ]
  
  # Just do the most recent monitoring date
  monitoring_dates = range(permit_df$Monitoring.Period.Date)
  last_monitor_date = monitoring_dates[length(monitoring_dates)]
  permit_df = permit_df[permit_df$Monitoring.Period.Date == last_monitor_date, ]
  facility = unique(discharge_df[discharge_df$`NPDES Permit Number` == discharge_permits[i], "Facility"]$Facility)
  
  units = unique(permit_df$DMR.Value.Unit)
  
  if(length(units) > 1){
    message(sprintf("More than one unit for permit %d", i))
    
    for(j in 1:length(units)){
      if(units[j] == "mg/L"){
        permit_df[permit_df$DMR.Value.Unit == "mg/L", "DMR.Value"] = 1000 * as.numeric(permit_df[permit_df$DMR.Value.Unit == "mg/L", "DMR.Value"])
        permit_df[permit_df$DMR.Value.Unit == "mg/L", "DMR.Value.Unit"] = "ug/L"
      }
      else if(units[j] == "ng/L"){
        permit_df[permit_df$DMR.Value.Unit == "ng/L", "DMR.Value"] = 0.001 * as.numeric(permit_df[permit_df$DMR.Value.Unit == "ng/L", "DMR.Value"])
        permit_df[permit_df$DMR.Value.Unit == "ng/L", "DMR.Value.Unit"] = "ug/L"
      }
      else{
        message("Error, confused")
      }
      
    }
  }
}

i = 4
permit_df = read_permits(index=i, discharge_permits=discharge_permits, 
                         discharge_df = discharge_df, permit_dir = permit_dir,
                         start_date = start_date, end_date = end_date)





permit = unique(discharge_df[discharge_df$Facility == discharge_facilities[1], "NPDES Permit Number"]$`NPDES Permit Number`)

num_2017_2020 = sum(as.numeric(permit_df$DMR.Value) > 0, na.rm=T) # 30 :(
num_2017_2023 = sum(as.numeric(permit_df$DMR.Value) > 0, na.rm=T) # 1729