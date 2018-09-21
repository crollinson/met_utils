# Function to calculate growing degree days for a given date based on a reference date
# Default reference date will be Jan 1 of the year provided.

calc.gdd <- function(path.met, gdd.date, ref.date=NULL, base.temp=5){
  # path.met = CF-standard meteorology in netcdf files; 1 file per year; preferably daily resolution
  #            This function is initially written to work with the format proudced by extract_Daymet.R
  # gdd.date = date for which you want the cumulative GDD (e.g. date of leaf out)
  # ref.date = reference date for which you want to start accumulating GDD; 
  #            if left blank, will defulat to Jan 1 of the year for which you're caluclating GDD
  # base.temp = the base temperature for GDD calculation; defaults to 5 C
  
  
  
}