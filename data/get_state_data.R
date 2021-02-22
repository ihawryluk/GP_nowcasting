library(openxlsx)
library(lubridate)
library(ggplot2)
library(tidyr)
library(repr)
library(ggpubr)
library(readr)
library(tidyverse)
options(repr.plot.width=12, repr.plot.height=12)

#################################################
# before running the script update the paths in lines 15-16, 40 and 50!!!!
#################################################

path <- paste0("C:/Users/name/Downloads/INFLUD-ALL/")
path2 <- paste0("C:/Users/name/Desktop/nowcasting/data/")

RELEASE_DATE=c(
  "07-07-2020","14-07-2020","21-07-2020",
  "29-07-2020","03-08-2020","10-08-2020",
  "17-08-2020","24-08-2020","31-08-2020",
  "07-09-2020","14-09-2020","21-09-2020",
  "28-09-2020","05-10-2020","12-10-2020",
  "19-10-2020","26-10-2020","02-11-2020",
  "10-11-2020","16-11-2020","23-11-2020",
  "30-11-2020","07-12-2020","14-12-2020", 
  "21-12-2020", "28-12-2020", 
  "04-01-2021", "11-01-2021","18-01-2021",
  "25-01-2021", "01-02-2021", "08-02-2021"
)

####################################################
# load all sivep releases

# 2020 releases
for (RELEASE_DATE_i in RELEASE_DATE){
  print(RELEASE_DATE_i)
  url = paste("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2020/INFLUD-", RELEASE_DATE_i, sep='')
  url = paste(url, ".csv", sep='')
  dest = paste("C:/Users/name/Downloads/INFLUD-ALL/INFLUD-", RELEASE_DATE_i, sep='')
  dest = paste(dest, ".csv", sep='')
  download.file(url, dest)
}

# 2021 releases
for (RELEASE_DATE_i in RELEASE_DATE[-1]){
  print(RELEASE_DATE_i)
  url = paste("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2021/INFLUD21-", RELEASE_DATE_i, sep='')
  url = paste(url, ".csv", sep='')
  dest = paste("C:/Users/name/Downloads/INFLUD-ALL/INFLUD21-", RELEASE_DATE_i, sep='')
  dest = paste(dest, ".csv", sep='')
  download.file(url, dest)
}
#####################################################

df_SIVEP_nowcast = NULL
df_SIVEP_padded_nowcast = NULL
for (RELEASE_DATE_i in RELEASE_DATE){
  print(RELEASE_DATE_i)
  
  if ( dmy(RELEASE_DATE_i) > dmy("10-01-2021") )
  {
    df_SIVEP_ORIGINAL20=read_csv2(paste0(path,"INFLUD-",RELEASE_DATE_i,".csv"))
    df_SIVEP_ORIGINAL21=read_csv2(paste0(path,"INFLUD21-",RELEASE_DATE_i,".csv"))  
    df_SIVEP_ORIGINAL <- rbind(df_SIVEP_ORIGINAL20,df_SIVEP_ORIGINAL21)
  }
  else if ( year(dmy(RELEASE_DATE_i)) < dmy("10-01-2021") )
  {
    df_SIVEP_ORIGINAL=read_csv2(paste0(path,"INFLUD-",RELEASE_DATE_i,".csv"))   
  }
  
  df_SIVEP_ORIGINAL->df_SIVEP
  df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),] #4 is unexplained
  df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
  df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
  df_SIVEP = df_SIVEP[,c("DT_EVOLUCA", "SG_UF")] #Date, and Federative unit patient residence.
  df_SIVEP$Deaths = 1
  df_SIVEP = aggregate(. ~SG_UF+DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
  dim(df_SIVEP)
  df_SIVEP$Release = dmy(RELEASE_DATE_i)
  colnames(df_SIVEP) = c("State","Date","Deaths","Release")
  df_SIVEP_nowcast = rbind(df_SIVEP_nowcast,df_SIVEP)
}


###################################################################################

head(df_SIVEP_nowcast)
dim(df_SIVEP_nowcast)
plot(df_SIVEP_nowcast[,c("Date","Deaths")])
lines(df_SIVEP_nowcast[,c("Date","Deaths")])

df_SIVEP_nowcast$Date_index = as.numeric(df_SIVEP_nowcast$Date - ymd("2020-01-01"))
df_SIVEP_nowcast$Release_index = as.numeric(df_SIVEP_nowcast$Release - dmy(RELEASE_DATE[[1]]))

write_csv(path = paste0(path, "df_SIVEP_nowcast_allStates_", tail(RELEASE_DATE,1), ".csv"), df_SIVEP_nowcast)
write_csv(path = paste0(path2, "df_SIVEP_nowcast_allStates_", tail(RELEASE_DATE,1), ".csv"), df_SIVEP_nowcast)

####################################################################################