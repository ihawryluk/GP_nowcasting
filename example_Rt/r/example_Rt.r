library(rstan)
library(matrixStats)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(tidyverse)
library(dplyr)
library(abind)
library(xtable)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(openxlsx)
library(loo)


###########################################################################
### CHANGE THE PATHS BEFORE RUNNING THE SCRIPT ####
###########################################################################

source("utils/geom-stepribbon.r")
source("utils/gammaAlt.r")
path <- paste0("/Users/name/Documents/nowcasting/example_Rt/")
df_SIVEP_original = read_csv2(paste0(path,"data/","INFLUD-23-11-2020.csv"))

if(!file.exists(paste0(path,"data/","INFLUD-23-11-2020.csv"))){
  download.file(paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2020/INFLUD-3-11-2020.csv"),df_SIVEP_original)
}

selected_state <- "Brazil"
job = "example_Rt" 
makeFigDir <- paste0("if [ ! -f figures/example_Rt ]; then mkdir ",path,"figures/",job," ; fi")
makeResDir <- paste0("if [ ! -f results/example_Rt ]; then mkdir ",path,"results/",job," ; fi")
system(makeFigDir)
system(makeResDir)
name="_longer-test5_"
final_date = ymd("2020-11-22")
ITER = 2000
WARM = 1000
CORES = 4
SEED = 999

get_deaths <- function(arg) {
    df_SIVEP_original -> df_SIVEP
    df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5) || which(df_SIVEP$CLASSI_FIN==4),]
    df_SIVEP = df_SIVEP[which(df_SIVEP$EVOLUCAO==2),]
    df_SIVEP$DT_EVOLUCA = dmy(df_SIVEP$DT_EVOLUCA)
    filter_date = head(sort(df_SIVEP[df_SIVEP$CLASSI_FIN==5,]$DT_EVOLUCA),1)
    df_SIVEP = df_SIVEP %>% filter(DT_EVOLUCA >= ymd(filter_date))
    df_SIVEP = df_SIVEP[,c("DT_EVOLUCA")]
    dim(df_SIVEP)
    df_SIVEP$Deaths = 1
    df_SIVEP = aggregate(. ~DT_EVOLUCA, data=df_SIVEP, sum, na.rm=TRUE)
    sum(df_SIVEP$Deaths)
    colnames(df_SIVEP) = c("DateRep","Deaths")
    df_SIVEP$region = selected_state
    df_SIVEP = rbind(df_SIVEP[df_SIVEP$DateRep < ymd(head(arg,1)$DateRep),],arg)
    sum(df_SIVEP$Deaths)
    regions <- unique(df_SIVEP$region)
    aux_df_zeros_brazil = NULL
    for( jj in 1:length(regions)){
      aux_sub = subset(df_SIVEP,region==regions[jj])
      aux_sub = aux_sub[order(aux_sub$DateRep),]
      last_date = aux_sub$DateRep[1]
      zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                         region=regions[jj],
                         Deaths=0)
      zeroes$DateRep <- ymd(zeroes$DateRep)
      aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
    }

    df_SIVEP_deaths=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)
    head(df_SIVEP_deaths) 
  return(df_SIVEP_deaths)
}

deaths_sources <- c("Raw4and5","Nowcast4and5","GroundTruth4and5")

out_all = NULL
for ( deaths_source in deaths_sources )
{
    onset_paras=read.csv(paste0(path,'data/',"statelevel_OnsetDeathParams.csv"))
    onset_paras = onset_paras[c("state","mean","cv")]  %>% drop_na()
    head(onset_paras,1)

    deaths_nowcast_original <- read.csv(paste0("data/","nowcasted_daily_Brazil_23112020_2021-02-20.csv"),stringsAsFactors = FALSE)
    deaths_nowcast <- deaths_nowcast_original
    deaths_nowcast$"DateRep" = ymd(as.Date(deaths_nowcast$"Date"))
    deaths_nowcast = deaths_nowcast[c("DateRep","Deaths")]
    deaths_nowcast$region = as.factor(selected_state)

    deaths_groundtruth_original <- read.csv(paste0("data/","true_daily_Brazil_23112020.csv"),stringsAsFactors = FALSE)
    deaths_groundtruth <- deaths_groundtruth_original
    deaths_groundtruth$"DateRep" = ymd(as.Date(deaths_groundtruth$"Date"))
    deaths_groundtruth = deaths_groundtruth[c("DateRep","Deaths")]
    deaths_groundtruth$region = as.factor(selected_state)

    deaths_raw_original <- read.csv(paste0("data/","pre-nowcast_daily_Brazil_23112020.csv"),stringsAsFactors = FALSE)
    deaths_raw <- deaths_raw_original
    deaths_raw$"DateRep" = ymd(as.Date(deaths_raw$"Date"))
    deaths_raw = deaths_raw[c("DateRep","Deaths")]
    deaths_raw$region = as.factor(selected_state)


    if ( deaths_source == "Nowcast4and5" )
    {
        df_SIVEP_deaths <- get_deaths(deaths_nowcast)
    } else if ( deaths_source == "GroundTruth4and5" )
    {
        df_SIVEP_deaths <- get_deaths(deaths_groundtruth)    
    } else if ( deaths_source == "Raw4and5" )
    {
        df_SIVEP_deaths <- get_deaths(deaths_raw)      
    }

    df_SIVEP_original -> df_SIVEP
    df_SIVEP = df_SIVEP[which(df_SIVEP$CLASSI_FIN==5),]
    df_SIVEP = df_SIVEP[,c("DT_SIN_PRI")]
    df_SIVEP$DT_SIN_PRI = dmy(df_SIVEP$DT_SIN_PRI)
    df_SIVEP = df_SIVEP[order(df_SIVEP$DT_SIN_PRI),] %>% filter(DT_SIN_PRI >= ymd('2020-02-15')) #painel index case
    df_SIVEP$Cases = 1
    df_SIVEP = aggregate(. ~DT_SIN_PRI, data=df_SIVEP, sum, na.rm=TRUE)
    colnames(df_SIVEP) = c("DateRep","Cases")
    df_SIVEP$region = selected_state
    regions <- selected_state
    aux_df_zeros_brazil = NULL
    for( jj in 1:length(regions)){
      aux_sub = subset(df_SIVEP,region==regions[jj])
      aux_sub = aux_sub[order(aux_sub$DateRep),]
      last_date = aux_sub$DateRep[1]
      zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(last_date-1),"days"),
                         region=regions[jj],
                         Cases=0)
      zeroes$DateRep <- ymd(zeroes$DateRep)
      aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
    }
    df_SIVEP_cases=rbind.data.frame(aux_df_zeros_brazil,df_SIVEP)


    regions <- unique(df_SIVEP_deaths$region)
    aux_df_zeros_brazil = NULL
    for( jj in 1:length(regions)){
      zeroes<-data.frame(DateRep=seq(as.Date('2020/01/01'),as.Date(tail(sort(unique(df_SIVEP_deaths$DateRep)),1)),"days"),
                         region=regions[jj],
                         AUX=0)
      zeroes$DateRep <- ymd(zeroes$DateRep)
      aux_df_zeros_brazil = rbind.data.frame(zeroes,aux_df_zeros_brazil)
    }
    aux_df_zeros_brazil$key = paste0(aux_df_zeros_brazil$DateRep,aux_df_zeros_brazil$region)

    df_SIVEP_deaths$key = paste0(df_SIVEP_deaths$DateRep,df_SIVEP_deaths$region)
    df_SIVEP_cases$key = paste0(df_SIVEP_cases$DateRep,df_SIVEP_cases$region)
    df_SIVEP_deaths = merge(aux_df_zeros_brazil, df_SIVEP_deaths, by = "key", all=TRUE)
    df_SIVEP_cases = merge(aux_df_zeros_brazil, df_SIVEP_cases, by = "key", all=TRUE)

    df_SIVEP = merge(df_SIVEP_deaths,df_SIVEP_cases, by = "key")
    df_SIVEP = df_SIVEP[,c("DateRep.y.x","region.y.x","Deaths","Cases")]
    colnames(df_SIVEP) = c("DateRep","region","Deaths","Cases")
    df_SIVEP[c("Cases")] = df_SIVEP[c("Cases")] %>% replace(is.na(.), 0)

    df_SIVEP = df_SIVEP %>% group_by(region) %>% mutate(cumulativeCases = cumsum(Cases))
    df_SIVEP = df_SIVEP %>% group_by(region) %>% mutate(cumulativeDeaths = cumsum(Deaths))


    df_pop <- NULL
    df_pop<-read.csv(paste0(path,'data/',"brazil_population.csv"),sep=";", stringsAsFactors = FALSE)
    df_pop<-df_pop[c("region","population")]
    df_SIVEP = merge(x = df_SIVEP, y = df_pop, by = "region", all = TRUE)
    df_SIVEP = df_SIVEP[order(df_SIVEP$DateRep),]

    df=df_SIVEP %>% filter(DateRep <= final_date)
    START_TIME=ymd(as.Date("2020-02-15"))
    END_TIME=final_date
    RANGE_TIME=seq(START_TIME,END_TIME,by = '1 day')
    StanModel = as.character(job)


    models = job
    writeDir = job
    serial.interval = read.csv("data/serial_interval.csv")
    d<-df
    countries = selected_state
    N2 = length(RANGE_TIME)

    dates = list()
    reported_cases = list()
    stan_data = list(M=length(countries),N=NULL,
                     deaths=NULL,deaths_combined=NULL,f=NULL,
                   N0=6,cases=NULL,SI=NULL,
                   EpidemicStart = NULL,
                   pop = NULL)

    deaths_by_country = list()
    deaths_by_country_combined = list()
    mean1 = 5.1; cv1 = 0.86;
    x1 = rgammaAlt(1e6,mean1,cv1)
    aux.epidemicStart = NULL
    Country = selected_state
    mean2 = onset_paras[onset_paras$state == Country,]$mean
    cv2 = onset_paras[onset_paras$state == Country,]$cv
    x2 = rgammaAlt(1e6,mean2,cv2) # onset-to-death distribution
    ecdf.saved = ecdf(x1+x2)


    IFR = 0.5/100
    d1=d[d$region==Country,]
    d1$DateRep = as.Date(d1$DateRep)
    d1_pop = df_pop[df_pop$region==Country,]
    d1 = d1[order(d1$DateRep),]

    index1 = which(cumsum(d1$Deaths)>=20)[1] # 5, 10
    index2 = index1-26 #26
    d1=d1[index2:nrow(d1),]

    aux.epidemicStart = c(aux.epidemicStart,d1$DateRep[index1+1-index2])
    stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
    stan_data$EpidemicStart = as.array(stan_data$EpidemicStart)
    stan_data$pop = c(stan_data$pop, d1_pop$population)
    stan_data$pop = as.array(stan_data$pop)
    dates[[Country]] = d1$DateRep
    N = length(d1$Cases)
    print(sprintf("%s has %d days of data",Country,N))
    forecast = N2 - N
    if(forecast < 0) {
      print(sprintf("%s: %d", Country, N))
      print("ERROR!!!! increasing N2")
      N2 = N
      forecast = N2 - N
    }

    convolution = function(u) (IFR * ecdf.saved(u))   
    f = rep(0,N2) # f is the probability of dying on day i given infection
    f[1] = (convolution(1.5) - convolution(0))
    for(i in 2:N2) 
    {
          f[i] = (convolution(i+.5) - convolution(i-.5))
    }


    cases = as.vector(as.numeric(d1$Cases))
    deaths = as.vector(as.numeric(d1$Deaths))
    stan_data$N = c(stan_data$N,N)
    stan_data$f = cbind(stan_data$f,f)
    stan_data$deaths = cbind(stan_data$deaths,deaths)
    stan_data$cases = cbind(stan_data$cases,cases)
    stan_data$N2=N2
    stan_data$x=1:N2
    if(length(stan_data$N) == 1) {
      stan_data$N = as.array(stan_data$N)
    }
    stan_data$W <- ceiling(stan_data$N2/7)
    stan_data$week_index <- matrix(1,stan_data$M,stan_data$N2)
        for(state.i in 1:stan_data$M) {
        stan_data$week_index[state.i,] <- rep(2:(stan_data$W+1),each=7)[1:stan_data$N2]
        last_ar_week = which(dates[[state.i]]==max(df$Date) -  7)
        stan_data$week_index[state.i,last_ar_week:ncol(stan_data$week_index)] <-
        stan_data$week_index[state.i,last_ar_week]
    }
    stan_data$AR_SD_MEAN = 0.2
    stan_data$P <- 0
    if(length(stan_data$deaths < N2))
    {
        stan_data$deaths = rbind(stan_data$deaths,as.matrix(rep(-1,N2-N)))
    }
    stan_data$SI <- serial.interval$fit[1:N2]



    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
    m = stan_model(paste0(path,'stan_models/',StanModel,'.stan'))




    fit = sampling(m,
                   data=stan_data,
                   iter=ITER,
                   warmup=WARM,
                   seed = SEED,
                   chains=CORES,
                   control = list(adapt_delta = 0.99, max_treedepth = 12))
    out = rstan::extract(fit)
    out_all[[deaths_source]] = out
}



get_plots <- function(out,date_last_plot) {
    prediction = out$prediction
    estimated.deaths = out$E_deaths
    estimated.immune = out$immune
    JOBID = as.character(abs(round(rnorm(1) * 1000000)))
    print(sprintf("Jobid = %s",JOBID))
    filename <- paste0(StanModel,'-',JOBID)

    region=c("DF","AC","AL","AP","AM","BA","CE","ES","GO","MA","MT","MS","MG","PR","PB","PA","PE","PI","RN","RS","RJ","RO","RR","SC","SE","SP","TO")
    county=c("Federal District",'State of Acre','State of Alagoas','State of Amapá','State of Amazonas','State of Bahia','State of Ceará','State of Espírito Santo','State of Goiás','State of Maranhão','State of Mato Grosso','State of Mato Grosso do Sul','State of Minas Gerais','State of Paraná','State of Paraíba','State of Pará','State of Pernambuco','State of Piauí','State of Rio Grande do Norte','State of Rio Grande do Sul','State of Rio de Janeiro','State of Rondônia','State of Roraima','State of Santa Catarina','State of Sergipe','State of São Paulo','State of Tocantins')
    df_region_codes = data.frame(region,county)

    table_paper = NULL
    i=1
    reported_cases = cases
    country <- countries[[i]]
    dates_forecast = dates[[i]]
    N3 = length(dates_forecast)
    reported_cases_padded<- c(reported_cases,rep(0,N3-N))
    deaths_by_country_padded<-c(deaths,rep(0,N3-N))
    predicted_cases <- colMeans(prediction[,1:N3,i])
    predicted_cases_li <- colQuantiles(prediction[,1:N3,i], probs=.025)
    predicted_cases_ui <- colQuantiles(prediction[,1:N3,i], probs=.975)
    predicted_cases_li2 <- colQuantiles(prediction[,1:N3,i], probs=.25)
    predicted_cases_ui2 <- colQuantiles(prediction[,1:N3,i], probs=.75)
    estimated_deaths <- colMeans(estimated.deaths[,1:N3,i])
    estimated_deaths_li <- colQuantiles(estimated.deaths[,1:N3,i], probs=.025)
    estimated_deaths_ui <- colQuantiles(estimated.deaths[,1:N3,i], probs=.975)
    estimated_deaths_li2 <- colQuantiles(estimated.deaths[,1:N3,i], probs=.25)
    estimated_deaths_ui2 <- colQuantiles(estimated.deaths[,1:N3,i], probs=.75)
    rt <- colMeans(out$Rt_adj[,1:N3,i])
    rt_li <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.025)
    rt_ui <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.975)
    rt_li2 <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.25)
    rt_ui2 <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.75)

    rt_unadj <- colMeans(out$Rt[,1:N3,i])
    rt_li_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.025)
    rt_ui_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.975)
    rt_li2_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.25)
    rt_ui2_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.75)

    data_country = NULL
    table_paper = NULL
    data_country <- data.frame("time" = as_date(as.character(dates_forecast)),
                           "country" = rep(country, N3),
                           "reported_cases" = reported_cases_padded,
                           "reported_cases_c" = cumsum(reported_cases_padded),

                            "predicted_cases_c" = cumsum(predicted_cases),
                            "predicted_min_c" = cumsum(predicted_cases_li),
                            "predicted_max_c" = cumsum(predicted_cases_ui),
                            "predicted_cases" = predicted_cases,
                            "predicted_min" = predicted_cases_li,
                            "predicted_max" = predicted_cases_ui,
                            "predicted_min2" = predicted_cases_li2,
                            "predicted_max2" = predicted_cases_ui2,

                            "deaths" = deaths_by_country_padded,
                            "deaths_c" = cumsum(deaths_by_country_padded),
                            "estimated_deaths_c" =  cumsum(estimated_deaths),
                            "death_min_c" = cumsum(estimated_deaths_li),
                            "death_max_c"= cumsum(estimated_deaths_ui),
                            "estimated_deaths" = estimated_deaths,

                            "death_min" = estimated_deaths_li,
                            "death_max"= estimated_deaths_ui,
                            "death_min2" = estimated_deaths_li2,
                            "death_max2"= estimated_deaths_ui2,

                            "rt" = rt,
                            "rt_min" = rt_li,
                            "rt_max" = rt_ui,
                            "rt_min2" = rt_li2,
                            "rt_max2" = rt_ui2,

                            "rt_unadj" = rt_unadj,
                            "rt_min_unadj" = rt_li_unadj,
                            "rt_max_unadj" = rt_ui_unadj,
                            "rt_min2_unadj" = rt_li2_unadj,
                            "rt_max2_unadj" = rt_ui2_unadj
                          )

        data_deaths_95 <- data.frame(data_country$time, 
                                     data_country$deaths,
                                     data_country$death_min,
                                     data_country$death_max)
        names(data_deaths_95) <- c("time", "deaths",
                                   "death_min",
                                    "death_max")
        data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$time))
        data_deaths_50 <- data.frame(data_country$time, 
                                     data_country$deaths,
                                     data_country$death_min2,
                                     data_country$death_max2)

        names(data_deaths_50) <- c("time", "deaths",
                                   "death_min",
                                    "death_max")

    data_deaths_50$key <- (rep("fifty", length(data_deaths_50$time)))
    data_deaths <- rbind(data_deaths_95, data_deaths_50)
    levels(data_deaths$key) <- c("ninetyfive", "fifty")+ coord_fixed(ratio = 10)
        data_deaths <- data_deaths[data_deaths$time > date_last_plot,]
    p2 <-  data_deaths %>% ggplot(aes(x = time)) +
          geom_point(aes(y = deaths, col = "coral4", alpha=1)) +
            geom_ribbon(aes(ymin = death_min, ymax = death_max, fill = key)) +
            scale_x_date(date_breaks = "4 weeks", labels = date_format("%e %b")) +
            scale_y_continuous(expand = c(0, 0), labels = comma) +
            scale_fill_manual(name = "", labels = c("50%", "95%"),
                            values = c(alpha("deepskyblue4", 0.55),
                                       alpha("deepskyblue4", 0.45))) +
            ylab("Daily number of deaths\n") +
            xlab("") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "None") +
            guides(fill=guide_legend(ncol=1))

    #################################################
        data_rt_95 <- data.frame(data_country$time,
                                 data_country$rt_min, 
                                 data_country$rt_max)
        names(data_rt_95) <- c("time", "rt_min", "rt_max")
        data_rt_95$key <- rep("nintyfive", length(data_rt_95$time))
        data_rt_50 <- data.frame(data_country$time, 
                                 data_country$rt_min2,
                                 data_country$rt_max2)
        names(data_rt_50) <- c("time", "rt_min", "rt_max")
        data_rt_50$key <- rep("fifty", length(data_rt_50$time))
        data_rt <- rbind(data_rt_95, data_rt_50)
        levels(data_rt$key) <- c("ninetyfive", "fifth")
    data_rt <- data_rt[data_rt$time > date_last_plot,]
    p3 <- data_rt %>% ggplot() +
          geom_ribbon(aes(x = time, ymin = rt_min, ymax = rt_max,
                                          group = key,
                                          fill = key)) +
          geom_hline(yintercept = 1, color = 'black', size = 0.1) +
          xlab("") +
          ylab(expression(R[t])) +
          scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b"),
                       limits = c(data_rt$time[1],
                                  data_rt$time[length(data_rt$time)])) +
          scale_fill_manual(name = "", labels = c("50%", "95%","50% v2", "95% v2"),
                            values = c(alpha("red", 0.75), alpha("red", 0.5),
                                       alpha("seagreen", 0.75), alpha("seagreen", 0.5))) +
        ylim(0.1,2.5)+
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(legend.position="right")+
        theme(text = element_text(size=22))
            
        print(p2)
        print(p3)
}


get_plots(out_all[[deaths_sources[[1]]]],ymd("2020-10-01"))
get_plots(out_all[[deaths_sources[[2]]]],ymd("2020-10-01"))
get_plots(out_all[[deaths_sources[[3]]]],ymd("2020-10-01"))


out = out_all[[deaths_sources[[3]]]]
date_last_plot = ymd("2020-10-01")
    prediction = out$prediction
    estimated.deaths = out$E_deaths
    estimated.immune = out$immune
    JOBID = as.character(abs(round(rnorm(1) * 1000000)))
    print(sprintf("Jobid = %s",JOBID))
    filename <- paste0(StanModel,'-',JOBID)

    region=c("DF","AC","AL","AP","AM","BA","CE","ES","GO","MA","MT","MS","MG","PR","PB","PA","PE","PI","RN","RS","RJ","RO","RR","SC","SE","SP","TO")
    county=c("Federal District",'State of Acre','State of Alagoas','State of Amapá','State of Amazonas','State of Bahia','State of Ceará','State of Espírito Santo','State of Goiás','State of Maranhão','State of Mato Grosso','State of Mato Grosso do Sul','State of Minas Gerais','State of Paraná','State of Paraíba','State of Pará','State of Pernambuco','State of Piauí','State of Rio Grande do Norte','State of Rio Grande do Sul','State of Rio de Janeiro','State of Rondônia','State of Roraima','State of Santa Catarina','State of Sergipe','State of São Paulo','State of Tocantins')
    df_region_codes = data.frame(region,county)

    table_paper = NULL
    i=1
    reported_cases = cases
    country <- countries[[i]]
    dates_forecast = dates[[i]]
    N3 = length(dates_forecast)
    reported_cases_padded<- c(reported_cases,rep(0,N3-N))
    deaths_by_country_padded<-c(deaths,rep(0,N3-N))
    predicted_cases <- colMeans(prediction[,1:N3,i])
    predicted_cases_li <- colQuantiles(prediction[,1:N3,i], probs=.025)
    predicted_cases_ui <- colQuantiles(prediction[,1:N3,i], probs=.975)
    predicted_cases_li2 <- colQuantiles(prediction[,1:N3,i], probs=.25)
    predicted_cases_ui2 <- colQuantiles(prediction[,1:N3,i], probs=.75)
    estimated_deaths <- colMeans(estimated.deaths[,1:N3,i])
    estimated_deaths_li <- colQuantiles(estimated.deaths[,1:N3,i], probs=.025)
    estimated_deaths_ui <- colQuantiles(estimated.deaths[,1:N3,i], probs=.975)
    estimated_deaths_li2 <- colQuantiles(estimated.deaths[,1:N3,i], probs=.25)
    estimated_deaths_ui2 <- colQuantiles(estimated.deaths[,1:N3,i], probs=.75)
    rt <- colMeans(out$Rt_adj[,1:N3,i])
    rt_li <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.025)
    rt_ui <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.975)
    rt_li2 <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.25)
    rt_ui2 <- colQuantiles(out$Rt_adj[,1:N3,i],probs=.75)

    rt_unadj <- colMeans(out$Rt[,1:N3,i])
    rt_li_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.025)
    rt_ui_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.975)
    rt_li2_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.25)
    rt_ui2_unadj <- colQuantiles(out$Rt[,1:N3,i],probs=.75)

    data_country = NULL
    table_paper = NULL
    data_country <- data.frame("time" = as_date(as.character(dates_forecast)),
                           "country" = rep(country, N3),
                           "reported_cases" = reported_cases_padded,
                           "reported_cases_c" = cumsum(reported_cases_padded),

                            "predicted_cases_c" = cumsum(predicted_cases),
                            "predicted_min_c" = cumsum(predicted_cases_li),
                            "predicted_max_c" = cumsum(predicted_cases_ui),
                            "predicted_cases" = predicted_cases,
                            "predicted_min" = predicted_cases_li,
                            "predicted_max" = predicted_cases_ui,
                            "predicted_min2" = predicted_cases_li2,
                            "predicted_max2" = predicted_cases_ui2,


                            "deaths" = deaths_by_country_padded,
                            "deaths_c" = cumsum(deaths_by_country_padded),
                            "estimated_deaths_c" =  cumsum(estimated_deaths),
                            "death_min_c" = cumsum(estimated_deaths_li),
                            "death_max_c"= cumsum(estimated_deaths_ui),
                            "estimated_deaths" = estimated_deaths,

                            "death_min" = estimated_deaths_li,
                            "death_max"= estimated_deaths_ui,
                            "death_min2" = estimated_deaths_li2,
                            "death_max2"= estimated_deaths_ui2,

                            "rt" = rt,
                            "rt_min" = rt_li,
                            "rt_max" = rt_ui,
                            "rt_min2" = rt_li2,
                            "rt_max2" = rt_ui2,

                            "rt_unadj" = rt_unadj,
                            "rt_min_unadj" = rt_li_unadj,
                            "rt_max_unadj" = rt_ui_unadj,
                            "rt_min2_unadj" = rt_li2_unadj,
                            "rt_max2_unadj" = rt_ui2_unadj
                          )

        data_deaths_95 <- data.frame(data_country$time, 
                                     data_country$deaths,
                                     data_country$death_min,
                                     data_country$death_max)
        names(data_deaths_95) <- c("time", "deaths",
                                   "death_min",
                                    "death_max")
        data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$time))
        data_deaths_50 <- data.frame(data_country$time, 
                                     data_country$deaths,
                                     data_country$death_min2,
                                     data_country$death_max2)

        names(data_deaths_50) <- c("time", "deaths",
                                   "death_min",
                                    "death_max")

    data_deaths_50$key <- (rep("fifty", length(data_deaths_50$time)))
    data_deaths <- rbind(data_deaths_95, data_deaths_50)
    levels(data_deaths$key) <- c("ninetyfive", "fifty")+ coord_fixed(ratio = 10)
        data_deaths <- data_deaths[data_deaths$time > date_last_plot,]
    p2 <-  data_deaths %>% ggplot(aes(x = time)) +
          geom_point(aes(y = deaths, col = "coral4", alpha=1)) +
            geom_ribbon(aes(ymin = death_min, ymax = death_max, fill = key)) +
            scale_x_date(date_breaks = "4 weeks", labels = date_format("%e %b")) +
            scale_y_continuous(expand = c(0, 0), labels = comma) +
            scale_fill_manual(name = "", labels = c("50%", "95%"),
                            values = c(alpha("deepskyblue4", 0.55),
                                       alpha("deepskyblue4", 0.45))) +
            ylab("Daily number of deaths\n") +
            xlab("") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "None") +
            guides(fill=guide_legend(ncol=1))

    #################################################
        data_rt_95 <- data.frame(data_country$time,
                                 data_country$rt_min, 
                                 data_country$rt_max)
        names(data_rt_95) <- c("time", "rt_min", "rt_max")
        data_rt_95$key <- rep("nintyfive", length(data_rt_95$time))
        data_rt_50 <- data.frame(data_country$time, 
                                 data_country$rt_min2,
                                 data_country$rt_max2)
        names(data_rt_50) <- c("time", "rt_min", "rt_max")
        data_rt_50$key <- rep("fifty", length(data_rt_50$time))
        data_rt <- rbind(data_rt_95, data_rt_50)
        levels(data_rt$key) <- c("ninetyfive", "fifth")
    data_rt <- data_rt[data_rt$time > date_last_plot,]
    p3 <- data_rt %>% ggplot() +
          geom_ribbon(aes(x = time, ymin = rt_min, ymax = rt_max,
                                          group = key,
                                          fill = key)) +
          geom_hline(yintercept = 1, color = 'black', size = 0.1) +
          xlab("") +
          ylab(expression(R[t])) +
          scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b"),
                       limits = c(data_rt$time[1],
                                  data_rt$time[length(data_rt$time)])) +
          scale_fill_manual(name = "", labels = c("50%", "95%","50% v2", "95% v2"),
                            values = c(alpha("red", 0.75), alpha("red", 0.5),
                                       alpha("seagreen", 0.75), alpha("seagreen", 0.5))) +
        ylim(0.1,1.8)+
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(legend.position="right")+
        theme(text = element_text(size=20))+
        theme(legend.text=element_text(size=20))


deaths_raw_original <- read.csv(paste0("data/","pre-nowcast_daily_Brazil_23112020.csv"),stringsAsFactors = FALSE)
    deaths_raw <- deaths_raw_original
    deaths_raw$"DateRep" = ymd(as.Date(deaths_raw$"Date"))
    deaths_raw = deaths_raw[c("DateRep","Deaths")]
    deaths_raw$region = as.factor(selected_state)
rawplt = data.frame(
    "date"=deaths_raw$DateRep,
    "deaths"=deaths_raw$Deaths
    )
rawplt = rawplt[rawplt$date >date_last_plot,]


    data_deaths_95 <- data.frame(data_country$time, 
                                 data_country$deaths,
                                 data_country$death_min,
                                 data_country$death_max)
    names(data_deaths_95) <- c("time", "deaths",
                               "death_min",
                                "death_max")
    data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$time))
    data_deaths_50 <- data.frame(data_country$time, 
                                 data_country$deaths,
                                 data_country$death_min2,
                                 data_country$death_max2)

    names(data_deaths_50) <- c("time", "deaths",
                               "death_min",
                                "death_max")

data_deaths_50$key <- (rep("fifty", length(data_deaths_50$time)))
data_deaths <- rbind(data_deaths_95, data_deaths_50)
levels(data_deaths$key) <- c("ninetyfive", "fifty")+ coord_fixed(ratio = 10)
    data_deaths <- data_deaths[data_deaths$time > date_last_plot,]
p2_data <-  data_deaths %>% ggplot(aes(x = time)) +
        geom_point(data = rawplt, aes(x = date, y = deaths,size =1.5,alpha=1)) +
        geom_point(aes(y = deaths, col = "coral4", alpha=1,size =1.5)) +
        scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b")) +
        scale_y_continuous(expand = c(0, 0), labels = comma) +
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                        values = c(alpha("deepskyblue4", 0.55),
                                   alpha("deepskyblue4", 0.45))) +
        ylab("Daily number of deaths\n") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "None") +
        guides(fill=guide_legend(ncol=1)) +
        theme(text = element_text(size=22))

p2_data


kl <- function(outB, outA) {
    mu1 = colMeans(outA$Rt_adj[,1:N3,i])
    sd1 = colSds(outA$Rt_adj[,1:N3,i])
    mu2 = colMeans(outB$Rt[,1:N3,i])
    sd2 = colSds(outB$Rt[,1:N3,i])
    kl = log(sd2/sd1) + (sd1*sd1 + (mu1 - mu2) * (mu1 - mu2))/(2 * sd2 * sd2) - 0.5
    kl
}
out1 = out_all[[deaths_sources[[1]]]]
out2 = out_all[[deaths_sources[[2]]]]
out3 = out_all[[deaths_sources[[3]]]]
plot(kl(out1,out3))
plot(kl(out2,out3))

kl_df = data.frame("dates" = dates[[1]], "kl_raw_true" = kl(out1,out3),"kl_pred_true" = kl(out2,out3))
tail(kl_df,70) %>% ggplot(aes(x = dates)) + 
        geom_line(aes(y = kl_raw_true),col="coral4",size=2) + 
        geom_line(aes(y = kl_pred_true),col="deepskyblue4",size=2) + 
        scale_x_date(date_breaks = "2 weeks", labels = date_format("%e %b")) +
        scale_y_continuous(expand = c(0, 0), labels = comma) +
        scale_fill_manual(name = "", labels = c("50%", "95%"),
                        values = c(alpha("deepskyblue4", 0.55),
                                   alpha("deepskyblue4", 0.45))) +
        ylab("KL(. || ground truth)\n") +
        xlab("") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "None") +
        guides(fill=guide_legend(ncol=1)) +
        theme(text = element_text(size=22))
