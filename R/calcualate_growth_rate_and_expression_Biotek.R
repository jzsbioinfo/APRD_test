
#' calcualate_growth_rate_and_expression_Biotek
#' #'
#' A function use output of plate reader to caluclate growth rate, if have fluorescence data, it can also calculate fluorescence expression level.
#' @param inputfile inputfile="./dat-001-001-sd-p121.123.3AT_SCD-H.xlsx";
#' @param skip_num skip_num=40, how many lines before the real data
#' @param plate plate=96, means 96 well plate, plate=384 means 384 well plate.
#' @param read_times_per_hour read_times_per_hour=6, how many points in one hour.
#' @param total_read_times total_read_times=217, total points in this experiment.
#' @param blank_wells blank_wells=c(1:12); the wells only with medimum as control. For example, A01:A12 (first row) is 1:12 and A01:H01 (first column) is c(1, 13, 25, 37, 49, 61, 73, 85)
#' @param X_hours_SW X_hours_SW=2 The time length (hours) of sliding window, need to be less than the total time. The total time points used is X_hours_SW*read_times_per_hour. Suggest to make X_hours_SW\*read_times_per_hour >= 10
#' @param medimum1 Defaults to medimum1="carbon"; used to label samples.
#' @param medimum2 Defaults to medimum2="nitrogen"; used to label samples.
#' @param outputdir outputdir="./20200303output/";  the folder you want to output
#' @param outputfile outputfile= "20200303saad.xlsx"; xlsx format, the growth rate result file
#' @param outputfile2 outputfile2 = "./20190731test96fi_2.xlsx"; xlsx format, the fluoresence result file, if no, set like this.
#' @keywords growth_rate
#' @export
#' @examples
#' calcualate_growth_rate(inputfile,skip_num, plate,read_times_per_hour,total_read_times, blank_wells,X_hours_SW,medimum1="carbon", medimum2="nitrogen",outputdir,outputfile,outputfile2)




# read in biotek data calculate growth rate and FI---------------------------
calcualate_growth_rate_and_expression_Biotek <- function(inputfile,skip_num, plate,
                                                  read_times_per_hour,
                                                  total_read_times, blank_wells,
                                                  X_hours_SW,
                                                  medimum1="carbon", medimum2="nitrogen",
                                                  outputdir,outputfile,outputfile2) {

  options(stringsAsFactors = F)
  # load packages -----
  library(tidyverse)
  library(dplyr)
  library(readxl)
  library(stringr)
  library(platetools)
  library(beeswarm)
  library(data.table)
  library(lubridate) # this is for time format
  library(chron)
  library(openxlsx)
  library(zoo)

  # read date, experimental file path, protocol file path------------
  date <- read_excel(inputfile,col_names=F,range="B7",col_types = c("text"))
  start_time <- read_excel(inputfile,col_names=F,range="B8",col_types = c("text"))
  experimental_file_path <- read_excel(inputfile,col_names=F,range="B4",col_types = c("text"))
  protocol_file_path <- read_excel(inputfile,col_names=F,range="B5",col_types = c("text"))

  date <- as.Date(as.numeric(date), origin = "1899-12-30")
  start_time <- times(as.numeric(start_time))
  experimental_file_path <- as.character(experimental_file_path)
  protocol_file_path <-  as.character(protocol_file_path)

  # read raw data and transfrom to data.frame---------------
  #if(shake=="Y") skip_num <- 39 else skip_num <- 37

  to_num <- function(input) {
    as.data.frame(lapply(input, function(x) {
      as.numeric(as.character(x))
    }))
  }

  rowmedain <-  function(x) {
    if (is.matrix(x) | is.data.frame(x) | is.data.table(x)) {
      apply(x, 1, median)
    } else {
      median(x)
    }

  }


  # if it is a 96 well plate-----------------------------
  if (plate==96) {
    library(growthcurver)
    biotek <- read_excel(inputfile,skip = skip_num,sheet = 1)
    biotek <- biotek[,-c(1,3)]
    biotek <- na.omit(biotek)
    # read total_read_times
    # -1 is the header
    biotek_OD600 <- biotek[1:total_read_times,]
    biotek_GFP <- biotek[(total_read_times+1):(2*total_read_times+1),][-1,]
    biotek_mCherry <- biotek[(2*total_read_times+2):(3*total_read_times+2),][-1,]
    biotek_wassabi <- biotek[(3*total_read_times+3):(4*total_read_times+3),][-1,]


    gap_time <- 60/read_times_per_hour


    biotek_OD600$Time <- seq(0,gap_time*total_read_times-gap_time,gap_time)/60
    biotek_GFP$Time <- seq(0,gap_time*total_read_times-gap_time,gap_time)/60
    biotek_mCherry$Time <- seq(0,gap_time*total_read_times-gap_time,gap_time)/60
    biotek_wassabi$Time <- seq(0,gap_time*total_read_times-gap_time,gap_time)/60


    biotek_OD600 <- to_num(biotek_OD600)
    biotek_GFP <- to_num(biotek_GFP)
    biotek_mCherry <- to_num(biotek_mCherry)
    biotek_wassabi <- to_num(biotek_wassabi)


    # minus the blank value
    if ( ! blank_wells=="N" ) { biotek_OD600[,-1] <- biotek_OD600[,-1]- rowmedain(biotek_OD600[,blank_wells+1]) }

    if ( ! blank_wells=="N" ) { biotek_GFP[,-1] <- biotek_GFP[,-1]- rowmedain(biotek_GFP[,blank_wells+1]) }

    if ( ! blank_wells=="N" ) { biotek_mCherry[,-1] <- biotek_mCherry[,-1]- rowmedain(biotek_mCherry[,blank_wells+1]) }

    if ( ! blank_wells=="N" ) { biotek_wassabi[,-1] <- biotek_wassabi[,-1]- rowmedain(biotek_wassabi[,blank_wells+1]) }

    # change colnames
    colnames(biotek_OD600)[-1] <- num_to_well(1:plate,plate = plate)
    colnames(biotek_GFP)[-1] <- num_to_well(1:plate,plate = plate)
    colnames(biotek_mCherry)[-1] <- num_to_well(1:plate,plate = plate)
    colnames(biotek_wassabi)[-1] <- num_to_well(1:plate,plate = plate)

    # change wells order from rows to columns
    biotek_OD600 <- biotek_OD600[,c(1,order(as.numeric(substring(colnames(biotek_OD600)[-1], 2)))+1)]

    biotek_GFP <- biotek_GFP[,c(1,order(as.numeric(substring(colnames(biotek_GFP)[-1], 2)))+1)]

    biotek_mCherry <- biotek_mCherry[,c(1,order(as.numeric(substring(colnames(biotek_mCherry)[-1], 2)))+1)]

    biotek_wassabi <- biotek_wassabi[,c(1,order(as.numeric(substring(colnames(biotek_wassabi)[-1], 2)))+1)]

    colnames(biotek_OD600)[1] <- "time"
    colnames(biotek_GFP)[1] <- "time"
    colnames(biotek_mCherry)[1] <- "time"
    colnames(biotek_wassabi)[1] <- "time"


    # output growthcurver figures

    biotek_OD600_plot <- SummarizeGrowthByPlate(biotek_OD600,
                                                plot_fit = TRUE,
                                                plot_file = paste(outputdir,
                                                                  "OD600_96well_growthcurver_fit_figure.pdf"))

    #dev.off()

    if (! is.na(biotek_GFP[1,2]) ) {
      biotek_GFP_plot <- SummarizeGrowthByPlate(biotek_GFP,plot_fit = TRUE,plot_file = paste(outputdir,"FI1_96well_growthcurver_fit_figure.pdf"))}

    if (! is.na(biotek_mCherry[1,2]) ) {
      biotek_mCherry_plot <- SummarizeGrowthByPlate(biotek_mCherry,plot_fit = TRUE,plot_file = paste(outputdir,"FI2_96well_growthcurver_fit_figure.pdf"))}

    if (! is.na(biotek_wassabi[1,2]) ) {
      biotek_mCherry_plot <- SummarizeGrowthByPlate(biotek_wassabi,plot_fit = TRUE,plot_file = paste(outputdir,"FI3_96well_growthcurver_fit_figure.pdf"))}
  }


  # if it is a 384 plate-----------------------------
  if (plate==384) {
    library(growthcurver384)
    biotek <- read_excel(inputfile,skip = skip_num,sheet = 1)
    biotek <- biotek[,-c(1,3)]
    biotek <- na.omit(biotek)
    # read total_read_times
    # -1 is the Time column
    biotek_OD600 <- cbind(biotek[1:total_read_times,],
                          biotek[(total_read_times+2):(2*total_read_times+1),][,-1],
                          biotek[(2*total_read_times+3):(3*total_read_times+2),][,-1],
                          biotek[(3*total_read_times+4):(4*total_read_times+3),][,-1])

    biotek_GFP <- cbind(biotek[(4*total_read_times+5):(5*total_read_times+4),],
                        biotek[((5*total_read_times+6)):(6*total_read_times+5),][,-1],
                        biotek[(6*total_read_times+7):(7*total_read_times+6),][,-1],
                        biotek[(7*total_read_times+8):(8*total_read_times+7),][,-1])

    biotek_mCherry <- cbind(biotek[(8*total_read_times+9):(9*total_read_times+8),],
                            biotek[((9*total_read_times+10)):(10*total_read_times+9),][,-1],
                            biotek[(10*total_read_times+11):(11*total_read_times+10),][,-1],
                            biotek[(11*total_read_times+12):(12*total_read_times+11),][,-1])

    biotek_OD600$Time <- seq(0,gap_time*total_read_times-gap_time,gap_time)/60
    biotek_GFP$Time <- seq(0,gap_time*total_read_times-gap_time,gap_time)/60
    biotek_mCherry$Time <- seq(0,gap_time*total_read_times-gap_time,gap_time)/60

    biotek_OD600 <- to_num(biotek_OD600)
    biotek_GFP <- to_num(biotek_GFP)
    biotek_mCherry <- to_num(biotek_mCherry)

    # minus the blank value
    if ( ! blank_wells=="N" ) { biotek_OD600[,-1] <- biotek_OD600[,-1]- rowmedain(biotek_OD600[,blank_wells+1]) }

    # change colnames
    colnames(biotek_OD600)[-1] <- num_to_well(1:plate,plate = plate)
    colnames(biotek_GFP)[-1] <- num_to_well(1:plate,plate = plate)
    colnames(biotek_mCherry)[-1] <- num_to_well(1:plate,plate = plate)


    # change wells order from rows to columns
    biotek_OD600 <- biotek_OD600[,c(1,order(as.numeric(substring(colnames(biotek_OD600)[-1], 2)))+1)]

    biotek_GFP <- biotek_GFP[,c(1,order(as.numeric(substring(colnames(biotek_GFP)[-1], 2)))+1)]

    biotek_mCherry <- biotek_mCherry[,c(1,order(as.numeric(substring(colnames(biotek_mCherry)[-1], 2)))+1)]


    colnames(biotek_OD600)[1] <- "time"
    colnames(biotek_GFP)[1] <- "time"
    colnames(biotek_mCherry)[1] <- "time"
    # output growthcurver figures



    biotek_OD600_plot <- SummarizeGrowthByPlate384(biotek_OD600,
                                                   plot_fit = TRUE,
                                                   plot_file = paste(Sys.Date(),
                                                                     "OD600_384well_growthcurver_fit_figure.pdf"))

    if (! is.na(biotek_GFP[1,2]) ) {
      biotek_GFP_plot <- SummarizeGrowthByPlate384(biotek_GFP,
                                                   plot_fit = TRUE,
                                                   plot_file = paste(Sys.Date(),
                                                                     "FI1_384well_growthcurver_fit_figure.pdf"))
    }

    if (! is.na(biotek_mCherry[1,2]) ) {
      biotek_mCherry_plot <- SummarizeGrowthByPlate384(biotek_mCherry,
                                                       plot_fit = TRUE,
                                                       plot_file = paste(Sys.Date(),
                                                                         "FI2_384well_growthcurver_fit_figure.pdf"))
    }




  }

  # X_hours_SW hours sliding window to calculate slopes-----------------
  # OD below 2^(-5)==0.03125 will be discarded
  get_slopes <- function(data) {

    data[data<=2^(-5)] <- 0.01

    # transfer to log2 scaled
    t_log2 <- function(x) {
      return(log2(x))
    }

    # return slopes of window size X hours, step is 1
    time <- data[,1]
    lmfit_slopes <- function(x) {
      temp <- vector()
      for (i in seq(1,length(x)-X_hours_SW*read_times_per_hour+1,1)) {
        re <- lm(col~time,data=data.frame("time"=time[i:(i+X_hours_SW*read_times_per_hour-1)],"col"=x[i:(i+X_hours_SW*read_times_per_hour-1)]))
        re <- as.numeric(re$coefficients[2])
        temp <- c(temp,re)
      }
      #temp <- temp[order(temp)][1:5]
      return(temp)
    }



    re <- apply(t_log2(data[,-1]), 2, lmfit_slopes)
    re <- as.data.frame(re)
    #re[1:(4*7),] <- 0
    return(re)
    # discard the data out of 95% interval

  }

  # use slopes to calculate grwoth rate
  calculate_growth_rate <- function(x) {

    remove_first_8_slopes <- function(remove) {
      for (i in 2:(length(remove)-2*read_times_per_hour+1)) {
        if (remove[i] > 0.01 & remove[i-1]<0.0001) {

          remove[i:(i+2*read_times_per_hour-1)] <- 0
          break
        }
      }

      return(remove)

    }


    x <- apply(x, 2, remove_first_8_slopes)


    top11_median <- function(x) {
      x <- x[order(x,decreasing=T)]
      x <- median(x[1:11])

    }


    median_slope <- apply(x, 2, top11_median)


    # return where is the highest median
    temp <- vector()
    for (i in 1:length(median_slope)) {
      re <- match(median_slope[i],x[,i])
      temp <- c(temp,re)
    }

    result <- median_slope
    result <- data.frame("sample"=names(result),"growth_rate"=as.numeric(result))

    custom_sort <- function(x){x[order(as.numeric(substring(x, 2)))]}

    result <- result[match(custom_sort(num_to_well(1:plate,plate = plate)),result$sample),]
    result$sample <- as.character(result$sample)
    return(result)
  }

  slopes <- as.data.frame(get_slopes(biotek_OD600))
  growth_rate <- calculate_growth_rate(slopes)


  # make the annotation to the result file-----------------

  # add columns and rows
  growth_rate$columns <- substring(growth_rate$sample,2,3)
  growth_rate$rows <- substring(growth_rate$sample,1,1)

  if (plate==384) {
    real_order <- function() {
      temp <- vector()
      a <- 1:24
      b <- 25:48
      for(i in 1:24) {

        temp <- c(temp,rep(a[i],8),rep(b[i],8))
      }
      return(temp)
    }
  } else {
    real_order <- function() {
      temp <- vector()
      for (i in 1:12) {
        temp <- c(temp,rep(i,8))
      }
      return(temp)
    }
  }

  growth_rate$real_order <- real_order()

  growth_rate <- growth_rate[order(growth_rate$real_order),]


  growth_rate$medimum1 <- medimum1[growth_rate$real_order]

  growth_rate$medimum2 <- medimum2[growth_rate$real_order]


  paste_corresponding <- function(x,y) {
    res <- vector()
    for (i in 1:length(x)) {
      temp <- paste0(x[i],"+")
      temp <- paste0(temp,y[i])
      res <- c(res,temp)
    }
    return(res)
  }


  growth_rate$medimun1_medimum2 <- paste_corresponding(growth_rate$medimum1,growth_rate$medimum2)

  growth_rate$doubling_time.h <- 1/growth_rate$growth_rate

  growth_rate$date <- date

  growth_rate$start_time <- start_time

  growth_rate$experimental_file_path <- experimental_file_path

  growth_rate$protocol_file_path <- protocol_file_path

  growth_rate$read_times_per_hour <- read_times_per_hour

  growth_rate$total_read_times <- total_read_times

  growth_rate$blank_wells <- paste(blank_wells,sep='_',collapse = '_')

  growth_rate <- growth_rate[,c(1,4,3,5,6,7,8,2,9,10:16)]

  write.xlsx(growth_rate, outputfile, sheetName=paste(Sys.Date(),"result"))

  # output result
  #write.xlsx(growth_rate, outputfile, sheetName=paste(Sys.Date(),"result"))

  #########################################################################################

  # remove blanks to plot--------------

  growth_rate$real_order <- factor(growth_rate$real_order,levels=unique(growth_rate$real_order))

  plot_growth_rate <- growth_rate

  if (! blank_wells=="N") {plot_growth_rate <- plot_growth_rate[-blank_wells,]}

  # plot it

  png(paste0(outputdir,"doubling_time.png"),width=3600, height=1800)

  oldpar <- par()
  par(bg = "gray")
  beeswarm(doubling_time.h ~ real_order, data = plot_growth_rate,
           ylim=c(0,6),
           cex = 2,
           log = F, pch = 16, col = rainbow(12),
           main = paste(Sys.Date(),"MP_biotek"),
           xlab="",
           cex.axis=1,
           cex.lab=1.5,
           cex.main=3)
  bxplot(doubling_time.h ~ real_order, data = growth_rate, add = TRUE)
  par(oldpar)

  dev.off()

  #return(growth_rate)



  # (1) use 9 point window to calculate slope to find the maximum first derivation------------------------------

  slope <- function(data) {

    time <- data$time

    rolling_slope <- function(x) {
      temp <- vector()
      # 9 point window
      for (i in seq(1,length(x)-9+1,1)) {
        re <- lm(col~time,data=data.frame("time"=time[i:(i+9-1)],"col"=x[i:(i+9-1)]))
        re <- as.numeric(re$coefficients[2])
        temp <- c(temp,re)
      }
      #temp <- temp[order(temp)][1:5]
      return(temp)
    }

    slope <- apply(data[,-1], 2, rolling_slope)
    return(slope)
  }

  slope_OD <- as.data.frame(slope(biotek_OD600))*10^2



  # the max slope time point will be the median of the 8 replicates
  max_slope_time_point <- apply(slope_OD,2,which.max)

  median_8 <- function(x) {
    for (i in 1:12) {
      x[(8*i-7):(8*i)] <- median(x[(8*i-7):(8*i)])
    }
    return(round(x,0))
  }


  max_slope_time_point <- median_8(max_slope_time_point) + 4 # middle in 9 points, + 4 is because the calculation of max_slope_time_point


  # save figures
  png(paste0(outputdir,"slope.png"),width=3600, height=1800)

  plot_slope <- function(data,x) {
    plot(rownames(data), data[,8*x-6],type = "o",main = x,col="purple")
    lines(rownames(data),data[,8*x-5],col="orange",type = "o")
    lines(rownames(data),data[,8*x-4],col="red",type = "o")
    lines(rownames(data),data[,8*x-3],col="green",type = "o")
    lines(rownames(data),data[,8*x-2],col="yellow",type = "o")
    lines(rownames(data),data[,8*x-1],col="blue",type = "o")
  }

  plot_slope(slope_OD,3)

  dev.off()






  # calculate interval value for the two doubling time area, with max slope time point as the middle point ( Keren, Zackay et al. 2013. GFP expression per cell per second) ---------


  # calculate the area--------------------------------------------
  integral_value <- function(data,from,to) {

    if ( is.na(from) ) { from <- 1; to <- 10 } # NAN ???

    time <- data$time[from:to]

    # function to calculate single integral value
    single_integral <- function(x) {
      sum(diff(time)*rollmean(x,2))
    }



    # calculate all dataframe
    re <- apply(data[from:to,], 2, single_integral )

    return(re[2])
  }


  # calculate GFP_t2 - GFP_t1 --------------------------------------
  diff_fi <-  function(data,from,to) {
    interval <- function(x) {
      x[length(x)]-x[1]
    }

    if ( is.na(from) ) { from <- 1; to <- 10 } # NAN ???

    re <- interval(data[from:to])
    return(re)
  }


  # use the OD600 data to calculate area (per cell per second)-----------------

  doubling_time <- growth_rate$doubling_time.h
  doubling_time[doubling_time>10] <- 0

  interval_area <- function(data,max_slope_time_point) {
    temp <- vector()
    for(i in 1:length(max_slope_time_point)) {
      re <- integral_value(data[,c(1,(i+1))],
                           round((max_slope_time_point[i]-doubling_time[i]*read_times_per_hour)),
                           round((max_slope_time_point[i]+doubling_time[i]*read_times_per_hour)))
      temp <- c(temp,re)
    }
    return(temp)
  }


  OD600_interval_area <- interval_area(biotek_OD600,max_slope_time_point)

  # use GFP data to calculate GFP_t1 - GFP_t2------------------------

  interval_fi <- function(data,max_slope_time_point) {
    temp <- vector()
    for(i in 1:length(max_slope_time_point)) {
      re <- diff_fi(data[,c((i+1))],
                    round((max_slope_time_point[i]-doubling_time[i]*read_times_per_hour)),
                    round((max_slope_time_point[i]+doubling_time[i]*read_times_per_hour)))
      temp <- c(temp,re)
    }
    return(temp)
  }


  GFP_diff_value <- interval_fi(biotek_GFP,max_slope_time_point)
  mCherry_diff_value <- interval_fi(biotek_mCherry,max_slope_time_point)
  wassabi_diff_value <- interval_fi(biotek_wassabi,max_slope_time_point)



  # let the magnitude at the same level

  OD600_interval_area <- OD600_interval_area/median(OD600_interval_area)*100

  GFP_diff_value <- GFP_diff_value/median(GFP_diff_value)*100

  mCherry_diff_value <- mCherry_diff_value/median(mCherry_diff_value)*100

  wassabi_diff_value <- wassabi_diff_value/median(wassabi_diff_value)*100


  GFP_expression <- GFP_diff_value/OD600_interval_area
  mCherry_expression <- mCherry_diff_value/OD600_interval_area
  wassabi_expression <- wassabi_diff_value/OD600_interval_area


  growth_rate$FI_1 <- GFP_expression

  growth_rate$FI_2 <- mCherry_expression

  growth_rate$FI_3 <- wassabi_expression

  write.xlsx(growth_rate, outputfile2, sheetName=paste(Sys.Date(),"result"))


}












