
# SHAMP  -----------------------------------------------------------------

#' @title Import and Check Data
#'
#' @description Imports the data for simulations using ergm.ego and SHAMP function.  Ego and alter data are checked to verify that all 
#'              of the required fields are present.  Egos and alters with missing data are dropped. 
#'
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{data.params}, which can be passed to
#' EpiModel function \code{params}
#' EpiModel function \code{setup}
#' .
#'
#' @keywords shamp
#'
#' @export
input_shamp <- function(data, data.params, immigration, msm.msmf) {
 

data.params$rates_m <-NULL 
data.params$durs_m <-NULL
data.params$rates_c <-NULL 
data.params$durs_c <-NULL
data.params$msm.frac <-NULL
data.params$msmf.frac <-NULL
  
##Check to make sure there is data for egos and all three alters in the expected data.frames.
names(data) 
x<- names(data)  
  if (any(x=="egos")==FALSE){
    stop("data list must include observed egos data.frame 'egos'",
         call. = FALSE)}

if (any(x=="altersMain")==FALSE){
  stop("data must include alters  data.frame 'altersMain'",
       call. = FALSE)}

if (any(x=="altersCasual")==FALSE){
  stop("data must include alters  data.frame 'altersCasual'",
       call. = FALSE)}

if (any(x=="altersOT")==FALSE){
  stop("data must include alters  data.frame 'egos'",
       call. = FALSE)}

##Check to make sure we have all of the expected attributes for egos.

fields<-c("weight","ego","sex","age","sqrt.age","sex.ident","immigrant","role.class","race","deg.main", "deg.pers")
fields.egos<-names(data$egos)

for (i in 1:length(fields)){
  if(any(fields.egos==fields[i]) == FALSE){
  stop("missing required field in egos",
       call. = FALSE)}}


  ##Check to make sure we have all of the expected attributes for altersMain, altersCasual, and alterOT
  
  fields<-c("ego","len")

  fields.altersMain<-names(data$altersMain)
  for (i in 1:length(fields)){
    if(any(fields.altersMain==fields[i]) == FALSE){
      stop("missing required field in altersMain",
           call. = FALSE)}}
    
  fields.altersCasual<-names(data$altersCasual)
  for (i in 1:length(fields)){
    if(any(fields.altersCasual==fields[i]) == FALSE){
      stop("missing required field in altersCasual",
           call. = FALSE)}}
    
  fields.altersOT<-names(data$altersOT)
  for (i in 1:length(fields)){
    if(any(fields.altersOT==fields[i]) == FALSE){
      stop("missing required field in altersOT'",
           call. = FALSE)}}
    
  
##Drop egos and alters with missing data.
  
##egos with missing data.
  
  fields<-c("weight","ego","sex","age","sqrt.age","sex.ident","immigrant","role.class","race","deg.main", "deg.pers")
  
  ids.list<-NULL
  for(i in 1:length(fields)){
  temp<-is.na(data$egos[i])
  ids<-which(temp==TRUE)
  ego.ids<-data$egos$ego[ids]
  ids.list<-c(ids.list,ego.ids)}
  
  if(length(ego.ids>0)){warning("Egos and associated alters deleted: missing required field",call. = FALSE)
    print(length(ids.list))}

  data$egos<-data$egos[data$egos$ego %in% ids.list == FALSE,]
  data$altersMain<-data$altersMain[data$altersMain$ego %in% ids.list == FALSE,]
  data$altersCasual<-data$altersCasual[data$altersCasual$ego %in% ids.list == FALSE,]
  data$altersOT<-data$altersOT[data$altersOT$ego %in% ids.list == FALSE,]

##Alters with missing data.
##Main
fields<-names(data$altersMain)

alter.list<-NULL
for(i in 1:length(fields)){
  temp<-is.na(data$altersMain[i])
  ids<-which(temp==TRUE)
  alter.ids<-data$altersMain$ego[ids]
  alter.list<-c(alter.list,alter.ids)}

if(length(alter.list>0)){warning("Main alters deleted: missing values",call. = FALSE)
    print(length(alter.list))}
    data$altersMain<-data$altersMain[data$altersMain$ego %in% alter.list == FALSE,]


##Casual.

fields<-names(data$altersCasual)

alter.list<-NULL
for(i in 1:length(fields)){
  temp<-is.na(data$altersCasual[i])
  ids<-which(temp==TRUE)
  alter.ids<-data$altersCasual$ego[ids]
  alter.list<-c(alter.list,alter.ids)}

if(length(alter.list>0)){warning("Casual alters deleted: missing values",call. = FALSE)
   print(length(alter.list))}
    data$altersCasual<-data$altersCasual[data$altersCasual$ego %in% alter.list == FALSE,]



##one time.

fields<-names(data$altersOT)

alter.list<-NULL
for(i in 1:length(fields)){
  temp<-is.na(data$altersOT[i])
  ids<-which(temp==TRUE)
  alter.ids<-data$altersOT$ego[ids]
  alter.list<-c(alter.list,alter.ids)}

if(length(alter.list>0)){warning("OT alters deleted: missing values",call. = FALSE)
  print(length(alter.list))}
    data$altersOT<-data$altersOT[data$altersOT$ego %in% alter.list == FALSE,]
    
##Data checks completed############################################

###  Calculate the paramters needed for estimation that are not line data.
  
  #Limit Main and casual to a 2 by 3 set (Main 0:1, Casual 0:2)  
  nsfg.egodata$egos$deg.pers<-ifelse(nsfg.egodata$egos$deg.pers > 1, 2, nsfg.egodata$egos$deg.pers)
  nsfg.egodata$egos$deg.main<-ifelse(nsfg.egodata$egos$deg.main > 0, 1, nsfg.egodata$egos$deg.main)
    
  
  # Mean durations
  data.params$durs_m <- (min(data$altersMain$len)+(2*(median(data$altersMain$len)))+max(data$altersMain$len))/4
  data.params$durs_m <- (data.params$durs_m/12)*365
  data.params$durs_c <- (min(data$altersCasual$len)+(2*(median(data$altersCasual$len)))+max(data$altersCasual$len))/4
  data.params$durs_c <- (data.params$durs_c/12)*365 
  
  ##If using immigration set Black and Hispanic immigrants to BI and HI.
  if(immigration==TRUE){
    data$egos$race<-ifelse(data$egos$race=="B" & data$egos$immigrant=="Yes","BI",
                                   ifelse(data$egos$race=="H" & data$egos$immigrant=="Yes","HI",
                                                 data$egos$race))
    
    data$altersMain$race<-ifelse(data$altersMain$race=="B" & data$altersMain$immigrant=="Yes","BI",
                                         ifelse(data$altersMain$race=="H" & data$altersMain$immigrant=="Yes","HI",
                                                       data$altersMain$race))
    
    data$altersCasual$race<-ifelse(data$altersCasual$race=="B" & data$altersCasual$immigrant=="Yes","BI",
                                           ifelse(data$altersCasual$race=="H" & data$altersCasual$immigrant=="Yes","HI",
                                                         data$altersCasual$race))
    
    data$altersOT$race<-ifelse(data$altersOT$race=="B" & data$altersOT$immigrant=="Yes","BI",
                                       ifelse(data$altersOT$race=="H" & data$altersOT$immigrant=="Yes","HI",
                                                     data$altersOT$race))
  }

  
  ##If using immigration set Black and Hispanic immigrants to BI and HI.
  if(msm.msmf==TRUE){

    data.params$msm.frac<-max(0,sum(data$egos$sex.ident=="msm")/sum(data$egos$sex=="M"))
    data.params$msmf.frac<-max(0,sum(data$egos$sex.ident=="msmf")/sum(data$egos$sex=="M"))
  }
  

  ##Make vector with age, race and sex specific proportions for demographic consistency in birth and deaths.
  
  data$egos$demog.cat<-rep(NA,length(data$egos$ego))
  
  sex.groups<-sort(unique(data$egos$sex))
  for (i in 1:(length(sex.groups))){
    data$egos$demog.cat<-ifelse(data$egos$sex==sex.groups[i],i*1000,data$egos$demog.cat)      
  }
  
  race.groups<-sort(unique(data$egos$race))
  for (i in 1:(length(race.groups))){
    data$egos$demog.cat<-ifelse(data$egos$race==race.groups[i],data$egos$demog.cat+(i*100),data$egos$demog.cat)      
  }
  
  
  age.groups<-sort(unique(data$egos$age))
  for (i in 1:(length(age.groups))){
    data$egos$demog.cat<-ifelse (data$egos$age==age.groups[i],data$egos$demog.cat+(age.groups[i]),data$egos$demog.cat)      
  }
  

  temp<-table(data$egos$demog.cat)
  data.params$demog.list<-as.numeric(names(temp))
  data.params$demog.dist<-as.numeric(temp)
  data.params$demog.dist<-data.params$demog.dist/(sum(data.params$demog.dist))
  
  #class(data.params) <- "data.params"
 
  return(list(data.params,data))

}


