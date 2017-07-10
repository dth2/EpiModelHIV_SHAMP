
#' @title Demography Check Module
#'
#' @description Module track demographic distribution over time.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
demog_shamp <- function(dat, at) {

 if(at==4){demog4<<-table(dat$attr$demog.cat)}
 if(at==52){demog52<<-table(dat$attr$demog.cat)}
 if(at==104){demog104<<-table(dat$attr$demog.cat)}
 if(at==208){demog208<<-table(dat$attr$demog.cat)}
 if(at==312){demog312<<-table(dat$attr$demog.cat)}
 if(at==416){demog416<<-table(dat$attr$demog.cat)}
 if(at==520){demog520<<-table(dat$attr$demog.cat)} 
 if(at==1040){demog1040<<-table(dat$attr$demog.cat)}
  return(dat)
}
