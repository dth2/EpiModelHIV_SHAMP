
#' @title Size test for memory allocation.
#'
#' @description This module measures the size of the elements of dat for diagnostic purposes.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This only checks object size in dat and not specific memory allocation functions.
#'
#' @return
#' This function returns the \code{dat} with no changes.
#'
#' @keywords diagnostics
#' @export
#'
size.test <- function(dat, at){

  if(at == 2){
  object.param.size <<-rep(NA,520)  
  object.init.size  <<-rep(NA,520)   
  object.control.size <<-rep(NA,520)
  object.attr.size <<-rep(NA,520)
  object.stats.size <<-rep(NA,520)
  object.temp.size <<-rep(NA,520)
  object.epi.size <<-rep(NA,520)
  object.el.size <<-rep(NA,520)
  object.p.size <<-rep(NA,520)
  object.nwparam.size <<-rep(NA,520)
  main.edges <<-rep(NA,520)
  casual.edges <<-rep(NA,520)
  onetime.edges <<-rep(NA,520)
  }
  
  if(at > 2){
    
    object.param.size[at] <<-object.size(dat$param)  
    object.init.size[at] <<-object.size(dat$init)  
    object.control.size[at] <<-object.size(dat$control)
    object.attr.size[at] <<-object.size(dat$attr)
    object.stats.size[at] <<-object.size(dat$stats)
    object.temp.size[at] <<-object.size(dat$temp)
    object.epi.size[at] <<-object.size(dat$epi)
    object.el.size[at] <<-object.size(dat$el)
    object.p.size[at] <<-object.size(dat$p)
    object.nwparam.size[at] <<-object.size(dat$nwparam)
    
    main.edges[at] <<-sum(get_degree(dat$el[[1]]))/2
    casual.edges[at] <<-sum(get_degree(dat$el[[2]]))/2
    onetime.edges[at] <<-sum(get_degree(dat$el[[3]]))/2
}
  
  return(dat)
}


