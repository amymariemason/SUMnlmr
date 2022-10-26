#' Gelman--Rubin (GR) Statistics
#' @description  getGRvalues calculates the the Gelman--Rubin (GR) statistics
#' for each strata (i.e. LACE-strata) as described in Haodong's paper on
#' doubly-ranked methods in Text S1
#'
#' @param X vector of coarsened exposures
#' @param Zstratum vector of pre-strata index column
#' @param NC number of chains - a lower NC is encouraged
#' @param tl threshold level, should be >1
#' @author Haodong Tian for underlying code, Amy Mason for conversion to R function <haodong.tian@mrc-bsu.cam.ac.uk>
#' @importFrom dplyr mutate group_by row_number arrange summarise n first
#' @import assertthat
#' @return max_GR is the maximum G-R statistic accross all strata
#' @return statified_strata is a list of the strata where the GR-statistic meets
#'  the threshold limit
#' @return GR_values is a table of the {a} and {b} values for each strata
#' @examples
#'  #Toy data example------------------------------------------------------------------
#'  N<-10000  #N is the total ssample size
#' Z<-rnorm( N,0,0.5   )  # instrument
#' U1<-rnorm(N,0,1)
#' U2<-rnorm(N,0,1)   #U2 is confounder

#' h_X<-function(Z){ 0.5*Z  }
#' X<-h_X(Z)+U2+U1

#' fitX<-lm( X~Z   )
#' RR<-as.numeric(summary(fitX)$r.squared[1])

#' No<-1000   #No is the number of pre-strata (i.e. IV-strata)
#' X_<-X  # true exact X
#' X<-round(X)  #coarsened by round
#' dat<-cbind(   Z, X, X_)  # X: coarsened exposure; X_: exact exposure
#' dat_order<-dat[ order(dat[,1]  ),  ]  #ordered by Z
#' dat_order<-cbind(  dat_order   ,     sort(rep( 1:No, N/No  )))
#' dat_order<-as.data.frame( dat_order )
#' names(dat_order)<-c(  'Z_order', 'X','X_' ,'ZStratum')
#' getGRvalues( X=dat_order$X, Zstratum = dat_order$ZStratum, NC=2, tl = 1.02  )
#' #--------------------------------------------------------------------------------------
#' @export
getGRvalues<-function(X=X,
                      Zstratum=Zstratum,
                      NC=2,
                      tl=1.02
){
  # define local vars
  strata1 <- x0q <- chains <- Lows<- Ups<- NULL
  cmL <- cmU <- wcvL <- wcvU<- tmU <-tmL<- NULL
  cmL <- cmU <- wcvL <- wcvU<- tmU <-tmL<- NULL
  WL <- WU <- BL <- BU<- gr_low <- gr_up<- NULL


  #checks

  stopifnot(
    "NC is not numeric" = is.numeric(NC),
    "tl is not numeric" = is.numeric(tl),
    "X is not numeric" = is.numeric(X),
    "ZStratum is not numeric" = is.numeric(Zstratum),
    "difference number of observations for the strata & exposure" =
      length(X) == length(Zstratum),
    "missing X values" = !anyNA(X),
    "missing Zstratum values" = !anyNA(Zstratum),
    "Inappropriate threshold value, should be >1"= (tl>1)
  )

  # internal values
   N<- length(X) #total sample size
  No<- max(Zstratum)  #number of pre-strata (i.e. IV-strata)
  strata_no<- ceiling(N/No) # number of final strata

  ###sort the data frame such that in each pre-strata the (coarsened) exposure
  #is ordered, first by the pre-strata and then by value of X within the strata

  id<- seq(X)
  temp<- data.frame(x=X,strata1=Zstratum,id=id)
  temp<-group_by(.data=temp, strata1)
  temp<- arrange(.data=temp, x, .by_group = TRUE)
  temp<-mutate(.data=temp, x0q= row_number())

  #obtain the two coefficients ({a} and {b}) (see Equation (A4)&(A5) of Text S1)
  temp<-mutate(.data=temp, Lows=lowfun(x))
  temp<-mutate(.data=temp, Ups=upfun(x))


  ###Gelman--Rubin (GR) statistics----------------------------------------------------------
  GR_low<-c()
  GR_up<-c()
    #split into chains
    temp2<-group_by(.data=temp, x0q)
    temp2<-mutate(.data=temp2, chains = floor((row_number()-1)/n()*NC)+1)
    temp2<-group_by(.data=temp2, chains, .add=TRUE) # creates equiv of cc
    # create summary of chain
    chain_summ<-summarise(.data=temp2,
                          .groups="keep",
                          cmL=mean(Lows),
                          wcvL=var(Lows),
                          cmU=mean(Ups),
                          wcvU=var(Ups))
    #cm is mean of each chain
    #wcv is within-chain variance
    chain_summ<-group_by(.data=chain_summ, x0q)
    chain_summ<-mutate(.data=chain_summ,
                       tmL=mean(cmL),
                       WL=mean(wcvL),
                       tmU=mean(cmU),
                       WU=mean(wcvU))
    #tm is total mean
    # B is between-chain variance
    chain_summ<-mutate(.data=chain_summ,
                       BL = No/NC*(1/(NC-1))*sum(  ( cmL-tmL )^2),
                      BU = No/NC*(1/(NC-1))*sum(  ( cmU-tmU )^2),
    # gr_low, gr_up
                      gr_low=(   ( (No/NC)-1 )/(  No/NC )*WL +
                                   BL/(  (No/NC) )        )/WL,
                      gr_up=(   ( (No/NC)-1 )/(  No/NC )*WU +
                                    BU/(  (No/NC) )        )/WU)


    #store GR values
    GR_low<-summarise(.data=chain_summ, gr_low=first(gr_low))$gr_low
    GR_up<-summarise(.data=chain_summ, gr_up=first(gr_up))$gr_up

  GR_low[is.na( GR_low  )  ]<-1
  #There are at least two constant-value chains,
  #so let their GR values to be 1
  GR_up[is.na( GR_up  )  ]<-1

#evaluate strata
    grlevel<-tl^2  #threshold level
    PP<-(1:(N/No))[  (GR_low<grlevel)&( GR_up<grlevel  ) ]
    #The strata (LACE-strata) passing the GR diagnosis
    max_GR=max(sqrt(GR_up), sqrt(GR_low))

# return warning if issue
    if (max_GR>tl){
      warning("Warning: Gelman-Rubin statistic is above 1.02 for at least one stratum.
      Please consider using fewer strata.")
    }

  return( list( 'max_GR'=max_GR,
                'satisfied_strata'=PP,
                'GR_values' = data.frame( Strata=1:ceiling((N/No))  ,
                                          a_GR=sqrt(GR_up),
                                          b_GR=sqrt(GR_low)  )

  ) )
}

#' Lower Coefficient value
#' @description function to obtain the lower coefficient values
#' (i.e. {b} in Text S1) for each pre-strata
#' @param vec input vector
#' @importFrom utils head
#' @export
#'
lowfun<-function( vec  ){
  rep(c(1,head(cumsum(as.numeric(table(vec)))+1,-1)),
      as.numeric(table(vec)))
}

#' Upper Coefficient value
#' @description function to obtain the upper coefficient values (i.e. {a} in Text S1) for
#' each pre-strata
#' @param vec input vector
#' @importFrom utils head
#' @export
upfun<-function( vec  ){
  rep(cumsum(as.numeric(table(vec))),
      as.numeric(table(vec)))
}





