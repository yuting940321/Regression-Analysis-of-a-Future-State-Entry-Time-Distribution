library(timereg)

#N_{jj'} count the number of state j to state j' transitions in the interval [t, t']. t is the transition time to state j, and t' is the transition time to state j'.
#Y_j is the number of individual who have risk to enter other state
#t01_true is the real transition time from state 0 to state 1; t02_true is the real transition time from state 0 to state 2; 
#t23_true is the real transition time from state 2 to state 3; t24_true is the real transition time from state 2 to state 4. 
#N01, N02 and Y0 is calculate by formula in Datta(2002) (page 3); N32new, N42new, Y02new is calculated by my regression shown in my manuscript.

#Knew is the function to calculate the inverse weighting estimation using the formula in (Datta and Satten 2001). K is the P(C_i>t).
Knew<-Vectorize(function(time,i){
  j<-sum(time>=time1)
  if(j==0){
    return(K[i,1])
  }else{
    return(K[i,j])
  }
})
#Knewpse is the function to calculate the inverse weighting estimation using the formula in (Datta and Satten 2001) in Pseudo-Value approach. K is the P(C_i>t). 
Knewpse<-Vectorize(function(time,i){
  j<-sum(time>=time1pse)
  if(j==0){
    return(Kpse[i,1])
  }else{
    return(Kpse[i,j])
  }
})
#N01 count the number that transit from state 0 to state 1. 
N01<-function(time){
  #a is the transition time to state 1 for individuals who is observed that transition time to state 1 is less than time. 
  a_which<-which(t01_true<=time & !is.na(t01_true) & is.na(t02_true) & censor>=t01_true)
  a<-t01_true[a_which]-1e-10
  #if nobody is observed that the transition time to state 1 is less than time, n01=0
  if (length(which(!is.na(a)))==0){
    n01<-0
  }
  #if some individuals are observed that the transition time to state 1 is less than time, n01 is calculated based on the formula 3.4 in Datta(2002)
  if (length(which(!is.na(a)))>0){
    KK<-Knew(a,a_which)
    n01<-sum(1/KK)
  }
  return(n01)
}
#N01pse count the number that transit from state 0 to state 1 in pseudo-value approach.
N01pse<-function(time){
  a_which<-which(t01_truepse<=time & !is.na(t01_truepse) & is.na(t02_truepse) & censorpse>=t01_truepse)
  a<-t01_truepse[a_which]-1e-10
  if (length(which(!is.na(a)))==0){
    n01<-0
  }
  if (length(which(!is.na(a)))>0){
    # KK<-sapply(a[which(!is.na(a))], Knew)
    KK<-Knewpse(a,a_which)
    n01<-sum(1/KK)
  }
  return(n01)
}

#N02 count the number that transit from state 0 to state 2. 
N02<-function(time){
  #a is the transition time to state 2 for individuals who is observed that transition time to state 2 is less than time.
  a_which<-which(t02_true<=time & !is.na(t02_true) & is.na(t01_true) & censor>=t02_true)
  a<-t02_true[a_which]-1e-10
  #if nobody is observed that the transition time to state 2 is less than time, n02=0
  if (length(which(!is.na(a)))==0){
    n02<-0
  }
  #if some individuals are observed that the transition time to state 2 is less than time, n02 is calculated based on the formula 3.4 in Datta(2002)
  if (length(which(!is.na(a)))>0){
    KK<-Knew(a,a_which)
    n02<-sum(1/KK)
  }
  return(n02)
}
#N02pse count the number that transit from state 0 to state 2 in pseudo value approach. 
N02pse<-function(time){
  a_which<-which(t02_truepse<=time & !is.na(t02_truepse) & is.na(t01_truepse) & censorpse>=t02_truepse)
  a<-t02_truepse[a_which]-1e-10
  if (length(which(!is.na(a)))==0){
    n02<-0
  }
  if (length(which(!is.na(a)))>0){
    # KK<-sapply(a[which(!is.na(a))], Knew)
    # n02<-sum(1/KK)
    KK<-Knewpse(a,a_which)
    n02<-sum(1/KK)
  }
  return(n02)
}
#N32new count the number who is observed to transit from state 2 to state 3
N32new<-function(time){
  ##phi_i2 is the probability of entering state 2 and is calculated by the formula in my manuscript
  #a is probability of entering state 2 for individual who is observed to transit to state 3 
  #b is transition time to state 3 for individual who is observed to transit to state 3 
  which_a<-which(t23_true<=time & delta23==1)
  a<-phi_i2new[which_a]
  b<-t23_true[which_a]
  #if nobody is observed that the transition time to state 3 is less than time, n32=0
  if (length(which(!is.na(b)))==0){
    n32new<-0
  }
  #if some individuals are observed that the transition time to state 3 is less than time, n32 is calculated based on the formula in my manuscript
  if (length(which(!is.na(b)))!=0){
    cc<-Knew(b,which_a)
    n32new<-sum(a/cc)
  }
  return(n32new)
}
#N32newpse count the number who is observed to transit from state 2 to state 3 in pseudo value approach
N32newpse<-function(time){
  which_a<-which(t23_truepse<=time & delta23pse==1)
  a<-phi_i2newpse[which_a]
  b<-t23_truepse[which_a]
  if (length(which(!is.na(b)))==0){
    n32new<-0
  }
  if (length(which(!is.na(b)))!=0){
    # cc<-sapply(b[which(!is.na(b))], Knew)
    cc<-Knewpse(b,which_a)
    n32new<-sum(a/cc)
  }
  return(n32new)
}
#N42new count the number who is observed to transit from state 2 to state 4
N42new<-function(time){
  #a is probability of entering state 2 for individual who is observed to transit to state 4 
  #b is transition time to state 4 for individual who is observed to transit to state 4
  which_a<-which(t24_true<=time & delta24==1)
  a<-phi_i2new[which_a]
  b<-t24_true[which_a]-1e-10
  #if nobody is observed that the transition time to state 4 is less than time, n42=0
  if (length(which(!is.na(b)))==0){
    n42new<-0
  }
  #if some individuals are observed that the transition time to state 4 is less than time, n42 is calculated based on the formula in my manuscript
  if (length(which(!is.na(b)))!=0){
    cc<-Knew(b,which_a)
    n42new<-sum(a/cc)
  }
  return(n42new)
}
#N42newpse count the number who is observed to transit from state 2 to state 4 in pseudo-value approach. 
N42newpse<-function(time){
  which_a<-which(t24_truepse<=time & delta24pse==1)
  a<-phi_i2newpse[which_a]
  b<-t24_truepse[which_a]-1e-10
  if (length(which(!is.na(b)))==0){
    n42new<-0
  }
  if (length(which(!is.na(b)))!=0){
    # cc<-sapply(b[which(!is.na(b))], Knew)
    cc<-Knewpse(b,which_a)
    n42new<-sum(a/cc)
  }
  return(n42new)
}
################################
#Y0 count the number of individual who have risk to enter state 1 or 2
Y0<-function(time){
  #if_a and if_b are the individuals whose transition time of state 1 and 2 are larger than time and censoring time are larger than time
  if_a<-time<=t01_true & censor>=time & !is.na(t01_true)
  if_b<-time<=t02_true & censor>=time & !is.na(t02_true)
  which_ab<-which(if_a|if_b)
  # formula 3.5 in Datta(2002)
  if (length(which_ab)>0){
    y0<-sum(1/Knew(time-1e-10,which_ab))
  }
  if (length(which_ab)==0){
    y0<-0
  }
  return(y0)
}
################################
#Y0pse count the number of individual who have risk to enter state 1 or 2 in pseudo-value approach
Y0pse<-function(time){
  if_a<-time<=t01_truepse & censorpse>=time & !is.na(t01_truepse)
  if_b<-time<=t02_truepse & censorpse>=time & !is.na(t02_truepse)
  which_ab<-which(if_a|if_b)
  if (length(which_ab)>0){
    y0<-sum(1/Knewpse(time-1e-10,which_ab))
  }
  if (length(which_ab)==0){
    y0<-0
  }
  return(y0)
}
################################
#Y02new count the number who have the risk to enter state 3 or 4
Y02new<-function(time){
  #if_a and if_b is indicator of transiting to state 3 or 4
  if_a<-t23_true>time & censor>time & !is.na(t23_true)
  if_b<-t24_true>time & censor>time & !is.na(t24_true)
  #which_ab represents individual who enter state 3 or 4
  which_ab<-which(if_a | if_b)
  #y02 use the formula Y_{0^*|2} in my manuscript 
  if (length(which_ab)>0){
    a<-phi_i2new[which_ab]
    y02<-sum((1/Knew(time-1e-10, which_ab))*a)
  }
  if (length(which_ab)==0){
    y02<-0
  }
  return(y02)
}
#################################
#Y02newpse count the number who have the risk to enter state 3 or 4 in pseudo-value approach
Y02newpse<-function(time){
  if_a<-t23_truepse>time & censorpse>time & !is.na(t23_truepse)
  if_b<-t24_truepse>time & censorpse>time & !is.na(t24_truepse)
  which_ab<-which(if_a | if_b)
  if (length(which_ab)>0){
    a<-phi_i2newpse[which_ab]
    y02<-sum((1/Knewpse(time-1e-10, which_ab))*a)
  }
  if (length(which_ab)==0){
    y02<-0
  }
  return(y02)
}

#state is an time-dependent covariate. It is used in Aalen's model to calculate the inverse weighting estimator K.
sta<-function(time){
  stage<-c()
  for(i in 1:nind){
    if(state1[i]==0 & state2[i]==0 & state3[i]==0 & state4[i]==0){
      stage[i]=ifelse(time<censor[i],0,NA)
    }else if(state1[i]==1 & state2[i]==0 & state3[i]==0 & state4[i]==0){
      stage[i]=ifelse(time<t01_true[i],0,NA)
    }else if(state1[i]==0 & state2[i]==1 & state3[i]==0 & state4[i]==0){
      stage[i]=ifelse(time<t02_true[i],0,ifelse(time<censor[i],1,0))
    }else if(state1[i]==0 & state2[i]==1 & state3[i]==1 & state4[i]==0){
      stage[i]=ifelse(time<t02_true[i],0,ifelse(time<t23_true[i],1,0))
    }else if(state1[i]==0 & state2[i]==1 & state3[i]==0 & state4[i]==1){
      stage[i]=ifelse(time<t02_true[i],0,ifelse(time<t24_true[i],1,0))
    }
  }
  return(stage)
}

#statepse is an time-dependent covariate. It is used in Aalen's model to calculate the inverse weighting estimator K in pseudo-value approach.
stapse<-function(time){
  stagepse<-c()
  for(i in 1:(nind-1)){
    if(state1pse[i]==0 & state2pse[i]==0 & state3pse[i]==0 & state4pse[i]==0){
      stagepse[i]=ifelse(time<censorpse[i],0,NA)
    }else if(state1pse[i]==1 & state2pse[i]==0 & state3pse[i]==0 & state4pse[i]==0){
      stagepse[i]=ifelse(time<t01_truepse[i],0,NA)
    }else if(state1pse[i]==0 & state2pse[i]==1 & state3pse[i]==0 & state4pse[i]==0){
      stagepse[i]=ifelse(time<t02_truepse[i],0,ifelse(time<censorpse[i],1,NA))
    }else if(state1pse[i]==0 & state2pse[i]==1 & state3pse[i]==1 & state4pse[i]==0){
      stagepse[i]=ifelse(time<t02_truepse[i],0,ifelse(time<t23_truepse[i],1,NA))
    }else if(state1pse[i]==0 & state2pse[i]==1 & state3pse[i]==0 & state4pse[i]==1){
      stagepse[i]=ifelse(time<t02_truepse[i],0,ifelse(time<t24_truepse[i],1,NA))
    }
  }
  return(stagepse)
}


#numbe of individual
nind<-50

#result is whether the p-value is smaller than 0,05. If p-value<0.05, result is reject; otherwise, result is accept
result<-rep(NA, 2000)
#censor_s1p is the censoring percentage at state 1
censor_s1p<-rep(NA, 2000)
#censor_s2p is the censoring percentage at state 2
censor_s2p<-rep(NA, 2000)
#hai32new is the Psi_{3|2}
hai32new<-rep(NA, 2000)
#hai32new0 is the Psi_{3|2} as covariate z=0
hai32new0<-rep(NA, 2000)
#hai32new0 is the Psi_{3|2} as covariate z=1
hai32new1<-rep(NA, 2000)
simulation_begin<-1
simulation_end<-2000

for (simulation in simulation_begin:simulation_end){
  set.seed(simulation)
  #alpha and bet is used in logistic regression to calculate the probability of entering state 3 and 4
  alpha=0
  beta=0.0
  #20% individual transit to state 1, 80% individuals transit to state 2
  group1<-rbinom(nind, 1, 0.2)
  #covariate is the baseline binary covariate
  covariate<-gender<-rbinom(nind, 1, 0.5)
  #t00_true is the transition time of state 0
  t00_true<-rep(0, nind)
  #t01_true is the transition time of state 1
  t01_true<-rep(NA, nind)
  #t02_true is the transition time of state 2
  t02_true<-rep(NA,nind)
  #tt is the transition time of state 1 or 2
  tt<-rweibull(nind, 4, 5)
  #if group1 = 1, individual transit to state 1; f group = 0, individual transit to state 2
  t01_true<-ifelse(group1==1, tt, NA)
  t02_true<-ifelse(group1==0, tt, NA)
  #after entering state 2, 50% individual transit to state 3, 50% individuals transit to state 4
  psi<-exp(alpha+beta*covariate)/(1+exp(alpha+beta*covariate))
  #group2 is the probability of entering state 3 and 4
  group2<-rbinom(nind, 1, psi)
  #t23_true is the transition time from state 2 to state 3
  t23_true<-rep(NA, nind)
  #t24_true is the transition time from state 2 to state 4
  t24_true<-rep(NA, nind)
  #tt1 is the transition time if state 3 or 4
  tt1<-rweibull(nind, 6,2)+t02_true
  #if individual entered state 2, as group2=1, they transit to state 3; otherwise, they transit to state 4
  t23_true[which(!is.na(t02_true) & group2==1)]=tt1[which(!is.na(t02_true) & group2==1)]
  t24_true[which(!is.na(t02_true) & group2==0)]=tt1[which(!is.na(t02_true) & group2==0)]
  #dependent-state censoring generation
  cc0<-rweibull(nind,5,7)
  texit0<-ifelse(is.na(t01_true),t02_true,t01_true)
  censor<-ifelse(cc0<=texit0,cc0,qweibull(pweibull(cc0,4,8)+runif(nind)*(1-pweibull(cc0,4,8)),4,8))
  #t1, t2, t3 and t4 are observed transiting time to state 1,2,3 and 4
  t1<-ifelse(t01_true<censor,t01_true,NA)
  t2<-ifelse(t02_true<censor,t02_true,NA)
  t3<-ifelse(t23_true<censor,t23_true,NA)
  t4<-ifelse(t24_true<censor,t24_true,NA)
  ##delta01 is the censoring indicator of state 1
  delta01<-t01_true<censor&!is.na(t01_true)
  #delta02 is the censoring indicator of state 2
  delta02<-t02_true<censor&!is.na(t02_true)
  #delta23 is the censoring indicator of state 3
  delta23<-t23_true<censor&!is.na(t23_true)
  #delta24 is the censoring indicator of state 4
  delta24<-t24_true<censor&!is.na(t24_true)
  #if statej=1, individual are observed to enter state j at time tj
  state1<-ifelse(delta01, 1, 0)
  state2<-ifelse(delta02, 1, 0)
  state3<-ifelse(delta23, 1, 0)
  state4<-ifelse(delta24, 1, 0)
  #tistar is the real last transition time
  tistar<-rep(NA, nind)
  tistar[which(!is.na(t01_true) & is.na(t02_true))]=t01_true[which(!is.na(t01_true) & is.na(t02_true))]
  tistar[which(is.na(t01_true) & !is.na(t02_true) & is.na(t23_true) & !is.na(t24_true))]=
    t24_true[which(is.na(t01_true) & !is.na(t02_true) & is.na(t23_true) & !is.na(t24_true))]
  tistar[which(is.na(t01_true) & !is.na(t02_true) & is.na(t24_true) & !is.na(t23_true))]=
    t23_true[which(is.na(t01_true) & !is.na(t02_true) & is.na(t24_true) & !is.na(t23_true))]
  #ti is the minimum of last transition time and censoring time
  ti<-pmin(tistar,censor)
  #tt3 is observed transition time of state 2
  tt3<-ifelse(censor>t02_true & !is.na(t02_true), t02_true, NA)
  #tinew is the combination of all transition time
  tinew<-sort(unique(c(ti, tt3[which(!is.na(tt3))])))
  #if censor<tistar, individual is censored and indicator censor=1; otherwise, indicator censor=0
  indicator_censor<-ifelse(censor<tistar,1,0)
  #censoring percentage of state 1
  censor_s1<-ifelse(censor<t01_true & !is.na(t01_true), 1, 0)
  censor_s1p[simulation]<-sum(censor_s1)/length(which(!is.na(t01_true)))
  #censoring percentage of state 2
  censor_s2<-ifelse(censor<t02_true & !is.na(t02_true), 1, 0)
  censor_s2p[simulation]<-sum(censor_s2)/length(which(!is.na(t02_true)))
  #reorganize dataset for Aalen's model. If invidual get censored at state 0 or 2, next state is "C"
  indi<-seq(1, nind, 1)
  data<-data.frame(indi, state1, state2, state3, state4, t00_true, t01_true, t02_true,t23_true,
                   t24_true, censor, tistar, ti, delta01, delta02, delta23, delta24, t1, t2, t3, t4)
  datalong<-data.frame()
  a<-which(state1==0 & state2==0 & state3==0 & state4==0)
  if (length(a)>0){
    start<-data[a,"t00_true"]
    end<-data[a,"ti"]
    newrow1<-data.frame(indi=data$indi[a],start=start,end=end,cov=covariate[a],next_state="C",stage=0)
    datalong<-rbind(datalong,newrow1)
  }
  b<-which(state1==1 & state2==0 & state3==0 & state4==0)
  if (length(b)>0){
    start<-data[b,"t00_true"]
    end<-data[b,"ti"]
    newrow2<-data.frame(indi=data$indi[b],start=start,end=end,cov=covariate[b],next_state="1",stage=0)
    datalong<-rbind(datalong,newrow2)
  }
  c<-which(state1==0 & state2==1 & state3==0 & state4==0)
  if (length(c)>0){
    start<-data[c,"t00_true"]
    end<-data[c,"t2"]
    newrow3<-data.frame(indi=data$indi[c],start=start,end=end,cov=covariate[c],next_state="2",stage=0)
    start<-data[c,"t2"]
    end<-data[c,"ti"]
    newrow4<-data.frame(indi=data$indi[c],start=start,end=end,cov=covariate[c],next_state="C",stage=2)
    datalong<-rbind(datalong,newrow3, newrow4)
  }
  d<-which(state1==0 & state2==1 & state3==1 & state4==0)
  if (length(d)>0){
    start<-data[d,"t00_true"]
    end<-data[d,"t2"]
    newrow5<-data.frame(indi=data$indi[d],start=start,end=end,cov=covariate[d],next_state="2",stage=0)
    start<-data[d,"t2"]
    end<-data[d,"ti"]
    newrow6<-data.frame(indi=data$indi[d],start=start,end=end,cov=covariate[d],next_state="3",stage=2)
    datalong<-rbind(datalong,newrow5, newrow6)
  }
  e<-which(state1==0 & state2==1 & state3==0 & state4==1)
  if (length(e)>0){
    start<-data[e,"t00_true"]
    end<-data[e,"t2"]
    newrow7<-data.frame(indi=data$indi[e],start=start,end=end,cov=covariate[e],next_state="2",stage=0)
    start<-data[e,"t2"]
    end<-data[e,"ti"]
    newrow8<-data.frame(indi=data$indi[e],start=start,end=end,cov=covariate[e],next_state="4",stage=2)
    datalong<-rbind(datalong,newrow7, newrow8)
  }
  #stage is the covariate used in Aalen's model since state affects the censoring.
  datalong$stage<-factor(datalong$stage,c(0,2))
  #Aalen's model is used to calculate the cumulative hazard function
  #We would calculate the inverse weighting estimator K using the formula in (Datta and Satten 2001). K is the P(C_i>t).
  fit2<-aalen(Surv(time=start,
                   time2=end,
                   event=next_state=="C")~stage,data=datalong)
  #cum is cumulative hazard function
  cum<-fit2$cum
  #time1 is the time that individual get censored
  time1<-cum[,"time"]
  #nnn is the number of indicator who get censored 
  nnn<-nrow(fit2$cum)
  if (nnn==1){
    K<-matrix(1, nrow=nind, ncol=nnn)
  }
  if (nnn>1){
    #beta0 and beta1 are the parameter used in hazard function
    beta0<-c(0,fit2$cum[2:nnn,"(Intercept)"]-fit2$cum[1:(nnn-1),"(Intercept)"]) ###
    beta1<-c(0,fit2$cum[2:nnn,"stage2"]-fit2$cum[1:(nnn-1),"stage2"]) ###
    Bhat0<-fit2$cum[1:nnn,"(Intercept)"]
    Bhat1<-fit2$cum[1:nnn,"stage2"]
    lambda<-matrix(NA, nrow = nind, ncol = nnn)
    cova<-matrix(NA, nind, nrow(cum))
    for (i in 1:nrow(cum)){
      time<-fit2$cum[i,"time"]
      cova[, i]<-sta(time)
    }
    #blambda is the hazard function
    blambda<-matrix(NA, nrow = nind, ncol = nrow(cum))
    for (i in 1:nind){
      blambda[i, ]<-beta0+cova[i, ] * beta1
    }
    blambda[is.na(blambda)]<-0
    #K is the inverse weighting estimator
    K<-matrix(NA, nrow=nind, ncol=nnn)
    K[,1]<-exp(-blambda[,1])
    for (j in 2:nnn){
      K[,j]<-exp(-blambda[,j])*K[,j-1]
    }
    K[abs(K)<1e-16]<-0
  }
  #n01 count the number of individual entering state 1 at tinew
  n01<-sapply(tinew, N01)
  #n02 count the number of individual entering state 2 at tinew
  n02<-sapply(tinew, N02)
  #n01 count the number of individual entering state 1 at tinew-
  n01_<-sapply(tinew-1e-10,N01)
  #n02 count the number of individual entering state 2 at tinew-
  n02_<-sapply(tinew-1e-10,N02)
  #y0 count the number of individual who have risk of entering state 1 and 2 at tinew
  y0<-sapply(tinew, Y0)
  #deltan01 and deltan02 are the numbers of individual entering state 1 and state 2 in time [tinew[tau_k],tinew[tau_{k+1}]]
  deltan01<-n01-n01_
  deltan02<-n02-n02_
  #Estimate the phi_i2 by formula in section 2.1 of my manuscript
  survival00<-ifelse(y0==0,1,1-(deltan01+deltan02)/y0)
  transition02<-ifelse(y0==0,0,deltan02/y0)
  
  phi_i2new<-rep(NA,nind)
  phi_i2new[which(delta02)]<-1
  phi_i2new[which(delta01)]<-0

  for(i in which(is.na(phi_i2new))){
    idx<-censor[i]<tinew
    aprod<-cumprod(c(1,head(survival00[idx],-1)))
    phi_i2new[i]<-sum(aprod*transition02[idx])
  }
  #n32new count the number of individual entering state 3 at tinew
  n32new<-sapply(tinew, N32new)
  #n42new count the number of individual entering state 4 at tinew
  n42new<-sapply(tinew, N42new)
  #y02new count the number who have the risk to enter state 3 or 4 at tinew 
  y02new<-sapply(tinew, Y02new)
  #n32new_ count the number of individual entering state 3 at tinew-
  n32new_<-sapply(tinew-1e-10, N32new)
  #y02new_ count the number who have the risk to enter state 3 or 4 at tinew- 
  y02new_<-sapply(tinew-1e-10, Y02new)
  #n42new_ count the number of individual entering state 3 at tinew-
  n42new_<-sapply(tinew-1e-10, N42new)
  #deltan32new and deltan42new are the number of individual entering state 3 and state 4
  deltan32new<-n32new-n32new_
  deltan42new<-n42new-n42new_
  #survival02 is the S_{0^*|2} in my manuscript
  survival02<-ifelse(y02new==0,0,(1-(deltan32new+deltan42new)/y02new))
  Survival02new_<-cumprod(c(1,head(survival02,-1)))
  #hai32 is the estimation of \psi_{3|2}
  hai32<-ifelse(y02new==0,0,Survival02new_*(deltan32new)/y02new)
  hai32new[simulation]<-sum(hai32)
  #pseudo value approach result
  pseu_hai32<-c()
  #remove ith individual from data and estimate the \psi_{3|2}
  for (i in 1:nind){
    group1pse<-group1[-i]
    t00_truepse<-t00_true[-i]
    t01_truepse<-t01_true[-i]
    t02_truepse<-t02_true[-i]
    t23_truepse<-t23_true[-i]
    t24_truepse<-t24_true[-i]
    censorpse<-censor[-i]
    t1pse<-t1[-i]
    t2pse<-t2[-i]
    t3pse<-t3[-i]
    t4pse<-t4[-i]
    covariatepse<-genderpse<-gender[-i]
    delta01pse<-delta01[-i]
    delta02pse<-delta02[-i]
    delta23pse<-delta23[-i]
    delta24pse<-delta24[-i]
    state1pse<-state1[-i]
    state2pse<-state2[-i]
    state3pse<-state3[-i]
    state4pse<-state4[-i]
    tistarpse<-tistar[-i]
    tipse<-ti[-i]
    datalongpse<-datalong[datalong$indi!=i, ]
    fit2pse<-aalen(Surv(time=start,
                        time2=end,
                        event=next_state=="C")~stage,data=datalongpse)
    cumpse<-fit2pse$cum
    time1pse<-cumpse[,"time"]
    nnnpse<-nrow(fit2pse$cum)
    if (nnnpse>1){
      beta0pse<-c(0,fit2pse$cum[2:nnnpse,"(Intercept)"]-fit2pse$cum[1:(nnnpse-1),"(Intercept)"]) ###
      beta1pse<-c(0,fit2pse$cum[2:nnnpse,"stage2"]-fit2pse$cum[1:(nnnpse-1),"stage2"]) ###
      ### I added the following things
      Bhat0pse<-fit2pse$cum[1:nnnpse,"(Intercept)"]
      Bhat1pse<-fit2pse$cum[1:nnnpse,"stage2"]
      covapse<-matrix(NA, nind-1, nrow(cumpse))
      for (j in 1:nrow(cumpse)){
        time<-fit2pse$cum[j,"time"]
        covapse[, j]<-stapse(time)
      }
      blambdapse<-matrix(NA, nrow = nind-1, ncol = nrow(cumpse))
      for (j in 1:(nind-1)){
        blambdapse[j, ]<-beta0pse+covapse[j, ] * beta1pse
      }
      blambdapse[is.na(blambdapse)]<-0
      Kpse<-matrix(NA, nrow=nind-1, ncol=nnnpse)
      Kpse[,1]<-exp(-blambdapse[,1])
      for (j in 2:nnnpse){
        Kpse[,j]<-exp(-blambdapse[,j])*Kpse[,j-1]
      }
      Kpse[abs(Kpse)<1e-16]<-0
    }
    if (nnnpse<=1){
      Kpse<-matrix(1, nrow=nind-1, ncol=nnnpse)
    }
    n01pse<-sapply(tinew,N01pse)
    n02pse<-sapply(tinew,N02pse)
    n01pse_<-sapply(tinew-1e-10,N01pse)
    n02pse_<-sapply(tinew-1e-10,N02pse)
    y0pse<-sapply(tinew, Y0pse)

    deltan01pse<-n01pse-n01pse_
    deltan02pse<-n02pse-n02pse_
    survival00pse<-ifelse(y0pse==0,1,1-(deltan01pse+deltan02pse)/y0pse)
    transition02pse<-ifelse(y0pse==0,0,deltan02pse/y0pse)
    
    phi_i2newpse<-rep(NA,nind-1)
    phi_i2newpse[which(delta02pse)]<-1
    phi_i2newpse[which(delta01pse)]<-0
    
    for(jj in which(is.na(phi_i2newpse))){
      idx<-censorpse[jj]<tinew
      aprod<-cumprod(c(1,head(survival00pse[idx],-1)))
      phi_i2newpse[jj]<-sum(aprod*transition02pse[idx])
    }

    
    n32newpse<-sapply(tinew, N32newpse)
    n42newpse<-sapply(tinew, N42newpse)
    y02newpse<-sapply(tinew, Y02newpse)
    n32newpse_<-sapply(tinew-1e-10, N32newpse)
    y02newpse_<-sapply(tinew-1e-10, Y02newpse)
    n42newpse_<-sapply(tinew-1e-10, N42newpse)
    deltan32newpse<-n32newpse-n32newpse_
    deltan42newpse<-n42newpse-n42newpse_
    
    survival02pse<-ifelse(y02newpse==0,0,(1-(deltan32newpse+deltan42newpse)/y02newpse))
    Survival02newpse_<-cumprod(c(1,head(survival02pse,-1)))
    hai32pse<-ifelse(y02newpse==0,0,Survival02newpse_*(deltan32newpse)/y02newpse)
    hai32newpse<-sum(hai32pse)
    
    pseu_hai32[i]=nind*hai32new[simulation]-(nind-1)*hai32newpse
  }
  #revert the pseudo value to either or 1 by the method in section 3.5
  response<-ifelse(pseu_hai32<0.5, 0, 1)
  #logistic regression model to test the significance of binary covariate
  mod<-glm(response~gender, family = "binomial")
  asumary<-summary(mod)
  #p-value of covariate in the model 
  pvalglmpse<-asumary$coefficients["gender","Pr(>|z|)"]
  result[simulation]<-ifelse(pvalglmpse<0.05, "reject", "accept")
  #alpha1 and beta1 are estimated parater from logistic regression
  alpha1<-asumary$coefficients["(Intercept)","Estimate"]
  beta1<-asumary$coefficients["gender","Estimate"]
  #estimating \psi_{3|2}(z) by the logistic regression
  hai32new0[simulation]<-exp(alpha1+beta1*0)/(1+exp(alpha1+beta1*0))
  hai32new1[simulation]<-exp(alpha1+beta1*1)/(1+exp(alpha1+beta1*1))
  print(simulation)
  save.image("SC.RData")
}

#power
(power<-sum(result=="reject")/2000)
#censoring percentage at state 1
(censor_probability1<-mean(censor_s1p))
#censoring percentage at state 2
(censor_probability2<-mean(censor_s2p))
#psi_{3|2}
(Psi32<-mean(hai32new))
#bias of psi_{3|2}
(bias<-mean(hai32new-0.5))
#psi_{3|2}(z)
(psi320<-mean(hai32new0))
(psi321<-mean(hai32new1))
#bias of psi_{3|2}(z)
(bias0<-mean(psi320-0.5))
(bias1<-mean(psi321-0.5))
#standard error of psi_{3|2}(z)
(se<-sqrt(mean((hai32new-Psi32)^2)))
(se0<-sqrt(mean((hai32new0[which(!is.na(hai32new0))]-mean(hai32new0[which(!is.na(hai32new0))]))^2)))
(se1<-sqrt(mean((hai32new1[which(!is.na(hai32new1))]-mean(hai32new1[which(!is.na(hai32new1))]))^2)))


save.image("SC.RData")
