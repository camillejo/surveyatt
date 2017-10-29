### Code for simulations for paper 1

#Contents:
# 1. phase.sim
# 2. glmm

#### phase.sim is a function for simulating data with 0-K change-points
phase.sim<-function(n,q,drop.rate,change.point=NULL,s=10000,seed=1234,...){
  #n is sample size
  #q is number of questions
  #change.point is vector of change points
  nphase<-length(change.point)+1 #number of phases
  #s is number of simulations
  #drop.rate is a vector of values for attrition
    #the length of number of phases
  if(length(drop.rate)!=nphase)
    stop("Wrong number of attrition amounts")
  #if(any(drop.rate)<0 | any(drop.rate)>1)
   # stop("Attrition values must be between 0 and 1")

  #set attrition rate for each question
  att.patt<-NULL
  if(nphase==1){
    #if only one phase, set constant attrition rate
    att.patt<-rep(drop.rate,q)
  }
  if(nphase>1){
    #if more than one phase
    #set the last phase
    att.patt[(change.point[length(change.point)]):q]<-drop.rate[length(drop.rate)]
    #set remaining phases
    if(nphase>2){
      p<-nphase-1
      while(p>1){
        l<-p-1
        att.patt[change.point[l]:(change.point[p]-1)]<-drop.rate[p]
        p<-p-1
      }
    }
    #set the first phase
    att.patt[1:(change.point[1]-1)]<-drop.rate[1]
  }

  sim.data <- NULL #to hold output
  set.seed(seed)
  ndrop<-matrix(0, nrow=s, ncol=q) #matrix of drops
  drops<-vector() #overall drop counterx

  #run simulation
  for(i in 1:s){
    drop<-0 #dropout counter
    out.all<-NULL
    #creating random environment
    t<-n*q
    ran<-runif(t)
    assign<-matrix(ran,nrow=n,ncol=q)

    #simulate dropout for each participant
    for(j in 1:n){
      out<-cbind(rep(j,q),c(1:q),rep(0,q))
      for(k in 2:q){
        drop.here<-0
        l=k-1
        if(out[l,3]==1) out[k,3]<-1
        if(out[k,3] != 1 && assign[j,k] <= att.patt[k]) {
          out[k,3]<-1
          drop.here<-1
        }
        ndrop[i,k]<-ndrop[i,k]+drop.here
      }
      if(sum(out[,3])>0) drop<-drop+1
      out.all<-rbind(out.all, out)
    }
    outdta<-data.frame(out.all)
    colnames(outdta)<-c("ID", "Q", "Drop")
    sim.data[[i]]<-outdta
    drops[i]<-drop
  }
  #sim.data is all data, drops is number of drops per simulation, ndrop is number
    #of people dropping at each question per simulation
  return(list(sim.data=sim.data,drops=drops, ndrop=ndrop))
}


#### glmm applies glmm and returns vector of pvalues, where first and last drop,
   # number of change-points, and whether there are no change-points
glmm.qbyqdrop<-function(dta, q, alpha=0.05,...){
  #load packages needed
  library(lme4) #GLMM package
  library(MASS) #contrasts

  #prepping the data
  #out<-data.frame(dta)
  #colnames(out)<-c("ID", "Q", "YorN")
  dta$ID<-factor(dta$ID)
  dta$Q<-factor(dta$Q)
  dta$Drop<-factor(dta$Drop)

  #Apply GLMM to entire survey
  contrasts(dta$Q)<-contr.sdif(q) #contrast
  m<-glmer(Drop~(Q)+(1|ID), data=dta, family=binomial, control=glmerControl(optimizer="bobyqa"), nAGQ=0)
  p<-summary(m)$coefficients[,4]
  pval<-p.adjust(p, method="fdr")[2:q]

  #Find first and last instances of significant dropout
  none<-0 #if no sig dropout
  nchange<-0 #number of change-points
  r<-q-1
  #prop<-matrix(nrow=s, ncol=r)
  q.d1<-rep(0,r)#comparison where first significant
  q.d2<-rep(0,r)#comparison where last significant
  for(j in 1:r){
    if(pval[j]<alpha){
      q.d1[j]<-1
      break
    }
  }
  if(sum(q.d1)==0){
     none<-1
  }
  for(k in r:1){
    if(pval[k]<alpha){
      q.d2[k]<-1
      break
    }
  }
  nchange<-sum(q.d1,q.d2) #add up start and end
  same<-sum(abs(q.d1-q.d2)) #is there only one point of change?
  if(sum(q.d1)==1 & same==0) nchange<-1
  if(same!=0 & q.d2[r]==1) nchange<-1 #trying to pick up one phase
  if(same!=0 & q.d1[1]==1) nchange<-1 #trying to pick up one phase
  if(q.d1[1]==1 & q.d2[r]==1) nchange<-0 #if constant throughout

  return(list(pval=pval, none=none, nchange=nchange, q.d1=q.d1, q.d2=q.d2))
}

#### user.specified
#user.specified<-function(dta,q){

#}
