# 22/7/2021 found an error in fit_sub part II
# 16/7 2021 replace the treatment rankings estimation in 
# function i) fit_subgroup_A
# function ii) fit_sub
# function iii) smallest_effect()
# with the following code
#if(min.eff>0){
#  ind_effect<-t_label[1] or reference treatment
#}else{ind_effect<-t_neg}



# 23/05/2021 created by Kim M Lee, King's College London
# This script contains data generating function and analysis methods.
# Please check the note on line 291 for a potential issue when using method D

library(stringr)
library(multiwayvcov)
require(sandwich)
# -------------------- #
# response probability #
# -------------------- #
res_probability<-function(phi,alpha=-1.36){
  exp(alpha+ phi)/(1+ exp(alpha+ phi) )
}

# logit for finding phi #
find_phi<-function(p, alpha=-1.36){log(p/(1-p)) - alpha } 

# ------------------------------- #
# generate data for each subgroup #
# ------------------------------- #
generate_subset_data<-function(k, size_pattern., pattern., res_probability_all.){
  # size_pattern.=size_pattern;   pattern.=pattern;  res_probability_all.=res_probability_all
  pattern_s<-pattern.[[k]]
  sp<-size_pattern.[k]
  res_p<-res_probability_all.[k, pattern_s]
  
  assigned_treatment<-t(rmultinom(sp, 1, 
                                  rep(1/length(pattern_s), length(pattern_s))))
  
  colnames(assigned_treatment)<-paste0("t",pattern_s)
  
  treatment_label<-apply(assigned_treatment, 1, 
                         function(x) x[which(x==max(x))]<-pattern_s[which(x==max(x))] )
  
  responses<-lapply(1:length(pattern_s), 
                    function(j)rbinom(sum(treatment_label==pattern_s[j]), 1, res_p[j]))
  assigned_treatment<-unlist(lapply(1:length(pattern_s),function(j)rep(pattern_s[j],length(responses[[j]])) ))
  
  assigned_pairwise<-unlist(lapply(1:length(pattern_s),function(j){rep(pattern_s[-j], sum(treatment_label==pattern_s[j])) } ) )
  
  rep_assigned_t<-rep(assigned_treatment, each=length(pattern_s)-1)
  assigned_pairwise<-paste(rep_assigned_t, assigned_pairwise, sep="_")
  
  pattern_lab<-rep(k,  sp)#unlist(sapply(1:no_pattern, function(j)rep(j,  sp)) )
  responses<-unlist(responses)
  return(list(responses=responses, 
              treatment_label=assigned_treatment, 
              assigned_pairwise=assigned_pairwise, pattern_lab=pattern_lab))
}


# -------------------------------------- #
# method_A: fit a model to each subgroup # 
# -------------------------------------- #

fit_subgroup_A<-function(alldata, no_p=no_pattern){
  fit_subA<-function(i){
    sub_data<-alldata[,i]
    sub_data$treatment<-as.factor(sub_data$treatment_label)
    my.glm<-myTryCatch(glm(responses~ treatment, data=sub_data,family="binomial"))
    if(is.null(my.glm$error) ) # if there is no error in model fit (there could be warning)
    { 
      sub_fit<-my.glm[[1]]
      t_label<-sort(unique(sub_data$treatment_label)) 
      
      max.eff<-coef(sub_fit)[1+which.max(coef(sub_fit)[-1])]
      min.eff<-coef(sub_fit)[1+which.min(coef(sub_fit)[-1])]
      t_neg<-t_label[1+which.min(coef(sub_fit)[-1] )]
      # 16/7/2021 change the following to if min is positive, then ind_effect=reference
      # if min is negative, select min
      if(min.eff>0){
        ind_effect<-t_label[1]
      }else{ind_effect<-t_neg}
      #if(max.eff>0){
      #  test_t1_better<-max.eff>abs(min.eff)
      #  if(test_t1_better==T){ind_effect<-t_label[1]}else{ind_effect<-t_neg}
      #}else{ind_effect<-t_neg}
      
    }else{ ind_effect<-my.glm$error[1]$message }
    
    
    if(is.character(ind_effect)){ 
      out<-sample(sub_data$treatment,1)
      not_fit<-1
    }else{
      out<-ind_effect
      not_fit<-0
    }
print(paste("fit_subA: best in pattern", i, "is treatment", out))    
print(my.glm)
    return(c(out,not_fit))
    
  }
  sapply(1:no_p, fit_subA)
}

# ---------------------------------------- #
# methods B fit a model for each pattern   #
# ---------------------------------------- #

# function to fit a model to data, adjust for pattern if there is more than one
fit_sub<-function(sub_data){
  #sub_data<-data_selected[[i]] is a matrix of data, #sub_data[,1]= response  #sub_data[,2]= treatment  #sub_data[,3]= pattern
  response<-sub_data[,1]; treatment_ind<-sub_data[,2]; pattern_var<-sub_data[,3]
  t_label<-sort(unique(treatment_ind)) 
  
  if(length(unique(pattern_var))==1){ # there is only one pattern
    
    my.glm<-myTryCatch( glm(response~ as.factor(treatment_ind),family="binomial"  ) )
    
    if(is.null(my.glm$error) ) # can fit the model
    { 
      sub_fit<-my.glm[[1]]
      max.eff<-coef(sub_fit)[which.max(coef(sub_fit)[2:length(t_label) ])+1]
      min.eff<-coef(sub_fit)[which.min(coef(sub_fit)[2:length(t_label) ])+1]
      t_neg<-t_label[which.min(coef(sub_fit)[2:length(t_label) ] )+1 ]
      
      # change the following to if min is positive, then ind_effect=reference
      # if min is negative, select min
      if(min.eff>0){
        ind_effect<-t_label[1]
      }else{ind_effect<-t_neg}
      
      # test which treatment has the most negative effect, if all coefficients are positive, baseline treatment is the best
      #if(max.eff>0){
      #  test_t1_better<-max.eff>abs(min.eff)
      #  if(test_t1_better==T){ind_effect<-t_label[1]}else{ind_effect<-t_neg}
      #}else{ind_effect<-t_neg} 
      
    } else{ # model can't fit
      ind_effect<-my.glm$error[1]$message }
print(paste("fit_sub for method B: best in pattern is treatment", ind_effect))
print(summary(my.glm))
print(my.glm)

    
  }else{ # if there are more than one pattern
    
    my.glm<-myTryCatch( glm(response~ as.factor(treatment_ind)+ as.factor(pattern_var), family = "binomial" ) )
    if(is.null(my.glm$error) ) 
    { 
      sub_fit<-my.glm[[1]]
      max.eff<-coef(sub_fit)[which.max(coef(sub_fit)[2:length(t_label) ])+1]
      min.eff<-coef(sub_fit)[which.min(coef(sub_fit)[2:length(t_label) ])+1]
      t_neg<-t_label[which.min(coef(sub_fit)[2:length(t_label) ] )+1 ]
      
      # 22/7/2021 update the following line
      if(min.eff>0){
        ind_effect<-t_label[1]
      }else{ind_effect<-t_neg}
      
      # test which treatment has the most negative effect, if all coefficients are positive, baseline treatment is the best
      #if(max.eff>0){
      #  test_t1_better<-max.eff>abs(min.eff)
      #  if(test_t1_better==T){ind_effect<-t_label[1]}else{ind_effect<-t_neg}
      #}else{ind_effect<-t_neg}
      
    }else{ ind_effect<-my.glm$error[1]$message}
    
    
    
  }
  
  if(is.character(ind_effect)){ 
    out<-sample(t_label,1)
    not_fit<-1
  }else{
    out<-ind_effect
    not_fit<-0
  }
  return(c(out,not_fit))
  
}


# method B include three approaches
methodB<-function(alldata=Alldata,  
                  no_p=no_pattern, 
                  size_p=size_pattern,
                  pat=pattern){
print("STARTING METHOD B")
  # generate individual data for method B
  id_v<-cumsum(size_p)
  ID<-c(1:id_v[1], unlist(sapply(1:(no_p-1),function(i){ (id_v[i]+1):id_v[i+1] } )) )
  ind_Data<-data.frame(y=unlist(alldata[1,]), treatment=unlist(alldata[2,]), pat_lab=unlist(alldata[4,]),ID=ID)
  patternI<-1:no_p
  
  # identify data that have same treatment to that of pattern k
  extract_f1<-function(j,k){which(ind_Data$treatment==pat[[k]][j]) }
  extract_f2<-function(K){unlist(sapply(1:length(pat[[K]]), function(J)extract_f1(J,K)))}
#print("methodB: ind_Data[extract_f2(i),]")
#print(ind_Data[extract_f2(i),])
  must_include<-function(i){
    exD<-ind_Data[extract_f2(i),]
    # identify if pattern k has treatment of pattern i, all true means contain all
    contain_t_f<-function(k)sapply(1:length(pat[[i]]), function(j){ length(which(pat[[k]]==pat[[i]][j]))==1 } )
    
    test_o<-sapply(patternI[-i], function(K){if( mean(contain_t_f(K))==1){K}else{NA} })
    incI<-c(i, test_o)
    incI<-incI[complete.cases(incI)]
    inc_ind<-unlist(sapply(incI, function(x) which(exD$pat_lab== x) ) )
    fit_sub(exD[inc_ind,])
  }
  
  minimal<-function(i){
    exD<-ind_Data[extract_f2(i),]
    # minimal
    test_first_e<-function(k){ if( length(which(pat[[k]]==pat[[i]][1])) ==1){k}else{NA} }
    incI<-c(i, sapply(patternI[-i],   test_first_e )  )
    incI<-incI[complete.cases(incI)]
    inc_ind<-unlist(sapply(incI, function(x) which(exD$pat_lab== x) ) )
    fit_sub(exD[inc_ind,])
print(paste("Summary of data for method B2, pattern", i))
print(dim(exD[inc_ind,]))
print(summary(exD[inc_ind,]))
print(table(exD[inc_ind,]$treatment,exD[inc_ind,]$pat_lab))
print(fit_sub(exD[inc_ind,]))
doodah
}
  
  all_direct<-function(i){
    exD<-ind_Data[extract_f2(i),]
    fit_sub(exD)
  }
  
  fit_partly_pattern_B1=sapply(1:no_p, must_include)
  fit_partly_minimal_B2=sapply(1:no_p, minimal)
  fit_partly_all_B3=sapply(1:no_p ,all_direct) 
  
  recommended_t<-rbind(fit_partly_pattern_B1[1,],
                       fit_partly_minimal_B2[1,],
                       fit_partly_all_B3[1,])
  fail_no<-rbind(fit_partly_pattern_B1[2,],
                 fit_partly_minimal_B2[2,],
                 fit_partly_all_B3[2,])
  return(list(recommended_t, fail_no))
}



# ----------------------------------------- #
# method C: fit one step model to all data  # 
# ----------------------------------------- #
fit_onestage_C<-function(alldata=Alldata, no_p=no_pattern,   q.val=qnorm(0.975), no_contrast=No_contrast){
  
  
  nma_data<-data.frame(y=unlist(alldata[1,]),
                       treatment=factor(unlist(alldata[2,])),
                       subgroup=factor(unlist(alldata[4,])))#patient_subgroup)))
  my.glm<-myTryCatch(glm(y~treatment+ subgroup,family="binomial",data=nma_data) )
  
  if(is.null(my.glm$error) ) #if do not have an error, model is fitted
  { 
    my.glm<-my.glm[[1]]
    mof<-summary(my.glm)
    std.err<-sqrt(diag(vcov(mof))[2:no_contrast]) 
    out<-cbind(Estimate=coefficients(mof)[2:no_contrast],
               model_var=std.err^2,
               z=coefficients(mof)[2:no_contrast]/std.err,
               LL=coefficients(mof)[2:no_contrast] - q.val  * std.err,
               UL=coefficients(mof)[2:no_contrast] + q.val  * std.err)
    out[which(abs(out[,1])>12),]<-NA #parameter not converged is set to NA 
  }else
  { # if there is error, do not fit model
    out<-matrix(rep(NA,(no_contrast-1)*5),nrow = no_contrast-1, ncol = 5 )
    out[1,5]<-my.glm$error[1]$message
    
  } 
  
  return(out)
}




# ---------------------------------------- #
# method D: fit duplicate data in one step # 
# ---------------------------------------- #
fit_robustSE_D<-function(alldata=Alldata, no_com=no_comparison,
                         no_p=no_pattern, 
                         no_t=no_treatment, 
                         size_p=size_pattern, 
                         no_contrast=No_contrast){
  #no_com=no_comparison; no_p=no_pattern; no_t=no_treatment; size_p=size_pattern
  dup_data<-sapply(1:no_p,function(i){cbind(rep(alldata[1,][[i]], each=no_com[i]),
                                            rep(alldata[2,][[i]], each=no_com[i]),
                                            alldata[3,][[i]] , rep(i, no_com[i]) ) } )
  dup_data<-as.data.frame(do.call(rbind, dup_data),make.names=F )
  colnames(dup_data)<-c("y", "treatment", "pairwise", "pattern")
  dup_data$treatment<-as.factor(dup_data$treatment)

  # create ID of patients
  id_v<-cumsum(size_p)
  ID<-sapply(1:(no_p-1),function(i){
    ids<-(id_v[i]+1):id_v[i+1]
    rep(ids,each=no_com[1+i])
  })
  ID<-c(rep(1:id_v[1], each=no_com[1]), unlist(ID))
  
  dup_data$id<-ID
  
  
  pairwisev<-sapply(1:no_t,function(i)paste(i, 1:no_t, sep="_"))
  v1=t(pairwisev)[upper.tri(pairwisev)]
  v2=pairwisev[upper.tri(t(pairwisev))]
  
  val<-dup_data$pairwise#val<-paste(dup_data$treatment, dup_data$pairwise, sep="_")
  for(z in 1:(no_t*(no_t-1)/2)){
    test<-which(val==v2[z])
    if(length(test)==0){}else{val[test]<-v1[z]}
  }
  dup_data$id_comparison<-factor(val)
  
  

  # note that some R version would convert the class of y to factor, if so, use the following line 292 instead of 294
  #fit_glm<-myTryCatch( glm(as.numeric(levels(y))[y]~treatment+id_comparison,family="binomial",data=dup_data) )
  
  fit_glm<-myTryCatch( glm(as.numeric(y)~treatment+id_comparison,family="binomial",data=dup_data) )
  
  if(is.null(fit_glm$error) ) #if there is no error, model is fitted 
  { 
    fit_glm<-fit_glm[[1]]
    # Calculate robust standard errors #
    cov.m1 <- cluster.vcov(fit_glm, dup_data$id)[2:no_contrast, 2:no_contrast]#vcovHC(fit_glm, type = "HC0")[2:no_contrast, 2:no_contrast]
    testcov<-myTryCatch( sqrt(diag(cov.m1)) )
    std.err <- if(is.null(testcov$warning) & is.null(testcov$error) ){testcov$value}else{rep(NaN,no_contrast-1)} 
    
    q.val <- qnorm(0.975)
    r.est <- cbind( Estimate = coef(fit_glm)[2:no_contrast], 
                    model_var=std.err^2, #"Robust SE" = std.err, 
                    z = (coef(fit_glm)[2:no_contrast]/std.err),
                    #"Pr(>|z|) "= 2 * pnorm(abs(coef(fit_glm)[2:no_contrast]/std.err), lower.tail = FALSE),
                    LL = coef(fit_glm)[2:no_contrast] - q.val  * std.err,
                    UL = coef(fit_glm)[2:no_contrast] + q.val  * std.err)
    r.est <-r.est[c(2:(no_contrast-1),1),]
    r.est[which(abs(r.est[,1])>12),]<-NA # unconverged parameter is set to NA
  } else{ # there is error in model fit
    r.est<-matrix(rep(NA,(no_contrast-1)*5),nrow = no_contrast-1, ncol = 5 )
    r.est[,5]<-fit_glm$error[1]$message
  }  
  
  return(r.est)
  
  
}

# ---------------------------------------- #
# to catch warning or error from model fit #
# ---------------------------------------- #
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# ------------------------------------- #
#   the following function is used to   # 
# select treatment from methods C and D #
# ------------------------------------- #
smallest_effect<-function(v, pat=pattern, no_p=no_pattern ){
  find.small<-function(i){ #v=output[[1]][,1];i=1
    if(all(is.na(v))){ # model is not fitted
      ind_effect<-sample(pat[[i]], 1)
      not_fit<-1
    }else{
      t_effect<-c(0,v) # changed to 0 on 16/7/2021
      t_in_pattern<-t_effect[pat[[i]]] 
      
      max.eff<-t_in_pattern[which.max(t_in_pattern)]
      min.eff<-t_in_pattern[which.min(t_in_pattern)]
      t_neg<-pat[[i]][which.min(t_in_pattern)]
      if(length(which(pat[[i]]==1))==1){ # check if t1 is in pattern i
        
        # change the following to if min is positive, then ind_effect=reference
        # if min is negative, select min
        if(min.eff>0){
          ind_effect<-1
        }else{ind_effect<-t_neg}
        
        #if(max.eff>0){
        #  test_t1_better<-max.eff>abs(min.eff)
        #  if(test_t1_better==T){ind_effect<-1}else{ind_effect<-t_neg}
        #}else{ind_effect<-t_neg}
        
      }else{ ind_effect<-t_neg}
      
      not_fit<-0
    }
    return(c(ind_effect, not_fit))
  }
  
  sapply(1:no_p, find.small)
}


method_Fail<-function(d)if(length(d)==0){Inf}else{d}

