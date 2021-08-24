source("1 PRaCTical_functions22072021_firth.R")

### START OF SETTINGS

no_treatment=10 

# treatment patterns
pattern1<-c(2,3,5,8,10)
pattern2<-1:7
pattern3<-c(1,2,4,9,10)
pattern4<-c(1,2,3,5,6,8,10)
pattern5<-c(1,2,3,4,6,7)
pattern6<-2:10
pattern7<-1:10
pattern8<-3:10
# pattern8<-6:10 # Change 2: only poor treatments in pattern 8

# treatment effect parameters
# scenario 1.2
alpha_1.2<-find_phi(seq(0.05, 0.6, length.out = 8), alpha=0)
alpha_1.2[3]<-find_phi(0.2,alpha=0)
#alpha_1.2 <- -2 - alpha_1.2 # Change 3: reverse pattern effects
phi_1.2<-find_phi(seq(0.1, 0.3, length.out = no_treatment), alpha=alpha_1.2[3])
#phi_1.2 <- -phi_1.2 # Change 4: reverse treatment effects
phi_1.2 <- 0.001*phi_1.2 # Change 5: tiny treatment effects
#phi_1.2 <- -0.001*phi_1.2 # Change 5: tiny & reversed treatment effects

# Change 6: exactly reverse the order of the treatments
# shouldn't affect results at all
# pattern1 <- rev(11-pattern1)
# pattern2 <- rev(11-pattern2)
# pattern3 <- rev(11-pattern3)
# pattern4 <- rev(11-pattern4)
# pattern5 <- rev(11-pattern5)
# pattern6 <- rev(11-pattern6)
# pattern7 <- rev(11-pattern7)
# pattern8 <- rev(11-pattern8)
# phi_1.2 <- rev(phi_1.2)

Nobs <- 1000
# Nobs <- 10000 # Change 1: large sample to avoid perfect prediction

Nreps <- 1000 # minimum 2, fails for 1

set.seed(1)

#res_rate1.2<-res_probability(phi_1.2,alpha_1.2[3])
res_rate_mat<-t(sapply(alpha_1.2, function(i)res_probability(phi_1.2,i)))
res_rate_mat<-res_rate_mat[c(3,1,2,4:8),] # I don't know why Kim did this but it doesn't matter

### END OF SETTINGS

patternV<-list(pattern1, pattern2, pattern3, pattern4, pattern5, pattern6, pattern7, pattern8)
pattern<-patternV

simulation<-function(N, phi_v, pattern, res_probability_all,
                     prob_pattern,  R=10){
  # N=1000; phi_v=phi_1; pattern=pattern3; response_prob_V=res_rate1; prob_pattern=rep(0.1, length(pattern3) ); R=5
  
  no_pattern<-length(pattern)
  
  no_comparison<-sapply(1:no_pattern, function(i){length(pattern[[i]])-1})
  
  no_treatment<-length(phi_v)
  
  #res_probability_all<-matrix(rep(response_prob_V, no_pattern), ncol = 10, byrow = T)
  colnames(res_probability_all)<-sapply(1:no_treatment, function(i){paste0("treatment_", i)} )
  rownames(res_probability_all)<-sapply(1:no_pattern, function(i){paste0("alpha_", i)} )
  
  
  # each person has prob_pattern to be allocated to one of the treatment patterns
  assigned_pattern<-t(rmultinom(N, size=1, prob_pattern))
  colnames(assigned_pattern)<-sapply(1:no_pattern, function(i){paste0("subgroup", i)} )
  
  # number of patients in each subgroup that is defined by the pattern
  size_pattern<-apply(assigned_pattern, 2, sum)
  lambda<-prob_pattern#size_pattern/N
  
  true.response.r<-lapply(1:no_pattern,function(i)res_probability_all[i, pattern[[i]]])
  true.mean.min<-lapply(1:no_pattern, function(i){
    v<-true.response.r[[i]]
    c("mean"=mean(v), "min"=min(v)) } )
  
  true.mean.min<-do.call(cbind, true.mean.min)
  No_contrast<-length(unique(unlist(pattern))) 
  
  gen.data<-function(j){    
# print(paste("Starting gen.data, replicate",j))
    Alldata<-sapply(1:no_pattern, function(i){
      generate_subset_data(i, size_pattern.=size_pattern, 
                           pattern.=pattern, res_probability_all.=res_probability_all)})
    feq_t_subgroup<-sapply(1:no_pattern, function(i)table(Alldata[2,][[i]]))
    feq_t<-table(unlist(Alldata[2,]))
# print("Data for pattern 1: outcome by treatment")
# print(table(Alldata[,1]$responses,Alldata[,1]$treatment_label))
# print("Data for pattern 2: outcome by treatment")
# print(table(Alldata[,2]$responses,Alldata[,2]$treatment_label))
# print("Data for pattern 3: outcome by treatment")
# print(table(Alldata[,3]$responses,Alldata[,3]$treatment_label))
# print("Data for pattern 4: outcome by treatment")
# print(table(Alldata[,4]$responses,Alldata[,4]$treatment_label))
# print("Data for pattern 5: outcome by treatment")
# print(table(Alldata[,5]$responses,Alldata[,5]$treatment_label))
# print("Data for pattern 6: outcome by treatment")
# print(table(Alldata[,6]$responses,Alldata[,6]$treatment_label))
# print("Data for pattern 7: outcome by treatment")
# print(table(Alldata[,7]$responses,Alldata[,7]$treatment_label))
# print("Data for pattern 8: outcome by treatment")
# print(table(Alldata[,8]$responses,Alldata[,8]$treatment_label))
    est_method_C<-fit_onestage_C(Alldata, no_p=no_pattern,   q.val=qnorm(0.975), no_contrast=No_contrast) # use original data
    est_method_D<-fit_robustSE_D(Alldata, no_com=no_comparison, # use duplicated data
                                 no_p=no_pattern, 
                                 no_t=no_treatment, 
                                 size_p=size_pattern,no_contrast=No_contrast)
    
    Identify_C=smallest_effect(est_method_C[,1], pat=pattern, no_p=no_pattern)
    Identify_D=smallest_effect(est_method_D[,1], pat=pattern, no_p=no_pattern)
# print("est_method_C")
# print(est_method_C)
# print("Identify_C")
# print(Identify_C)    
# print("est_method_D")
# print(est_method_D)
# print("Identify_D")
# print(Identify_D)    

    method_A_f<-fit_subgroup_A(Alldata, no_p=no_pattern)   # fit each subgroup
# print("method_A_f")
# print(method_A_f)    
    
    Identify_method_B<-methodB(Alldata, no_p=no_pattern, size_p=size_pattern, pat=pattern)    
# print("Identify_method_B")
# print(Identify_method_B)    
    
    identified_best_t<-rbind(method_A=method_A_f[1,],
                             method_B1=Identify_method_B[[1]][1,],
                             method_B2=Identify_method_B[[1]][2,],
                             method_B3=Identify_method_B[[1]][3,],
                             method_C=Identify_C[1,],
                             method_D=Identify_D[1,] )
    
# print("identified_best_t")
# print(identified_best_t)

    identify_bestR<-function(k){
      #v<-sapply(1:6, function(m){ which( pattern [[k]]==identified_best_t[m,k]) } ) 
      v<-sapply(1:6, function(m){ 
        o1<-which( pattern [[k]]==identified_best_t[m,k])
        if(length(o1)==0){NA}else{o1}
      } ) 
      t.rate<-true.response.r[[k]]
      t.rate<-t.rate[v]
      names(t.rate)<-rownames(identified_best_t)
      return(t.rate)
    }
    
    identify_best_rate<-sapply(1:no_pattern,identify_bestR)
    
    mortality_gain<-t(sapply(1:6, function(m){identify_best_rate[m,]-true.mean.min[1,] }) )
    better_treatment_I<-mortality_gain<0
    
    
    
    diff_min<-t(sapply(1:6, function(m){ identify_best_rate[m,]-true.mean.min[2,]  }) )
    
    best_treatment_I<-diff_min==0
    
    nearbest_treatment_5<-diff_min-0.05 <= 0
    nearbest_treatment_10<-diff_min-0.1 <= 0
    
    rownames(mortality_gain)<-rownames(better_treatment_I)<-rownames(identified_best_t)
    rownames(diff_min)<-rownames(best_treatment_I)<-rownames(identified_best_t)
    
    rownames(nearbest_treatment_5)<-rownames(nearbest_treatment_10)<-rownames(identified_best_t)
    
    estimand2<-list(mortality_gain=mortality_gain,
                    better_treatment_I=better_treatment_I,
                    best_treatment_I=best_treatment_I,
                    nearbest_treatment_5=nearbest_treatment_5,
                    nearbest_treatment_10=nearbest_treatment_10,
                    diff_min=diff_min )
    #names(measure)<-c("pattern1", "pattern2", "pattern3", "pattern4", "pattern5", "pattern6", "pattern7", "pattern8")
    identify_fail<-rbind(method_A_f[2,], Identify_method_B[[2]], Identify_C[2,], Identify_D[2,])
    
    list(identified_best_t=identified_best_t,
         est_method_C=est_method_C,
         est_method_D=est_method_D,
         performance_m=estimand2,
         identify_fail=identify_fail, 
         #methodA_notfit=method_A_f[2,], 
         feq_t_subgroup=feq_t_subgroup, feq_t=feq_t)
    
  }
  output_replication<-lapply(1:R, function(k){
    print(k)
    gen.data(k)})
  
  methodA_fail_no<-rbind(Method_A=sapply(1:no_pattern, function(k){ mean(sapply(1:R, function(z){output_replication[[z]]$identify_fail[1,k]}))} ),
                         Method_B1=sapply(1:no_pattern, function(k){ mean(sapply(1:R, function(z){output_replication[[z]]$identify_fail[2,k]}))} ),
                         Method_B2=sapply(1:no_pattern, function(k){ mean(sapply(1:R, function(z){output_replication[[z]]$identify_fail[3,k]}))} ),
                         Method_B3=sapply(1:no_pattern, function(k){ mean(sapply(1:R, function(z){output_replication[[z]]$identify_fail[4,k]}))} ),
                         Method_C=sapply(1:no_pattern, function(k){ mean(sapply(1:R, function(z){output_replication[[z]]$identify_fail[5,k]}))} ),
                         Method_D=sapply(1:no_pattern, function(k){ mean(sapply(1:R, function(z){output_replication[[z]]$identify_fail[6,k]}))} ) )
  
  estimator_method_C<-do.call(rbind, lapply(1:R, function(z){output_replication[[z]]$est_method_C[,1]}))
  estimator_method_D<-do.call(rbind, lapply(1:R, function(z){output_replication[[z]]$est_method_D[,1]}))
  estimator_all<-list(estimator_method_C=estimator_method_C, estimator_method_D=estimator_method_D)
  
  model_var_method_C<-do.call(rbind, lapply(1:R, function(z){output_replication[[z]]$est_method_C[,2]}))
  model_var_method_D<-do.call(rbind, lapply(1:R, function(z){output_replication[[z]]$est_method_D[,2]}))
  model_var_all<-list(model_var_method_C=model_var_method_C, model_var_method_D=model_var_method_D)
  
  # compute the property of estimator
  com_property<-function(out_one, q){
    if(all( is.na(out_one[,1]) ) ){rep(NA,6)}else{
      val<- out_one[,1]
      val<- as.numeric(val[complete.cases(val)])
      t.diff<-(phi_v[q+1]-phi_v[1])
      bias<-mean(val-t.diff); var.s<-var(val)
      
      meanv2<-mean(as.numeric(out_one[,2]), na.rm = T)
      v2<-out_one[,-c(1,3)]
      #meanv2<-apply(v2, 2, function(x){
      #  if(all(is.numeric(x))){mean(x)}else{
      #    indx<-which(is.na(str_extract(x, "[0-9]+")))
      #    mean(as.numeric(x[-indx])) } }  )
      pw<-v2[which(is.na(str_extract(v2[,3], "[0-9]+"))), 3]

      coverage_ind<-rbind(out_one[,5] >=t.diff, out_one[,4] <=t.diff)
      coverage_count<-apply(coverage_ind, 2, sum)
      coverage_prob<-length(which(coverage_count==2)) / length(which(is.na(coverage_count)==F))
      
      MSE<-mean((val-t.diff)^2)
      MCSE_mse<-sqrt( var((val-t.diff)^2)/R )
      list(pw, c(bias=bias, empirical_var=var.s, coverage_prob=coverage_prob, mse=MSE, 
                 MCSE_bias=sqrt(var.s/R), 
                 MCSE_cov_p=sqrt(coverage_prob*(1-coverage_prob)/R), 
                 MCSE_MSE=MCSE_mse,
                 ex_model_var=meanv2,  
                 fail.no=length(which(is.na(out_one[,1])))  ) )
    }
  }
  
  estimator_prop<-function(q){ #q=1; out_one=method_C
    method_C<-do.call(rbind, lapply(1:R, function(z){output_replication[[z]]$est_method_C[q,]}))
    method_D<-do.call(rbind, lapply(1:R, function(z){output_replication[[z]]$est_method_D[q,]}))
    method_Co=com_property(method_C, q)
    method_Do=com_property(method_D, q)
    list( method_C_warning=method_Co[[1]],  method_D_warning=method_Do[[1]] ,
          property=cbind(method_Cp=method_Co[[2]],  method_Dp=method_Do[[2]]  ) )
    
    
  }
  
  estimator_property<-lapply(1:(No_contrast-1),estimator_prop)
  names(estimator_property)<-sapply(2:No_contrast, function(i)paste0("phi",i))
  method_C_property<-sapply(1:(No_contrast-1), function(i)estimator_property[[i]][[3]][,1])
  method_D_property<-sapply(1:(No_contrast-1), function(i)estimator_property[[i]][[3]][,2])
  
  method_c_property<-t(method_C_property)#[,c(1,2,5,6)]
  method_d_property<-t(method_D_property)#[,c(1,2,5,6)]
  
  # identify the suggested treatment
  suggested_treatment<-function(q){
    all_out<-do.call(cbind,lapply(1:R, function(z){output_replication[[z]]$identified_best_t[,q]}))
    t(apply(all_out, 1,function(x){
      if(all(is.na(x))){rep(NA,3)}else{
        quantile(x, probs =c(0.25,0.5,0.75), type = 1, na.rm = T ) } } ) )
  }
  suggested_treatment_each<-lapply(1:no_pattern, suggested_treatment)
  names(suggested_treatment_each)<-sapply(1:no_pattern, function(i)paste0("pattern",i))

  # performance of each method
  # performance of each method
  ex_performance<-function(q,k){
    mat_all<-do.call(cbind,lapply(1:R, function(z){
      output_replication[[z]]$performance_m[[q]][,k] } ) )
    apply(mat_all, 1, function(x)mean(x, na.rm = T))
  }
  
  ex_performance_out<-lapply(names(output_replication[[1]]$performance_m)[-6],
                             function(j)sapply(1:no_pattern, function(i)ex_performance(j,i) ) )
  names(ex_performance_out)<-names(output_replication[[1]]$performance_m)[-6]
  
  estimand2<-do.call(cbind,lapply(ex_performance_out, function(x){apply(x,1, function(y){sum(y*lambda)})}) )
  
  estimand2_MCSE<-sqrt(estimand2[,-1]*(1-estimand2[,-1])/R )
  
  all_diff_min<-lapply(1:R, function(z){output_replication[[z]]$performance_m$diff_min })  
  
  mortality_gain<-do.call(rbind, lapply(1:R, function(z){
    apply(output_replication[[z]]$performance_m$mortality_gain, 1, function(x){sum(x*lambda)}) } )   )
  
  estimand2_MCSE<-cbind(Mortaliy=apply(mortality_gain, 2,function(x){ 
    sqrt(var(x[complete.cases(x)])/length(x[complete.cases(x)])) }), estimand2_MCSE)
  
  size_per_arm<-size_pattern/(no_comparison+1)
  
  each_t<-function(k){
    
    com_size_each<-function(i){
      
      v<-pattern[[i]]
      m<-cbind(v,rep(size_per_arm[i],length(v)) )
      
      #m<-size_per_gp[[i]]
      
      tv<-which(m[,1]==k)
      if(length(tv)==0){ 0 }else{ m[tv,2]}
      
    }
    sum(sapply(1:length(pattern), com_size_each))
  }
  
  ex_arm_size<-sapply(1:no_treatment, each_t)
  
  
  t_feq<-t(sapply(1:R, function(r)output_replication[[r]]$feq_t))
  feq_treatment<-apply(t_feq, 2, summary)[c(1,3,4,6), ]
  
  comp_no_t_mean<-function(i){
    #all_feq<-lapply(1:R, function(r)as.data.frame(output[[r]]$feq_t_subgroup[[i]]))
    
    count_val<-function(r){
      v<-as.data.frame( output_replication[[r]]$feq_t_subgroup[[i]] )
      count_v<-sapply(1:length(pattern[[i]]), function(j){
        test_arg<-which(v==pattern[[i]][j] )
        if(length(test_arg)==1){v[test_arg,2 ]}else{0} } )
      return(count_v)}
    mean_count<-apply(sapply(1:R, count_val), 1, mean)
    names(mean_count)<-pattern[[i]]
    return(mean_count)
  }
  
  
  feq_treatment_subgroup<-sapply(1:no_pattern, comp_no_t_mean)
  
  out<-list(method_C_property=method_c_property,
            method_D_property=method_d_property,#estimator_property=estimator_property, 
            ex_performance_out=ex_performance_out,
            suggested_treatment_each=suggested_treatment_each,
            estimator_all=estimator_all,
            model_var_all=model_var_all,
            all_diff_min=all_diff_min,
            method_fail_no=methodA_fail_no,
            estimand2=estimand2,
            estimand2_MCSE=estimand2_MCSE, 
            ex_arm_size=ex_arm_size,
            overall_size=feq_treatment, sub_size=feq_treatment_subgroup,
            Pattern=pattern)
  return(out)
line <- readline(prompt="Press [enter] to continue")
}


# 100 reps take ~50s
lambda = c(0.2, 0.2, rep(0.1, 6))
scenario_out<-simulation(N=Nobs, phi_v=phi_1.2, 
                         pattern=patternV, 
                         res_probability_all=res_rate_mat, 
                         prob_pattern= lambda, R=Nreps)
filename<-paste0("scenario_A_3.RData")
save(scenario_out, file=filename) # task_id=1

########################################################################
# Ian's addition to view output
Nobs
Nreps

# patterns
patternV

# pattern prevalences
lambda

# treatment effects
phi_1.2

# p(y|row=pattern, col=treatment) 
res_rate_mat

# overall PMs
scenario_out$estimand2 
scenario_out$estimand2_MCSE

# Better treatment by pattern
scenario_out$ex_performance_out$better_treatment_I
round(sqrt(scenario_out$ex_performance_out$better_treatment_I*(1-scenario_out$ex_performance_out$better_treatment_I)/Nreps),4)

# Best treatment by pattern
scenario_out$ex_performance_out$best_treatment_I
round(sqrt(scenario_out$ex_performance_out$best_treatment_I*(1-scenario_out$ex_performance_out$best_treatment_I)/Nreps),4)
########################################################################
