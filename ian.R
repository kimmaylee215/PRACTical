# taken from Kim's file 2
# 5 reps of N=1000
# i haven't yet worked out how to view one simulated data set

 N=1000; phi_v=phi_1; pattern=patternV; response_prob_V=res_rate1; prob_pattern=c(0.2, 0.2, rep(0.1, 6)); R=5
  
  no_pattern<-length(pattern)
  
  no_comparison<-sapply(1:no_pattern, function(i){length(pattern[[i]])-1})
  
  no_treatment<-length(phi_v)
  
  res_probability_all<-matrix(rep(response_prob_V, no_pattern), ncol = no_treatment, byrow = T)
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
    
    Alldata<-sapply(1:no_pattern, function(i){
      generate_subset_data(i, size_pattern.=size_pattern, 
                           pattern.=pattern, res_probability_all.=res_probability_all)})

    feq_t_subgroup<-sapply(1:no_pattern, function(i)table(Alldata[2,][[i]]))
    feq_t<-table(unlist(Alldata[2,]))
    
    est_method_C<-fit_onestage_C(Alldata, no_p=no_pattern,   q.val=qnorm(0.975), no_contrast=No_contrast) # use original data
    est_method_D<-fit_robustSE_D(Alldata, no_com=no_comparison, # use duplicated data
                                 no_p=no_pattern, 
                                 no_t=no_treatment, 
                                 size_p=size_pattern,no_contrast=No_contrast)
    
    Identify_C=smallest_effect(est_method_C[,1], pat=pattern, no_p=no_pattern)
    Identify_D=smallest_effect(est_method_D[,1], pat=pattern, no_p=no_pattern)
    
    
    method_A_f<-fit_subgroup_A(Alldata, no_p=no_pattern)   # fit each subgroup
    
    
    Identify_method_B<-methodB(Alldata, no_p=no_pattern, size_p=size_pattern, pat=pattern)    
    
    identified_best_t<-rbind(method_A=method_A_f[1,],
                             method_B1=Identify_method_B[[1]][1,],
                             method_B2=Identify_method_B[[1]][2,],
                             method_B3=Identify_method_B[[1]][3,],
                             method_C=Identify_C[1,],
                             method_D=Identify_D[1,] )
    
    
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
  #rownames(method_c_property)<- rownames(method_d_property)<- sapply(2:10, function(i){paste0("phi_", i, "-phi_1")} )
  
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
            method_D_property=method_d_property,
            #ex_performance_out=ex_performance_out,
            #suggested_treatment_each=suggested_treatment_each,
            #estimator_all=estimator_all,
            #model_var_all=model_var_all,
            #all_diff_min=all_diff_min,
            method_fail_no=methodA_fail_no,
            estimand2=estimand2,
            estimand2_MCSE=estimand2_MCSE, 
            #ex_arm_size=ex_arm_size,
            overall_size=feq_treatment, #sub_size=feq_treatment_subgroup,
            Pattern=pattern)
