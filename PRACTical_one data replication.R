source('1 PRaCTical_functions.R')
set.seed(3127) 

# 10/8/2021 added lines 81 to 92
# to print response rates within each pattern & overall response rate

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
patternV<-list(pattern1, pattern2, pattern3, pattern4, pattern5, pattern6, pattern7, pattern8)

pattern<-patternV
# treatment effect parameters
# scenario 1
alpha_1<-find_phi(0.2, alpha=0)
phi_1<-find_phi(seq(0.1, 0.3, length.out = no_treatment), alpha=alpha_1)
res_rate1<-res_probability(phi_1,alpha_1)


N=1000 # total number of patients
phi_v=phi_1 # true parameters of treatment effect
pattern=patternV # personalized randomization lists
response_prob_V=res_rate1 # response probability
prob_pattern= c(0.2, 0.2, rep(0.1, 6)) # prevelance rate of patterns
R=5 # number of trial replications

# ----- the following is copied from 2 simulation_replication.R ----- #
no_pattern<-length(pattern) # number of randomization lists

no_comparison<-sapply(1:no_pattern, function(i){length(pattern[[i]])-1}) 
# for each randomization list, the number of pairwise comparisons fixing a reference treatment. 

no_treatment<-length(phi_v) # number of treatments

res_probability_all<-matrix(rep(response_prob_V, no_pattern), ncol = no_treatment, byrow = T)
colnames(res_probability_all)<-sapply(1:no_treatment, function(i){paste0("treatment_", i)} )
rownames(res_probability_all)<-sapply(1:no_pattern, function(i){paste0("alpha_", i)} )
# response rate: row= pattern, column=treatment. All rows have same values for this scenario

# each person has prob_pattern to be allocated to one of the treatment patterns
assigned_pattern<-t(rmultinom(N, size=1, prob_pattern))
colnames(assigned_pattern)<-sapply(1:no_pattern, function(i){paste0("subgroup", i)} )

# number of patients in each subgroup that is defined by the pattern
size_pattern<-apply(assigned_pattern, 2, sum)
lambda<-prob_pattern # true prevalence rate of patterns

true.response.r<-lapply(1:no_pattern,function(i)res_probability_all[i, pattern[[i]]])
# response rates of the treatments in each pattern


true.mean.min<-lapply(1:no_pattern, function(i){
  v<-true.response.r[[i]]
  c("mean"=mean(v), "min"=min(v)) } )

true.mean.min<-do.call(cbind, true.mean.min) 
# compute the mean (and minimum value) of the treatments in each pattern 
# will be used for the performance measures about the treatment decisions

No_contrast<-length(unique(unlist(pattern))) 
# this is the total number of treatments
# to be used in method C and D as index for the nine treatment contrasts



# ** the following come from the loop of the function, exclude the list of output ** #

Alldata<-sapply(1:no_pattern, function(i){
  generate_subset_data(i, size_pattern.=size_pattern, 
                       pattern.=pattern, res_probability_all.=res_probability_all)})
# generate one dataset

# 10/8/2021 check response rate per subgroup #
sapply(1:8,function(k)table(Alldata[,k]$treatment_label, Alldata[,k]$responses)[,2]/table(Alldata[,k]$treatment_label) )
# show the response rate per arm within each pattern

table(unlist(Alldata[1,]), unlist(Alldata[2,]))
# responses to each treatment

arm.size<-apply(table(unlist(Alldata[1,]), unlist(Alldata[2,])), 2, sum)
# compute arm size

table(unlist(Alldata[1,]), unlist(Alldata[2,]))[2,]/arm.size
# overall response rate to each treatment
# 10/8/2021 end #

feq_t_subgroup<-sapply(1:no_pattern, function(i)table(Alldata[2,][[i]]))
# show how many have been randomized to a treatment arm within a pattern

feq_t<-table(unlist(Alldata[2,]))
# show how many have been randomized to each treatment arm

est_method_C<-fit_onestage_C(Alldata, no_p=no_pattern,   q.val=qnorm(0.975), no_contrast=No_contrast) # use original data
# estimates of treatment contrasts from method C

est_method_D<-fit_robustSE_D(Alldata, no_com=no_comparison, # use duplicated data
                             no_p=no_pattern, 
                             no_t=no_treatment, 
                             size_p=size_pattern,no_contrast=No_contrast)
# estimates of treatment contrasts from method D

Identify_C=smallest_effect(est_method_C[,1], pat=pattern, no_p=no_pattern)
# for method C
# row 1: the best treatment for each pattern (column)
# row 2: check if the model fails to fit, 1= fail, 0= model is fitted
# when row two=1 a random treatment is drawn for the treatment decision

Identify_D=smallest_effect(est_method_D[,1], pat=pattern, no_p=no_pattern)
# for method D
# row 1: the best treatment for each pattern (column)
# row 2: check if the model fails to fit, 1= fail, 0= model is fitted
# when row two=1 a random treatment is drawn for the treatment decision


method_A_f<-fit_subgroup_A(Alldata, no_p=no_pattern)   # fit each subgroup
# method A: fit a model to each pattern and identify the best treatment
# row 1 and 2 carry the same interpretation as methods C and D

Identify_method_B<-methodB(Alldata, no_p=no_pattern, size_p=size_pattern, pat=pattern)    
# method B with different pooling: fit a model to each pattern and identify the best treatment
# list 1 correspond to estimated best decision
# list 2 check if the model fails to fit, 1= fail, 0= model is fitted
# each row within a list corresponds to an approach of pooling

identified_best_t<-rbind(method_A=method_A_f[1,],
                         method_B1=Identify_method_B[[1]][1,],
                         method_B2=Identify_method_B[[1]][2,],
                         method_B3=Identify_method_B[[1]][3,],
                         method_C=Identify_C[1,],
                         method_D=Identify_D[1,] )
# combine estimated best treatments from all methods, row= methods, column= pattern


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
# true response rate of the estimated best treatment for each pattern(column) from each method (row)


mortality_gain<-t(sapply(1:6, function(m){identify_best_rate[m,]-true.mean.min[1,] }) )
# compute mortality reduction for each pattern(column) from each method (row)

mortality_gain # show mortality reduction for each pattern (column) from each method (row)
sum(lambda*mortality_gain[4,]) # average mortality reduction from method B3
sum(lambda*mortality_gain[5,]) # average mortality reduction from method C

# the following performance measures will be hard to see from one replication of the data
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
# combine the indicator of a model in a method fails to fit