# Ian's variant of Kim's file to run one rep
# run this with N=100 and seed=5 to demonstrate that method C mis-handles perfect prediction
#   as a treatment with 0 deaths is never selected as best
# N=1000 and seed=2 is a good example of B3 beating C

source('1 PRaCTical_functions22072021i.R')

### SETTINGS
set.seed(5)

no_treatment=10 
N=1000 # total number of patients

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
### END OF SETTINGS (MOSTLY)

pattern<-patternV
# treatment effect parameters
# scenario 1
alpha_1<-find_phi(0.2, alpha=0)
phi_1<-find_phi(seq(0.1, 0.3, length.out = no_treatment), alpha=alpha_1)
res_rate1<-res_probability(phi_1,alpha_1)


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


# Ian's additions to fit model C
# nma_data<-data.frame(y=unlist(Alldata[1,]),
#                      treatment=factor(unlist(Alldata[2,])),
#                      subgroup=factor(unlist(Alldata[4,])))#patient_subgroup)))
# table(nma_data$subgroup,nma_data$treatment)
# glm<-glm(y~factor(treatment)+factor(subgroup),data=nma_data,family=binomial(logit))
# summary(glm)



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


### REST OF PROGRAM DELETED

### IAN'S CODE TO EXPLORE WHAT WE'VE DONE ###
#  tabulate data
nma_data<-data.frame(y=unlist(Alldata[1,]),
                     treatment=factor(unlist(Alldata[2,])),
                     subgroup=factor(unlist(Alldata[4,])))#patient_subgroup)))
table(nma_data$y,nma_data$treatment,nma_data$subgroup)
# Model fit for method C:
est_method_C
# -> treatment 4 is the best (most negative)

# Best treatments for each method:
identified_best_t
