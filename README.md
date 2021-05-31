The following are the main Rscripts to run a simulation of a PRaCTical design:  

- “1 PRaCTical_functions.R” contains function that generate a data set and different analysis approaches. 

- “2 simulation_replication.R” contains a function to simulate a study with replications. 
It calls functions from the previous file. The output focuses on estimands one and two and their properties.




There are three types of scenarios:

1. Without interaction terms, all subgroups have the same baseline risk. 
The presented scenarios are S1, S1.1, S2, S3.1, S3.2, and S4.

2. Without interaction terms, subgroups have different baseline risks. The presented scenario is S1.2.

3. With interaction terms, subgroups either have the same or different baseline risks. 
The presented scenarios are S5.1, S5.2, S6.1, and S6.2.



The following files are used to run simulation with different patterns. 
Depends on the number of replications, the computation can take hours.
High performance computing facility has been used to run the following scenarios.

- “s1.R”, “s2.R”, “s3.R” and “s4.R” use both Rscripts.

- “s1_2.R” and “s_interaction.R” use only “1 PRaCTical_functions.R” because the response rates to each treatment vary across subgroups




Once having the outputs from different scenarios. The following Rscripts are for plotting the figures in the paper:

- PLOT figure 2 network.R
- PLOT figures 1 and 3.R
- PLOT figure 4.R

