########################################################################################################
# Title:          CEA model eptinezumab v. erenumab
# Programmer:     Mark Bounthavong
# Version:        1.11
# Date:           28 December 2023
# Updated:        09 January 2024
# Updated by:     Mark Bounthavong
########################################################################################################

#################################################################
##### Description
# Building decision tree model for Eptinezumab CEA using R
# 
# Note (12-29-2023): 
# - Started to add indirect costs (~line 96)
# - PSA added
# 
# NOTE (12-30-2023):
# Connection issue with "rt" error: Changed life table folder to 
# read director from "urfile" object earlier in the code sequence.
# 
# NOTE (12-31-2023):
# Changed erenumb cost to NC price
# 
# NOTE (01-08-2024):
# Change the cost of office-based, ER, and inpatient visits 
# to reflect monthly costs
# Strategy names are reversed? It's because the order of 
# preferred strategy goes with the lowest cost
# item. Hence, the authors of the package used the 
# SOC or reference group as the first one.
#################################################################

##### Clear memory
rm(list=ls())


## Windows:
# Set working directory - Windows
setwd("C:\\Users\\mbounthavong\\Dropbox\\Projects\\CEA Eptinezumab migraines\\R codes")

# Set working directory - MAC
setwd("/Users/mbounthavong/Library/CloudStorage/Dropbox/Projects/CEA Eptinezumab migraines/R codes")


##### Load functions needed for plotting a CE frontier
source("https://raw.githubusercontent.com/mbounthavong/Cost-Effectiveness-Analysis-of-Eptinezumab-and-Erenumab-for-Migraines/main/R%20Functions/CEA_functions2.R")


#### Load libraries
library("darthtools")
library("pacman")
p_load("dplyr", "tidyr", "reshape2", "devtools", "scales", "ellipse", "ggplot2", "ggrepel", "gridExtra", "lazyeval", "igraph", "truncnorm", "ggraph", "reshape2", "patchwork", "knitr", "stringr", "diagram", "dampack", "tidyverse", "psych")   



##################################################################################################################
################################################ Decision tree ###################################################
##################################################################################################################

################################
## EM model
################################

strategies <- c("EPTI", "EREN")


### Probabilities
p.reduction_epti = 0.66                         # Probability of a reduction of 50% or greater = EPTI
p.no_reduction_epti = 1 - p.reduction_epti      # Probability of a reduction of < 50% - EPTI

p.reduction_eren = 0.27                         # Probability of a reduction of 50% or greater = EREN
p.no_reduction_eren = 1 - p.reduction_eren      # Probability of a reduction of < 50% - EREN


### Drug costs
c.eptinezumab = 198.3867  # Cost of eptinezumab (1 month)
c.admin = 121.65          # Cost of administering epti

c.erenumab = 136.90       # Cost of erenumab (1 month)


### Utility scores
u.response_epti = 0.7938 + (-0.0189 * (8.7 - 3.9))
u.no_response_epti = 0.7938 + (-0.0189 * 8.7)

u.response_eren = 0.7938 + (-0.0189 * (8.3 - 3.2))
u.no_response_eren = 0.7938 + (-0.0189 * 8.3)

### Decision tree duration
tree_length = 6         # Length of the decision tree (duration in months)



### Total Costs
tc.epti_1 = tree_length*(c.eptinezumab + c.admin)      # Total costs for EPTI - branch 1
tc.epti_2 = tree_length*(c.eptinezumab + c.admin)      # Total costs for EPTI - branch 2

tc.eren_1 = tree_length*(c.erenumab)                   # Total costs for EREN - branch 1
tc.eren_2 = tree_length*(c.erenumab)                   # Total costs for EREN - branch 2


### Total QALYs
te.epti_1 = (u.response_epti / 2)      # EPTI - Total QALYs; divide by 2 for 6 months
te.epti_2 = (u.no_response_epti / 2)   # EPTI - Total QALYs; divide by 2 for 6 months

te.eren_1 = (u.response_eren / 2)      # EREN - Total QALYs; divide by 2 for 6 months
te.eren_2 = (u.no_response_eren / 2)   # EREN - Total QALYs; divide by 2 for 6 months

### WTP
wtp = 50000                            # Willingness-to-pay


### Indirect costs (Societal perspective)
pc.indirect <- 136         # Work productivity cost (assume 4-hours lost per migraine attach)
pc.epti     <- (8.7 - 3.9) * pc.indirect
pc.eren     <- (8.3 - 3.2) * pc.indirect





##################################################
#### OPEN TREE Method -- Decision Tree ####
##################################################

### Expected costs
C.EPTI  <- c(prod(c(p.reduction_epti)),
             prod(c(p.no_reduction_epti))) %*%
            c(tc.epti_1, 
              tc.epti_2)

C.EREN  <- c(prod(c(p.reduction_eren)),
             prod(c(p.no_reduction_eren))) %*%
           c(tc.eren_1, 
             tc.eren_2)

### Expected QALYs
E.EPTI  <- c(prod(c(p.reduction_epti)),
             prod(c(p.no_reduction_epti))) %*%
           c(te.epti_1, 
             te.epti_2)

E.EREN  <- c(prod(c(p.reduction_eren)),
             prod(c(p.no_reduction_eren))) %*%
           c(te.eren_1, 
             te.eren_2)

### Total costs and QALYs
TC <- c(C.EPTI, C.EREN)
TE <- c(E.EPTI, E.EREN)

names(TC) <- names(TE) <- c("EPTI", "EREN")
TC
TE


### Incremental costs
DC <- TC - TC[2]
DC

### Incremental effects
DE <- TE - TE[2]
DE


### Incremental Cost-Effectiveness Ratio
ICER <- DC / DE
ICER


### Net Monetary Benefit
NMB <- TE * wtp - TC
NMB

### Create Full ICER table
table.icer <- cbind(TC, TE, DC, DE, round(ICER, 2))
table.icer <- as.data.frame(table.icer)
colnames(table.icer) <- c("Costs", "QALYs", "Inc. Costs", "Inc. QALYs", "ICER")
rownames(table.icer) <- strategies
table.icer



##################################################################################################################
################################################ Markov model ####################################################
##################################################################################################################

## General setup 
cycle_length <- 1/12   # cycle length equal to one year (use 1/12 for monthly)
n_age_init   <- 40     # age at baseline
n_age_max    <- 50    # maximum age of follow up
n_cycles     <- (n_age_max - n_age_init)/cycle_length   # time horizon, number of cycles


# Age labels 
v_age_names  <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                      1:(1/cycle_length), 
                      sep = ".")

# the 3 health states of the model:
v_names_cycles  <- paste("cycle", 0:n_cycles)    # cycle names

v_names_states <- c("EM0",  # EM off preventative treatment (EM0)
                    "EM1",  # EM on preventative treatment (EM1)
                    "D")    # Death (D)

n_states <- length(v_names_states)   # number of health states 


### Discounting factors 
d_c <- d_e <- 0.03  # annual discount rate for costs and QALY


### Strategies 
v_names_str <- c("Eptinezumab",      # store the strategy names
                 "Erenumab") 
n_str       <- length(v_names_str)   # number of strategies


## Within-cycle correction (WCC) using Simpson's 1/3 rule 
v_wcc  <- gen_wcc(n_cycles = n_cycles, method = "Simpson1/3")


### Transition probabilities
tp.EM0_EM1      <- 0
tp.EM0_D        <- 0.0002826
tp.EM0_EM0      <- 1 - (tp.EM0_EM1 + tp.EM0_D)


tp.EM1_EM0_EREN <- 0.0055  # Discontinuation due to AE
tp.EM1_D        <- 0.0002826
tp.EM1_EM1      <- 1 - (tp.EM1_EM0_EREN + tp.EM1_D)

tp.D            <- 1 

tp.EM1_EM0_EPTI <- 0.0091  # Discontinuation due to AE


### US Life Table 2020
urlfile = ("https://raw.githubusercontent.com/mbounthavong/Cost-Effectiveness-Analysis-of-Eptinezumab-and-Erenumab-for-Migraines/main/Data/2020_life_table.csv")
life_table <- read.csv(urlfile, header = TRUE)

### Age-dependent mortality rates
v_r_mort_by_age <- life_table %>% 
  dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
  dplyr::select(Total) %>%
  as.matrix()


### Markov model parameters
p.female = 0.843             # Proportion of population that are female
p.male = 1 - p.female        # Proportion of population that are male

c.symptom_relief = 9.10      # Cost of symptom relief

p.office_visit_m = 0.031     # Probability of migraine-associated office-visit (male)
p.office_visit_f = 0.064     # Probability of migraine-associated office-visit (female)
c.office_visit = 253.33      # Cost of migraine-associated office-visit

p.hospitalization_m = 0.003  # Probability of migraine-associated hospitalization (male)
p.hospitalization_f = 0.002  # Probability of migraine-associated hospitalization (female)
c.hospitalization = 7364     # Cost of migraine-associated hospitalization

p.ed_visit_m = 0.007         # Probability of migraine-associated ED visit (male)
p.ed_visit_f = 0.013         # Probability of migraine-associated ED visit (female)
c.ed_visit = 1114.55         # Cost if migraine-associated ED visit

### State costs
C_EM0 <- p.female*((p.office_visit_f*c.office_visit) + (p.hospitalization_f*c.hospitalization) + (p.ed_visit_f*c.ed_visit)) + p.male*((p.office_visit_m*c.office_visit) + (p.hospitalization_m*c.hospitalization) + (p.ed_visit_m*c.ed_visit)) + c.symptom_relief

C_EM1 <- p.female*((p.office_visit_f*c.office_visit) + (p.hospitalization_f*c.hospitalization) + (p.ed_visit_f*c.ed_visit)) + p.male*((p.office_visit_m*c.office_visit) + (p.hospitalization_m*c.hospitalization) + (p.ed_visit_m*c.ed_visit)) 

C_D <- 0


### State utilities - EPTI
U_EM0_epti <-  u.no_response_epti / 12    # Adjusted from 12 months to 1 month

U_EM1_epti <-  u.response_epti / 12       # Adjusted from 12 months to 1 month 

U_D <- 0


### State utilities - EREN
U_EM0_eren <-  u.no_response_eren / 12    # Adjusted from 12 months to 1 month

U_EM1_eren <-  u.response_eren / 12       # Adjusted from 12 months to 1 month

U_D <- 0


### Discount weight for costs and effects 
v_dwc     <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
v_dwe     <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))


# Process model inputs 
## Age-specific transition rates to the Dead state for all cycles 
v_r_Dage  <- rep(v_r_mort_by_age, each = 1/cycle_length)
# Name age-specific mortality vector 
names(v_r_Dage) <- v_age_names

### Rate to probability
v_p_Dage   <- rate_to_prob(v_r_Dage,  t = cycle_length) # Age-specific mortality risk All states to Death 

### Construct Markov model figure
m_P_diag <- matrix(0, nrow = n_states, ncol = n_states, dimnames = list(v_names_states, v_names_states))
m_P_diag["EM0", "EM0" ] = ""
m_P_diag["EM0", "EM1" ] = "" 
m_P_diag["EM0", "D" ] = ""
m_P_diag["EM1", "EM1" ] = ""
m_P_diag["EM1", "EM0" ] = ""
m_P_diag["EM1", "D" ] = ""
m_P_diag["D", "D" ] = ""
layout.fig <- c(2, 1)
plotmat(t(m_P_diag), t(layout.fig), self.cex = 0.5, curve = 0, arr.pos = 0.8,  
        latex = T, arr.type = "curved", relsize = 0.85, box.prop = 0.8, 
        cex = 0.8, box.cex = 0.7, lwd = 1)


# All starting healthy
v_m_init_EPTI <- c("EM0" = p.no_reduction_epti, "EM1" = p.reduction_epti, "D" = 0)   # Change based on decision tree
v_m_init_EPTI

v_m_init_EREN <- c("EM0" = p.no_reduction_eren, "EM1" = p.reduction_eren, "D" = 0)   # Change based on decision tree
v_m_init_EREN


### Initialize cohort trace for EREN and EPTI (These will be different based on the decision tree)
m_M_EPTI <- matrix(0, 
                   nrow = (n_cycles + 1), ncol = n_states, 
                   dimnames = list(v_names_cycles, v_names_states))

m_M_EREN <- matrix(0, 
                  nrow = (n_cycles + 1), ncol = n_states, 
                  dimnames = list(v_names_cycles, v_names_states))


# Store the initial state vector in the first row of the cohort trace
m_M_EPTI[1, ] <- v_m_init_EPTI
m_M_EREN[1, ] <- v_m_init_EREN


## Create transition probability arrays for strategy EREN  
### Initialize transition probability array for strategy EREN 
# All transitions to a non-death state are assumed to be conditional on survival
a_P_EPTI <- array(0,  # Create 3-D array
                  dim = c(n_states, n_states, n_cycles),
                  dimnames = list(v_names_states, v_names_states, 
                                  v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 

a_P_EREN <- array(0,  # Create 3-D array
                 dim = c(n_states, n_states, n_cycles),
                 dimnames = list(v_names_states, v_names_states, 
                                 v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 


### Fill in array
## EPTI
# From EM0
a_P_EPTI["EM0", "EM0", ] <- (1 - v_p_Dage) * (1 - tp.EM0_EM1)
a_P_EPTI["EM0", "EM1", ] <- (1 - v_p_Dage) *      tp.EM0_EM1
a_P_EPTI["EM0", "D", ]   <-      v_p_Dage

# From EM1
a_P_EPTI["EM1", "EM1", ] <- (1 - v_p_Dage) * (1 - tp.EM1_EM0_EPTI)
a_P_EPTI["EM1", "EM0", ] <- (1 - v_p_Dage) *      tp.EM1_EM0_EPTI
a_P_EPTI["EM1", "D", ]   <-      v_p_Dage

# From Dead
a_P_EPTI["D", "D", ] <- 1


## EREN
# From EM0
a_P_EREN["EM0", "EM0", ] <- (1 - v_p_Dage) * (1 - tp.EM0_EM1)
a_P_EREN["EM0", "EM1", ] <- (1 - v_p_Dage) *      tp.EM0_EM1
a_P_EREN["EM0", "D", ]   <-      v_p_Dage

# From EM1
a_P_EREN["EM1", "EM1", ] <- (1 - v_p_Dage) * (1 - tp.EM1_EM0_EREN)
a_P_EREN["EM1", "EM0", ] <- (1 - v_p_Dage) *      tp.EM1_EM0_EREN
a_P_EREN["EM1", "D", ]   <-      v_p_Dage

# From Dead
a_P_EREN["D", "D", ] <- 1


### Check if transition array and probabilities are valid
# Check that transition probabilities are in [0, 1]
check_transition_probability(a_P_EPTI, verbose = TRUE)
check_transition_probability(a_P_EREN,  verbose = TRUE)

# Check that all rows sum to 1
check_sum_of_transition_array(a_P_EPTI, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(a_P_EREN,  n_states = n_states, n_cycles = n_cycles, verbose = TRUE)


###########################################
### Run the Markov model
###########################################
# Iterative solution of age-dependent cSTM
for(t in 1:n_cycles){
  ## Fill in cohort trace
  # For EPTI
  m_M_EPTI[t + 1, ]  <- m_M_EPTI[t, ]  %*% a_P_EPTI[, , t]
  # For EREN
  m_M_EREN[t + 1, ]  <- m_M_EREN[t, ]  %*% a_P_EREN[, , t]
}

## Store the cohort traces in a list 
l_m_M <- list(EPTI =  m_M_EPTI,
              EREN =  m_M_EREN)
names(l_m_M) <- v_names_str


### Plot the traces for both strategies
plot_trace(m_M_EPTI) 
plot_trace(m_M_EREN) 


##########################################
## OS Survival Curves
##########################################
### Print the Overall Survival - EPTI
v_os_EPTI <- 1 - m_M_EPTI[, "D"]    # calculate the overall survival (OS) probability
v_os_EPTI <- rowSums(m_M_EPTI[, 1:2])  # alternative way of calculating the OS probability   

plot(v_os_EPTI, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")  # create a simple plot showing the OS

# add grid 
grid(nx = n_cycles, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), 
     equilogs = TRUE) 


### Print the Overall Survival - EREN
v_os_EREN <- 1 - m_M_EREN[, "D"]    # calculate the overall survival (OS) probability
v_os_EREN <- rowSums(m_M_EREN[, 1:2])  # alternative way of calculating the OS probability   

plot(v_os_EREN, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")  # create a simple plot showing the OS

# add grid 
grid(nx = n_cycles, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), 
     equilogs = TRUE) 



### Life Expectancy
le_EPTI <- sum(v_os_EPTI)  # summing probability of OS over time  (i.e. life expectancy)
le_EPTI

le_EREN <- sum(v_os_EREN)  # summing probability of OS over time  (i.e. life expectancy)
le_EREN




####################################
## Prevalence of EM1
####################################
### EPTI
v_prev_EPTI <- m_M_EPTI[, "EM1"]/v_os_EPTI
plot(v_prev_EPTI,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

### EREN
v_prev_EREN <- m_M_EREN[, "EM1"]/v_os_EREN
plot(v_prev_EREN,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")




####################################
## State Payoffs
####################################
## Scale by the cycle length

# EPTI
# vector of state QALYs accrued each cycle
v_u_EPTI   <- c(EM0  = u.no_response_epti, 
                EM1  = u.response_epti, 
                D    = U_D) * cycle_length
# vector of state costs accrued each cycle
v_c_EPTI   <- c(EM0  = C_EM0, 
                EM1  = C_EM1 + c.eptinezumab + c.admin, 
                D    = C_D) * cycle_length


# EREN
# vector of state QALYs accrued each cycle
v_u_EREN    <- c(EM0  = u.no_response_eren, 
                 EM1  = u.response_eren,
                 D    = U_D) * cycle_length
# vector of state costs accrued each cycle
v_c_EREN    <- c(EM0  = C_EM0, 
                 EM1  = C_EM1 + c.erenumab,
                 D    = C_D) * cycle_length



## Store state rewards 
# Store the vectors of state utilities for each strategy in a list 
l_u   <- list(EPTI = v_u_EPTI,
              EREN = v_u_EREN)
l_u
# Store the vectors of state cost for each strategy in a list 
l_c   <- list(EPTI = v_c_EPTI,
              EREN = v_c_EREN)
l_c



####################################
## Compute expected outcomes
####################################

# Create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

## Loop through each strategy and calculate total utilities and costs 
for (i in 1:n_str) {
  v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
  v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
  
  ### Expected QALYs and costs per cycle 
  ## Vector of QALYs and Costs
  # Apply state rewards 
  v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
  v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
  
  ### Discounted total expected QALYs and Costs per strategy and apply within-cycle correction if applicable
  # QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
  # Costs
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
}



###################################
## CEA
###################################

## Incremental cost-effectiveness ratios (ICERs) 
df_cea <- calculate_icers(cost       = v_tot_cost + TC,    # Add the Decision Tree Costs
                          effect     = v_tot_qaly + TE,    # Add the Decision Tree Effects
                          strategies = v_names_str)
df_cea


## CEA table in proper format 
table_cea <- format_table_cea(df_cea) 
table_cea


plot(df_cea, label = "all", txtsize = 14) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.3))




########################################
## Probabilistic Sensitivity Analysis 
########################################
calculate_ce_out <- function(l_params_all, n_wtp = 50000, verbose = verbose){ # User defined
  with(as.list(l_params_all), {
    
    
    ##################################################################################################################
    ################################################ Decision tree ###################################################
    ##################################################################################################################
    
    ################################
    ## EM model
    ################################
    
    strategies <- c("EPTI", "EREN")
    
    
    ### Probabilities
    p.reduction_epti = p.reduction_epti             # Probability of a reduction of 50% or greater = EPTI
    p.no_reduction_epti = 1 - p.reduction_epti      # Probability of a reduction of < 50% - EPTI
    
    p.reduction_eren = p.reduction_eren             # Probability of a reduction of 50% or greater = EREN
    p.no_reduction_eren = 1 - p.reduction_eren      # Probability of a reduction of < 50% - EREN
    
    
    ### Drug costs
    c.eptinezumab = c.eptinezumab  # Cost of eptinezumab (1 month)
    c.admin = c.admin              # Cost of administering epti
    
    c.erenumab = c.erenumab        # Cost of erenumab (1 month)
    
    
    ### Utility scores
    u.response_epti = u.response_epti
    u.no_response_epti = u.no_response_epti
    
    u.response_eren = u.response_eren
    u.no_response_eren = u.no_response_eren
    
    ### Decision tree duration
    tree_length = tree_length         # Length of the decision tree (duration in months)
    
    
    
    ### Total Costs
    tc.epti_1 = tree_length*(c.eptinezumab + c.admin)      # Total costs for EPTI - branch 1
    tc.epti_2 = tree_length*(c.eptinezumab + c.admin)      # Total costs for EPTI - branch 2
    
    tc.eren_1 = tree_length*(c.erenumab)                   # Total costs for EREN - branch 1
    tc.eren_2 = tree_length*(c.erenumab)                   # Total costs for EREN - branch 2
    
    
    ### Total QALYs
    te.epti_1 = (u.response_epti / 2)      # EPTI - Total QALYs; divide by 2 for 6 months
    te.epti_2 = (u.no_response_epti / 2)   # EPTI - Total QALYs; divide by 2 for 6 months
    
    te.eren_1 = (u.response_eren / 2)      # EREN - Total QALYs; divide by 2 for 6 months
    te.eren_2 = (u.no_response_eren / 2)   # EREN - Total QALYs; divide by 2 for 6 months
    
    ### WTP
    wtp = wtp                            # Willingness-to-pay
    

    
    ##################################################
    #### OPEN TREE Method -- Decision Tree ####
    ##################################################
    
    ### Expected costs
    C.EPTI  <- c(prod(c(p.reduction_epti)),
                 prod(c(p.no_reduction_epti))) %*%
               c(tc.epti_1, 
                 tc.epti_2)
    
    C.EREN  <- c(prod(c(p.reduction_eren)),
                 prod(c(p.no_reduction_eren))) %*%
               c(tc.eren_1, 
                 tc.eren_2)
    
    ### Expected QALYs
    E.EPTI  <- c(prod(c(p.reduction_epti)),
                 prod(c(p.no_reduction_epti))) %*%
               c(te.epti_1, 
                 te.epti_2)
    
    E.EREN  <- c(prod(c(p.reduction_eren)),
                 prod(c(p.no_reduction_eren))) %*%
               c(te.eren_1, 
                 te.eren_2)
    
    ### Total costs and QALYs
    TC <- c(C.EPTI, C.EREN)
    TE <- c(E.EPTI, E.EREN)
    
    names(TC) <- names(TE) <- c("EPTI", "EREN")
    TC
    TE
    
    
    ### Incremental costs
    DC <- TC - TC[2]
    DC
    
    ### Incremental effects
    DE <- TE - TE[2]
    DE
    
    
    ### Incremental Cost-Effectiveness Ratio
    ICER <- DC / DE
    ICER
    
    
    ### Net Monetary Benefit
    NMB <- TE * wtp - TC
    NMB
    
    ### Create Full ICER table
    table.icer <- cbind(TC, TE, DC, DE, round(ICER, 2))
    table.icer <- as.data.frame(table.icer)
    colnames(table.icer) <- c("Costs", "QALYs", "Inc. Costs", "Inc. QALYs", "ICER")
    rownames(table.icer) <- strategies
    table.icer
    
    
    
    ##################################################################################################################
    ################################################ Markov model ####################################################
    ##################################################################################################################
    
    ## General setup 
    cycle_length <- 1/12   # cycle length equal to one year (use 1/12 for monthly)
    n_age_init   <- 40     # age at baseline
    n_age_max    <- 50    # maximum age of follow up
    n_cycles     <- (n_age_max - n_age_init)/cycle_length   # time horizon, number of cycles
    
    
    # Age labels 
    v_age_names  <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                          1:(1/cycle_length), 
                          sep = ".")
    
    # the 3 health states of the model:
    v_names_cycles  <- paste("cycle", 0:n_cycles)    # cycle names
    
    v_names_states <- c("EM0",  # EM off preventative treatment (EM0)
                        "EM1",  # EM on preventative treatment (EM1)
                        "D")    # Death (D)
    
    n_states <- length(v_names_states)   # number of health states 
    
    
    ### Discounting factors 
    d_c <- d_e <- 0.03  # annual discount rate for costs and QALY
    
    
    ### Strategies 
    v_names_str <- c("Eptinezumab",      # store the strategy names
                     "Erenumab") 
    n_str       <- length(v_names_str)   # number of strategies
    
    
    ## Within-cycle correction (WCC) using Simpson's 1/3 rule 
    v_wcc  <- gen_wcc(n_cycles = n_cycles, method = "Simpson1/3")
    
    
    ### Transition probabilities
    tp.EM0_EM1      <- 0
    tp.EM0_D        <- 0.0002826
    tp.EM0_EM0      <- 1 - (tp.EM0_EM1 + tp.EM0_D)
    
    
    tp.EM1_EM0_EREN <- tp.EM1_EM0_EREN  # Discontinuation due to AE
    tp.EM1_D        <- 0.0002826
    tp.EM1_EM1      <- 1 - (tp.EM1_EM0_EREN + tp.EM1_D)
    
    tp.D            <- 1 
    
    tp.EM1_EM0_EPTI <- tp.EM1_EM0_EPTI  # Discontinuation due to AE
    
    
    ### Age-dependent mortality rates (Use the life table)
    v_r_mort_by_age <- life_table %>% 
      dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
      dplyr::select(Total) %>%
      as.matrix()
    
    
    ### Markov model parameters
    p.female = 0.843             # Proportion of population that are female
    p.male = 1 - p.female        # Proportion of population that are male
    
    c.symptom_relief = 9.10      # Cost of symptom relief
    
    p.office_visit_m = 0.031     # Probability of migraine-associated office-visit (male)
    p.office_visit_f = 0.064     # Probability of migraine-associated office-visit (female)
    c.office_visit = 153.33         # Cost of migraine-associated office-visit
    
    p.hospitalization_m = 0.003  # Probability of migraine-associated hospitalization (male)
    p.hospitalization_f = 0.002  # Probability of migraine-associated hospitalization (female)
    c.hospitalization = 7364     # Cost of migraine-associated hospitalization
    
    p.ed_visit_m = 0.007         # Probability of migraine-associated ED visit (male)
    p.ed_visit_f = 0.013         # Probability of migraine-associated ED visit (female)
    c.ed_visit = 1114.55         # Cost if migraine-associated ED visit
    
    ### State costs
    C_EM0 <- p.female*((p.office_visit_f*c.office_visit) + (p.hospitalization_f*c.hospitalization) + (p.ed_visit_f*c.ed_visit)) + p.male*((p.office_visit_m*c.office_visit) + (p.hospitalization_m*c.hospitalization) + (p.ed_visit_m*c.ed_visit)) + c.symptom_relief
    
    C_EM1 <- p.female*((p.office_visit_f*c.office_visit) + (p.hospitalization_f*c.hospitalization) + (p.ed_visit_f*c.ed_visit)) + p.male*((p.office_visit_m*c.office_visit) + (p.hospitalization_m*c.hospitalization) + (p.ed_visit_m*c.ed_visit)) 
    
    C_D <- 0
    
    
    ### State utilities - EPTI
    U_EM0_epti <-  u.no_response_epti / 12    # Adjusted from 12 months to 1 month
    
    U_EM1_epti <-  u.response_epti / 12       # Adjusted from 12 months to 1 month
    
    U_D <- 0
    
    
    ### State utilities - EREN
    U_EM0_eren <-  u.no_response_eren / 12    # Adjusted from 12 months to 1 month
    
    U_EM1_eren <-  u.response_eren / 12       # Adjusted from 12 months to 1 month
    
    U_D <- 0
    
    
    ### Discount weight for costs and effects 
    v_dwc     <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
    v_dwe     <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
    
    
    # Process model inputs 
    ## Age-specific transition rates to the Dead state for all cycles 
    v_r_Dage  <- rep(v_r_mort_by_age, each = 1/cycle_length)
    # Name age-specific mortality vector 
    names(v_r_Dage) <- v_age_names
    
    ### Rate to probability
    v_p_Dage   <- rate_to_prob(v_r_Dage,  t = cycle_length) # Age-specific mortality risk All states to Death 
    
    # All starting healthy
    v_m_init_EPTI <- c("EM0" = p.no_reduction_epti, "EM1" = p.reduction_epti, "D" = 0)   # Change based on decision tree
    v_m_init_EPTI
    
    v_m_init_EREN <- c("EM0" = p.no_reduction_eren, "EM1" = p.reduction_eren, "D" = 0)   # Change based on decision tree
    v_m_init_EREN
    
    
    ### Initialize cohort trace for EREN and EPTI (These will be different based on the decision tree)
    m_M_EPTI <- matrix(0, 
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(v_names_cycles, v_names_states))
    
  
    m_M_EREN <- matrix(0, 
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(v_names_cycles, v_names_states))
    

    # Store the initial state vector in the first row of the cohort trace
    m_M_EPTI[1, ] <- v_m_init_EPTI
    m_M_EREN[1, ] <- v_m_init_EREN
    
    
    ## Create transition probability arrays for strategy EREN  
    ### Initialize transition probability array for strategy EREN 
    # All transitions to a non-death state are assumed to be conditional on survival
    a_P_EPTI <- array(0,  # Create 3-D array
                      dim = c(n_states, n_states, n_cycles),
                      dimnames = list(v_names_states, v_names_states, 
                                      v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 
    
    a_P_EREN <- array(0,  # Create 3-D array
                      dim = c(n_states, n_states, n_cycles),
                      dimnames = list(v_names_states, v_names_states, 
                                      v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 
    

    ### Fill in array
    ## EPTI
    # From EM0
    a_P_EPTI["EM0", "EM0", ] <- (1 - v_p_Dage) * (1 - tp.EM0_EM1)
    a_P_EPTI["EM0", "EM1", ] <- (1 - v_p_Dage) *      tp.EM0_EM1
    a_P_EPTI["EM0", "D", ]   <-      v_p_Dage
    
    # From EM1
    a_P_EPTI["EM1", "EM1", ] <- (1 - v_p_Dage) * (1 - tp.EM1_EM0_EPTI)
    a_P_EPTI["EM1", "EM0", ] <- (1 - v_p_Dage) *      tp.EM1_EM0_EPTI
    a_P_EPTI["EM1", "D", ]   <-      v_p_Dage
    
    # From Dead
    a_P_EPTI["D", "D", ] <- 1
        
    ## EREN
    # From EM0
    a_P_EREN["EM0", "EM0", ] <- (1 - v_p_Dage) * (1 - tp.EM0_EM1)
    a_P_EREN["EM0", "EM1", ] <- (1 - v_p_Dage) *      tp.EM0_EM1
    a_P_EREN["EM0", "D", ]   <-      v_p_Dage
    
    # From EM1
    a_P_EREN["EM1", "EM1", ] <- (1 - v_p_Dage) * (1 - tp.EM1_EM0_EREN)
    a_P_EREN["EM1", "EM0", ] <- (1 - v_p_Dage) *      tp.EM1_EM0_EREN
    a_P_EREN["EM1", "D", ]   <-      v_p_Dage
    
    # From Dead
    a_P_EREN["D", "D", ] <- 1
    

    
    ### Check if transition array and probabilities are valid
    # Check that transition probabilities are in [0, 1]
    check_transition_probability(a_P_EPTI, verbose = TRUE)
    check_transition_probability(a_P_EREN,  verbose = TRUE)
    
    # Check that all rows sum to 1
    check_sum_of_transition_array(a_P_EPTI, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    check_sum_of_transition_array(a_P_EREN, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    
    
    ###########################################
    ### Run the Markov model
    ###########################################
    # Iterative solution of age-dependent cSTM
    for(t in 1:n_cycles){
      ## Fill in cohort trace
      # For EPTI
      m_M_EPTI[t + 1, ]  <- m_M_EPTI[t, ]  %*% a_P_EPTI[, , t]      
      # For EREN
      m_M_EREN[t + 1, ]  <- m_M_EREN[t, ]  %*% a_P_EREN[, , t]
    }
    
    ## Store the cohort traces in a list 
    l_m_M <- list(EPTI =  m_M_EPTI,
                  EREN =  m_M_EREN)
    names(l_m_M) <- v_names_str
    

    
    ### Life Expectancy
    le_EPTI <- sum(v_os_EPTI)  # summing probability of OS over time  (i.e. life expectancy)
    le_EPTI
    
    le_EREN <- sum(v_os_EREN)  # summing probability of OS over time  (i.e. life expectancy)
    le_EREN
    

    
    ####################################
    ## State Payoffs
    ####################################
    ## Scale by the cycle length

    # EPTI
    # vector of state QALYs accrued each cycle
    v_u_EPTI   <- c(EM0  = u.no_response_epti, 
                    EM1  = u.response_epti, 
                    D    = U_D) * cycle_length
    # vector of state costs accrued each cycle
    v_c_EPTI   <- c(EM0  = C_EM0, 
                    EM1  = C_EM1 + c.eptinezumab + c.admin, 
                    D    = C_D) * cycle_length
    
        
    # EREN
    # vector of state QALYs accrued each cycle
    v_u_EREN    <- c(EM0  = u.no_response_eren, 
                     EM1  = u.response_eren,
                     D    = U_D) * cycle_length
    # vector of state costs accrued each cycle
    v_c_EREN    <- c(EM0  = C_EM0, 
                     EM1  = C_EM1 + c.erenumab,
                     D    = C_D) * cycle_length
    
    
    
    ## Store state rewards 
    # Store the vectors of state cost for each strategy in a list 
    l_c   <- list(EPTI = v_c_EPTI,
                  EREN = v_c_EREN)
    l_c
    
    # Store the vectors of state utilities for each strategy in a list 
    l_u   <- list(EPTI = v_u_EPTI,
                  EREN = v_u_EREN)
    l_u

    
    
    ####################################
    ## Compute expected outcomes
    ####################################
    
    # Create empty vectors to store total utilities and costs 
    v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
    names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
    
    ## Loop through each strategy and calculate total utilities and costs 
    for (i in 1:n_str) {
      v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
      v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
      
      ### Expected QALYs and costs per cycle 
      ## Vector of QALYs and Costs
      # Apply state rewards 
      v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
      v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
      
      ### Discounted total expected QALYs and Costs per strategy and apply within-cycle correction if applicable
      # QALYs
      v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
      # Costs
      v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
    }
    
    
    
    ###################################
    ## CEA
    ###################################

    ## Incremental cost-effectiveness ratios (ICERs) 
    df_cea <- calculate_icers(cost       = v_tot_cost + TC,    # Add the Decision Tree Costs
                              effect     = v_tot_qaly + TE,    # Add the Decision Tree Effects
                              strategies = v_names_str)
    return(df_cea)
  }
  )
}

    
    



##############################
## Parameters - PSA
##############################

l_params_all <- list(
  ### Probabilities
  p.reduction_epti = 0.66,                         # Probability of a reduction of 50% or greater = EPTI
  p.no_reduction_epti = 1 - p.reduction_epti,      # Probability of a reduction of < 50% - EPTI
  
  p.reduction_eren = 0.27,                         # Probability of a reduction of 50% or greater = EREN
  p.no_reduction_eren = 1 - p.reduction_eren,      # Probability of a reduction of < 50% - EREN
  
  
  ### Drug costs
  c.eptinezumab = 198.3867,  # Cost of eptinezumab (1 month)
  c.admin = 121.65,          # Cost of administering epti
  
  c.erenumab = 136.90,       # Cost of erenumab (1 month)
  
  
  ### Utility scores
  u.response_epti = 0.7938 + (-0.0189 * (8.7 - 3.9)),
  u.no_response_epti = 0.7938 + (-0.0189 * 8.7),
  
  u.response_eren = 0.7938 + (-0.0189 * (8.3 - 3.2)),
  u.no_response_eren = 0.7938 + (-0.0189 * 8.3),

  
  ### Decision tree duration
  tree_length = 6,           # Length of the decision tree (duration in months)
  
  ### Transition probabilities
  tp.EM1_EM0_EREN = 0.0055,  # Discontinuation due to AE
  tp.EM1_EM0_EPTI = 0.0091,  # Discontinuation due to AE 
  
  ### WTP
  wtp = 50000                            # Willingness-to-pay
)


calculate_ce_out(l_params_all = l_params_all )  



generate_psa_params <- function(n_sim = 1, seed = 12345){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(

    p.reduction_epti     = rbeta(n_sim, 5.9224, 3.050933),
    p.reduction_eren     = rbeta(n_sim, 5.0517, 13.6583),
    
    c.eptinezumab        = rgamma(n_sim, 4, 0.02016264),
    c.admin              = rgamma(n_sim, 4.89214, 0.04021488),
    
    c.erenumab           = rgamma(n_sim, 2.313779, 0.01690123),
    
    
    u.response_epti      = rbeta(n_sim, 9.489554, 4.007565),
    u.no_response_epti   = rbeta(n_sim, 11.5036, 6.774363),
    
    u.response_eren      = rbeta(n_sim, 9.523002, 4.13181),
    u.no_response_eren   = rbeta(n_sim, 11.53581, 6.575773),
    
    tp.EM1_EM0_EREN      = rbeta(n_sim, 3.9725, 718.3002),  
    tp.EM1_EM0_EPTI      = rbeta(n_sim, 5.119427, 557.4549)
  )
  return(df_psa)
}
    


# Store the parameter names into a vector
v_names_params <- names(l_params_all)

## Test functions to generate CE outcomes and PSA dataset 
# Test function to compute CE outcomes
calculate_ce_out(l_params_all) 

# Test function to generate PSA input dataset
generate_psa_params(10) 

## Generate PSA dataset 
# Number of simulations
n_sim <- 10000

# Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)

# First six observations
head(df_psa_input)

### Histogram of parameters 
ggplot(melt(df_psa_input, variable.name = "Parameter"), aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = ..density..)) +
  ylab("") +
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()) 



# Initialize data.frames with PSA output 
# data.frame of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str

# data.frame of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str

# Conduct probabilistic sensitivity analysis
# Run Markov model on each parameter set of PSA input dataset
n_time_init_psa_series <- Sys.time()

for (i in 1:n_sim) { # i <- 1
  l_psa_input <- update_param_list(l_params_all, df_psa_input[i,])
  
  # Outcomes
  l_out_ce_temp  <- calculate_ce_out(l_psa_input)
  df_c[i, ]  <- l_out_ce_temp$Cost  
  df_e[i, ]  <- l_out_ce_temp$Effect
# Display simulation progress
if (i/(n_sim/100) == round(i/(n_sim/100), 0)) { # display progress every 5%
  cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}
n_time_end_psa_series <- Sys.time()
n_time_total_psa_series <- n_time_end_psa_series - n_time_init_psa_series
print(paste0("PSA with ", scales::comma(n_sim), " simulations run in series in ", 
             round(n_time_total_psa_series, 2), " ", 
             units(n_time_total_psa_series)))



### Create PSA object 
l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = c("Erenumab", "Eptizenumab")) # Strategies are in reverse order

l_psa$strategies <- c("Erenumab", "Eptizenumab") # Strategies are in reverse order
colnames(l_psa$effectiveness) <- c("Erenumab", "Eptizenumab") # Strategies are in reverse order
colnames(l_psa$cost) <- c("Erenumab", "Eptizenumab") # Strategies are in reverse order

# Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 100000, by = 5000)


### Export to Excel to create PSA scatter
scatter_outcome <- cbind(df_c, df_e)
write.csv(scatter_outcome, "C:\\Users\\mbounthavong\\Dropbox\\Projects\\CEA Eptinezumab migraines\\PSA results\\psa_results.csv", row.names = TRUE)




#########################################
### ICER from PSA with 95% CrI
#########################################
## Use these functions, but it'll break make_psa_obj; need to reload "dampack" to fix.
##### Load functions needed for the PSA Credible Interval calculations
source("https://raw.githubusercontent.com/mbounthavong/Cost-Effectiveness-Analysis-of-Eptinezumab-and-Erenumab-for-Migraines/main/R%20Functions/icers.R")
source("https://raw.githubusercontent.com/mbounthavong/Cost-Effectiveness-Analysis-of-Eptinezumab-and-Erenumab-for-Migraines/main/R%20Functions/psa.R")
### ICER from PSA with 95% CrI
calculate_icers_psa(l_psa, uncertainty = TRUE)



### Cost-Effectiveness Scatter plot 
txtsize <- 13

gg_scattter <- plot_psa(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_scattter




### Incremental cost-effectiveness ratios (ICERs) with probabilistic output 
# Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)
df_cea_psa <- calculate_icers(cost       = df_out_ce_psa$meanCost, 
                              effect     = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa



### Plot cost-effectiveness frontier with probabilistic output 
plot_icers(df_cea_psa, label = "all", txtsize = txtsize) +
  expand_limits(x = max(df_cea_psa$Effect) + 0.1) +
  theme(legend.position = c(0.8, 0.3))




### Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) 
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)

# Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)

# CEAC & CEAF plot
gg_ceac <- plot_ceac(ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.25))
gg_ceac



### Expected Loss Curves (ELCs) 
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

# ELC plot
gg_elc <- plot_exp_loss(elc_obj, log_y = FALSE, 
                        txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
                        col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7),)
gg_elc



### Expected value of perfect information (EVPI) 
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
# EVPI plot
gg_evpi <- plot_evpi(evpi, effect_units = "QALY", 
                     txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000)
gg_evpi




