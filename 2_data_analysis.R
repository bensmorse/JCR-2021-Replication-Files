#########################################################################
# Replication file for 
# "Violence, Displacement, and Support for Internally Displaced Persons:
# Evidence from Syria"
# Authors: Alexandra Hartman, Benjamin Morse, Sigrid Weber
# Last update: 15/02/2021
# Script: Analysis
#########################################################################

# This script analyzes the cleaned and anonymized data 
# (conjoint experiment and observational analysis)

# load required packages -------------------------------------------------------
library(estimatr)
library(tidyverse)
library(texreg)
library(stargazer)
library(dotwhisker)
library(sf)
library(tmap)

# read data if necessary -------------------------------------------------------
# obs_data <- read.csv("observational_data.csv",row.names=1)
# conjoint_data <- read.csv("experimental_data.csv",row.names=1)

# descriptive map of the sample areas in the North and South of Syria ----------

sampled_areas <- c("Harim","Atareb","Heish", "Dana", "A'zaz", "Kafr Takharim",
                   "Kafr Nobol", "Daret Azza", "Qourqeena", "Mare", 
                   "Ma'arrat An Nu'man", "Aghtrin", "Suran", "Nawa", "Dar'a",
                   "Da'el","Al-Khashniyyeh","Jizeh", "Hrak", "Mzeireb", "Jasim",
                   "Mseifra", "Quneitra", "Busra Esh-Sham", "Khan Arnaba", "Fiq",
                   "As-Sanamayn", "Kherbet Ghazala", "Mashnaf")

syria <- read_sf("syria_shape") %>% 
  mutate(sampled_areas = ifelse(NAME_EN %in% sampled_areas,1,0))

# sampled areas (figure 1 in manuscript)
tm_shape(syria)+
  tm_fill("sampled_areas", palette = c("gray90", "gray35"))+
  tm_borders()+
  tm_compass(position= c("left","top"))+
  tm_scale_bar(position=c("right","bottom"))+
  tm_legend(show=F)
  
        
# CONJOINT ANALYSIS ############################################################

# fit conjoint models ----------------------------------------------------------
m1 <- lm_robust(idp_selected ~  
                 idp_hh_single+
                 idp_occ_farmer+
                 idp_child_sick+
                 idp_ethnic_outgroup+
                 idp_religion_christian,
               clusters = submission_id,
               se_type = "stata",
               data = conjoint_data) # baseline model 

m2 <- conjoint_data %>%
  split(.$violence_index_dum) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id,
                   se_type = "stata",
                   data = .)) # effect of violence 

m3 <- 
  conjoint_data %>%
  split(list(.$ethnic,.$violence_index_dum)) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+ 
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id, 
                   se_type = "stata",
                   data = .)) # effect of violence splitted by ethnic group 

# display conjoint graphs ------------------------------------------------------
label <- c("Christian","Kurdish","Sick children","Farmer","Female HH")

# baseline of conjoint experiment (figure 2)
dwplot(m1,
       vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
       dot_args = list(size=3,color="black"),
       whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(text = element_text(size=15))

# conjoint experiment by experience of violence (figure 3)
dwplot(m2%>%
         map(tidy) %>%
         bind_rows(.id = "violence_index_dum")%>%
         mutate(model=ifelse(violence_index_dum==0,"Low violence","High violence")),
       vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
       dot_args = list(aes(shape = model),
                       size=3,color="black"),
       whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        text = element_text(size=15))

# conjoint experiment by experience of violence and ethnic group (figure 4)
dwplot(m3%>%
         map(tidy) %>%
         bind_rows(.id = "ethnic_violence_index_dum") %>%
         mutate(model=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                               ethnic_violence_index_dum=="Syriac/Assyrian (Christian).0",
                             "Low violence","High violence"),
                ethnic=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                                ethnic_violence_index_dum=="Arab (Sunni Muslim).1",
                              "Arab (Sunni Muslim)","Syriac/Assyrian (Christian)")),
       vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
       dot_args = list(aes(shape = model),
                       size=3,color="black"),
       whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  facet_wrap(.~ethnic)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))

# OBSERVATIONAL ANALYSIS #######################################################

# mean imputation of predictors ------------------------------------------------
obs_data <- obs_data %>% 
  mutate(violence_index=replace(violence_index, is.na(violence_index), mean(violence_index, na.rm=TRUE)),
         education_before_conflict=replace(education_before_conflict, is.na(education_before_conflict),"Primary school"),
         work_before_conflict=replace(work_before_conflict,is.na(work_before_conflict),"Salaried occupation"),
         prev_res_land_type=replace(prev_res_land_type,is.na(prev_res_land_type),"Urban"),
         prev_res_who_own=replace(prev_res_who_own,is.na(prev_res_who_own),"In family hands"),
         prev_res_struc_type=replace(prev_res_struc_type,is.na(prev_res_struc_type),"Standalone house"),
         ethnic=replace(ethnic, is.na(ethnic),"Arab (Sunni Muslim)"),
         respondent_gender=replace(respondent_gender,is.na(respondent_gender),"Male"),
         age=replace(age,is.na(age),mean(age, na.rm=TRUE)))

# relevel factors & z-transformation of violence index -------------------------
obs_data <- obs_data %>%
  mutate(work_before_conflict=relevel(factor(work_before_conflict),ref="Unemployed"),
         violence_index_std=(violence_index-mean(violence_index,na.rm=TRUE))/
      mean(violence_index,na.rm=TRUE))

# fit observational models ----------------------------------------------------- 
obs_data$currently_hosting <- ifelse(obs_data$currently_hosting=="Yes",1,0)

m4 <- lm_robust(currently_hosting ~ 
                  violence_index_std,
                clusters = current_village, 
                se_type ="stata",
                data = obs_data)

m5 <- lm_robust(currently_hosting~ 
                  violence_index_std+
                  education_before_conflict+
                  work_before_conflict+
                  prev_res_land_type+
                  prev_res_who_own+
                  prev_res_struc_type+
                  ethnic+
                  respondent_gender+
                  age,
                clusters = current_village,
                se_type ="stata",
                data = obs_data)

m6 <- lm_robust(hosted_number_duration ~ 
                  violence_index_std,
                clusters = current_village,
                se_type ="stata",
                data = obs_data)

m7 <- lm_robust(hosted_number_duration ~  
                  violence_index_std+
                  education_before_conflict+
                  work_before_conflict+
                  prev_res_land_type+
                  prev_res_who_own+
                  prev_res_struc_type+
                  ethnic+
                  respondent_gender+
                  age,
                clusters = current_village,
                se_type ="stata",
                data = obs_data)

# display regression output for observational analysis -------------------------
coef_names_obs <- c("Violence index (standardized)",
                    "Postgraduate degree","Primary school",
                    "Secondary school","University degree",
                    "Agricultural occupation","Domestic work",
                    "Informal work","Private sector employee", 
                    "Salaried occupation","Small business owner",
                    "Urban residence",
                    "In family hands","Other private owner","Own property",
                    "Standalone house","Tent",
                    "Syriac-Assyrian",
                    "Gender (Male)",
                    "Age")

# main results table for observational analysis (table 2)
screenreg(list(m4,m5,m6,m7),
          include.ci = FALSE,
          omit.coef="(Intercept)",
          custom.coef.names =coef_names_obs,
          scalebox = 0.8,
          float.pos="hptb!",
          custom.model.names = 
            c("Currently hosting IDPs", " ","Total # of months hosted (Family size x duration)"," "),
          caption = "Effect of the experience of violence on observable hosting behaviour"
)

# Benjamini & Hochber adjusted p values
p_values <- c(measure1= m4$p.value["violence_index_std"],
              measure1 = m5$p.value["violence_index_std"],
              measure2 = m6$p.value["violence_index_std"],
              measure2 =m7$p.value["violence_index_std"])

p_values_adjusted <- c(without_controls = p.adjust(p_values[c(1,3)],"BH"),
                       with_controls = p.adjust(p_values[c(2,4)],"BH"))

p_values_adjusted <- p_values_adjusted[c(1,3,2,4)] # bring in right order
p_values_adjusted 

  
# FIGURES AND TABLES FOR APPENDIX ##############################################

# Figure A.1: Distribution of exposure to violence in data
distribution1 <- ggplot(obs_data %>% filter(violence_index!=mean(violence_index)))+
  geom_bar(aes(4*violence_index), col = "gray50", fill ="lightgray")+
  theme_bw()+
  xlab("Index of exposure to violence (from 0 to 4)")+
  ylab("Count of respondents")+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))

distribution2 <- ggplot(obs_data %>% filter(is.na(violence_index_dum)==F))+
  geom_bar(aes(factor(violence_index_dum)), col = "gray50", fill ="lightgray")+
  theme_bw()+
  xlab("Dummy of exposure to violence")+
  ylab("Count of respondents")+ 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("Below average", "Above average"))+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))

gridExtra::grid.arrange(distribution1, distribution2, ncol=2)

# Figure A.2. Distribution of exposure to violence per ethnic group
distribution1_ethnic <- ggplot(obs_data %>% filter(violence_index!=mean(violence_index)))+
  geom_bar(aes(y = ..prop..,4*violence_index), col = "gray50", fill ="lightgray")+
  theme_bw()+
  xlab("Index of exposure to violence (from 0 to 4)")+
  ylab("Frequency")+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))+
  facet_wrap(~ethnic)

distribution2_ethnic<- ggplot(obs_data)+
  geom_bar(aes(y = ..prop..,violence_index_dum), col = "gray50", fill ="lightgray")+
  theme_bw()+
  scale_x_continuous(breaks=c(0,1),
                     labels = c("Below average","Above average"))+
  xlab("Dummy of exposure to violence")+
  ylab("Frequency")+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))+
  facet_wrap(~ethnic)

gridExtra::grid.arrange(distribution1_ethnic, distribution2_ethnic, nrow=2)

# Table A.1: Descriptive summary statistics ------------------------------------

summarise_variables <- function(x, na.rm=TRUE){
  result <- c(Mean=mean(x, na.rm=na.rm),
              SD=sd(x, na.rm=na.rm),
              Min=min(x, na.rm=na.rm),
              Max=max(x, na.rm=na.rm),
              Missing=sum(is.na(x))
  )
}

df <- obs_data[, c("displaced_during_crisis","hh_deaths_conflict",
                   "prewar_res_damaged_destroyed","prewar_biz_damaged_destroyed",
                   "violence_index","violence_index_std","violence_index_dum",
                   "currently_hosting","hosted_number_duration","age")]

covariates <- data.frame(
  model.matrix( ~ education_before_conflict - 1 + 
                  work_before_conflict - 1 +
                  prev_res_land_type - 1 +                        
                  prev_res_who_own - 1 + 
                  prev_res_struc_type - 1 + 
                  ethnic - 1 + 
                  respondent_gender - 1, 
                data=obs_data))

df<-sapply(df,summarise_variables)
covariates <- sapply(covariates,summarise_variables)
df <- cbind(df, covariates)

colnames(df) <- c("Experienced displacement",
                  "Experienced death in household",
                  "Experienced residence destruction",
                  "Experienced business destruction", 
                  "Violence index",
                  "Violence index (standardized)",
                  "Experienced more violence than mean",
                  "Hosted any IDP",
                  "Months hosted x family size",
                  "Age","No education before conflict",
                  "Post-graduate degree before conflict",
                  "Primary school before conflict",
                  "Secondary school before conflict",
                  "University degree before conflict",
                  "Agricultural occupation","Domestic work",
                  "Informal work","Private sector emploee", 
                  "Salaried occupation","Small business owner",
                  "Urban residence","Residence in family hands",
                  "Residence owned by private owner",
                  "Residence owned in own hands","Standalone house",
                  "Tent", "Syriac-Assyrian","Male")
stargazer::stargazer(df,header=FALSE,flip = TRUE,digits=2, type ="text")

# Table A.2: Regression table for conjoint (main results in paper) -------------

m1_with_interaction <- conjoint_data %>%
  lm_robust(idp_selected ~
              idp_hh_single +
              idp_occ_farmer +
              idp_child_sick +
              idp_ethnic_outgroup +
              idp_religion_christian +
              idp_hh_single:violence_index_dum +
              idp_occ_farmer:violence_index_dum +
              idp_child_sick:violence_index_dum +
              idp_ethnic_outgroup:violence_index_dum +
              idp_religion_christian:violence_index_dum,
            clusters = submission_id,
            se_type = "stata",
            data = .)

m1_with_interaction_and_baseterm <- conjoint_data %>%
  lm_robust(idp_selected ~
              idp_hh_single +
              idp_occ_farmer +
              idp_child_sick +
              idp_ethnic_outgroup +
              idp_religion_christian +
              idp_hh_single:violence_index_dum +
              idp_occ_farmer:violence_index_dum +
              idp_child_sick:violence_index_dum +
              idp_ethnic_outgroup:violence_index_dum +
              idp_religion_christian:violence_index_dum+
              violence_index_dum,
            clusters = submission_id,
            se_type = "stata",
            data = .)

m_with_interaction_covariates <- conjoint_data %>%
  lm_robust(idp_selected ~
              idp_hh_single +
              idp_occ_farmer +
              idp_child_sick +
              idp_ethnic_outgroup +
              idp_religion_christian +
              
              idp_hh_single:violence_index_dum +
              idp_occ_farmer:violence_index_dum +
              idp_child_sick:violence_index_dum +
              idp_ethnic_outgroup:violence_index_dum +
              idp_religion_christian:violence_index_dum +
              
              idp_hh_single:respondent_gender + 
              idp_occ_farmer:respondent_gender + 
              idp_child_sick:respondent_gender + 
              idp_ethnic_outgroup:respondent_gender + 
              idp_religion_christian:respondent_gender +
              
              idp_hh_single:prev_res_struc_type + 
              idp_occ_farmer:prev_res_struc_type + 
              idp_child_sick:prev_res_struc_type +
              idp_ethnic_outgroup:prev_res_struc_type + 
              idp_religion_christian:prev_res_struc_type+
              
              idp_hh_single:prev_res_who_own + 
              idp_occ_farmer:prev_res_who_own + 
              idp_child_sick:prev_res_who_own +
              idp_ethnic_outgroup:prev_res_who_own + 
              idp_religion_christian:prev_res_who_own +
              
              idp_hh_single:education_before_conflict + 
              idp_occ_farmer:education_before_conflict + 
              idp_child_sick:education_before_conflict + 
              idp_ethnic_outgroup:education_before_conflict + 
              idp_religion_christian:education_before_conflict +
              
              idp_hh_single:work_before_conflict + 
              idp_occ_farmer:work_before_conflict + 
              idp_child_sick:work_before_conflict + 
              idp_ethnic_outgroup:work_before_conflict + 
              idp_religion_christian:work_before_conflict +
              
              idp_hh_single:age + 
              idp_occ_farmer:age + 
              idp_child_sick:age + 
              idp_ethnic_outgroup:age + 
              idp_religion_christian:age,
            clusters = submission_id,
            se_type = "stata",
            data = .)

conjoint_with_interaction <- c("Female HH","Farmer","Sick child",
                               "Kurdish","Christian", 
                               "Female HH : High violence",
                               "Farmer : High violence",
                               "Sick child : High violence",
                               "Kurdish : High violence",
                               "Christian : High violence",
                               "High violence")


screenreg(list(m1,m2[[1]],m2[[2]], 
               m1_with_interaction, 
               m1_with_interaction_and_baseterm,
               m_with_interaction_covariates),
          include.ci = FALSE,
          omit.coef="(Intercept)|(work)|(age)|(education)|(prev)|(gender)",
          custom.coef.names = conjoint_with_interaction,
          custom.model.names=c("Baseline model",
                               "Low violence",
                               "High violence",
                               "Interaction", 
                               "Interaction (+ Base term)",
                               "Full data with interaction (+ Interaction with Control)"),
          scalebox = 0.6,
          float.pos="hptb!",
          caption = "Full regression results for conjoint analysis")


# Table A.3: Full results of observational analysis + additional models --------

# first estimate models for number of hosted and duration of hosted separately
number_hosted <- lm_robust(ifelse(is.na(host_how_many_crisis)==T,0,host_how_many_crisis) ~  
                             violence_index_std+
                             education_before_conflict+
                             work_before_conflict+
                             prev_res_land_type+
                             prev_res_who_own+
                             prev_res_struc_type+
                             ethnic+
                             respondent_gender+
                             age,
                           clusters = current_village,
                           se_type ="stata",
                           data = obs_data)

duration_hosted <- lm_robust(ifelse(is.na(host_how_long)==T,0,host_how_long) ~  
                               violence_index_std+
                               education_before_conflict+
                               work_before_conflict+
                               prev_res_land_type+
                               prev_res_who_own+
                               prev_res_struc_type+
                               ethnic+
                               respondent_gender+
                               age,
                             clusters = current_village,
                             se_type ="stata",
                             data = obs_data)

# get mean vale of outcome variables
means <- c(round(mean(model.frame(m4)$currently_hosting),3),
           round(mean(model.frame(m5)$currently_hosting),3),
           round(mean(model.frame(m6)$hosted_number_duration),3),
           round(mean(model.frame(m7)$hosted_number_duration),3),
           round(mean(model.frame(number_hosted)$"ifelse(is.na(host_how_many_crisis) == T, 0, host_how_many_crisis)"),3),
           round(mean(model.frame(duration_hosted)$"ifelse(is.na(host_how_long) == T, 0, host_how_long)",na.rm=T),3))

# then add this to the main observational models presented in the paper
screenreg(list(m4,m5,m6,m7, number_hosted,duration_hosted),
       include.ci = FALSE,
       omit.coef="(Intercept)",
       custom.coef.names = coef_names_obs,
       scalebox = 0.8,
       float.pos="hptb!",
       custom.model.names = c("Currently hosting"," ",
                              "Total \\# of months hosted","",
                              "Number of hosted IDPs", 
                              "Duration of hosting"),
       caption = "Effect of violence on number of IDPs hosted and 
          duration of IDP hosting")

# Table A.4: Selection into violence -------------------------------------------
m8 <- lm_robust(violence_index_dum ~ 
                  education_before_conflict+
                  work_before_conflict+
                  prev_res_land_type+
                  prev_res_who_own+
                  prev_res_struc_type+
                  ethnic+ 
                  respondent_gender+
                  age,
                clusters = submission_id,
                se_type = "stata",
                data = obs_data)


m9 <- lm_robust(violence_index_std ~  
                  education_before_conflict+
                  work_before_conflict+
                  prev_res_land_type+
                  prev_res_who_own+
                  prev_res_struc_type+
                  ethnic+ 
                  respondent_gender+
                  age,
                clusters = submission_id,
                se_type = "stata",
                data = obs_data)

coef_names_balance <- c("Postgraduate degree","Primary school",
                        "Secondary school","University degree",
                        "Agricultural occupation","Domestic work",
                        "Informal work","Private sector employee", 
                        "Salaried occupation","Small business owner",
                        "Urban residence","In family hands",
                        "Other private owner","Own property",
                        "Standalone house","Tent",
                        "Syriac-Assyrian","Gender (Male)","Age")

screenreg(list(m8,m9),
          include.ci = FALSE,
          omit.coef="(Intercept)",
          scalebox = 0.8,
          custom.coef.names =coef_names_balance, 
          float.pos="hptb!", 
          custom.model.names=c("Violence (binary)","Violence (continuous)" ),
          caption="Balance regression: prediction of violence"
)

# Sensitivity analysis for unobserved confounders: see stata script ------------

# the dta file used for the analysis in stata can be written with these commands:
# library(foreign)
# write.dta(obs_data,"stata_file.dta")


# Figure A.3+4 Reduced violence index without displacement ---------------------
conjoint_data <- conjoint_data%>%
  mutate(violence_reduced_index = rowMeans(
    .[,c("hh_deaths_conflict",
         "prewar_res_damaged_destroyed",
         "prewar_biz_damaged_destroyed")]),
    violence_reduced_std = 
      (violence_reduced_index-mean(violence_reduced_index,na.rm=TRUE))/
      mean(violence_reduced_index,na.rm=TRUE))

conjoint_data <- conjoint_data %>% 
  mutate(violence_reduced_dum = ifelse(violence_reduced_index >= 
                                         mean(violence_reduced_index,na.rm = TRUE), 1,0))

m2_reduced <- conjoint_data %>%
  split(.$violence_reduced_dum) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id,
                   se_type = "stata",
                   data = .)) # effect of violence 

m2_displacement <- conjoint_data %>%
  split(.$displaced_during_crisis) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id,
                   se_type = "stata",
                   data = .)) # effect of violence 
m3_reduced <- 
  conjoint_data %>%
  split(list(.$ethnic,.$violence_reduced_dum)) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+ 
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id, 
                   se_type = "stata",
                   data = .)) # effect of violence splitted by ethnic group 

m3_displacement <- 
  conjoint_data %>%
  split(list(.$ethnic,.$displaced_during_crisis)) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+ 
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id, 
                   se_type = "stata",
                   data = .)) # effect of violence splitted by ethnic group 


m2_reduced <- dwplot(m2_reduced%>%
                       map(tidy) %>%
                       bind_rows(.id = "violence_reduced_dum")%>%
                       mutate(model=ifelse(violence_reduced_dum==0,"Low violence","High violence")),
                     vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
                     dot_args = list(aes(shape = model),
                                     size=3,color="black"),
                     whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        text = element_text(size=15))

m2_displacement <- dwplot(m2_displacement%>%
                            map(tidy) %>%
                            bind_rows(.id = "displaced_during_crisis")%>%
                            mutate(model=ifelse(displaced_during_crisis==0,"No displacement","Displacement")),
                          vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
                          dot_args = list(aes(shape = model),
                                          size=3,color="black"),
                          whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        text = element_text(size=15))

gridExtra::grid.arrange(m2_reduced,m2_displacement, nrow=1)

figure3_reduced <- dwplot(m3_reduced%>%
                            map(tidy) %>%
                            bind_rows(.id = "ethnic_violence_index_dum") %>%
                            mutate(model=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                                                  ethnic_violence_index_dum=="Syriac/Assyrian (Christian).0",
                                                "Low violence","High violence"),
                                   ethnic=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                                                   ethnic_violence_index_dum=="Arab (Sunni Muslim).1",
                                                 "Arab (Sunni Muslim)","Syriac/Assyrian (Christian)")),
                          vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
                          dot_args = list(aes(shape = model),
                                          size=3,color="black"),
                          whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  facet_wrap(.~ethnic)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))

figure3_displacement <- dwplot(m3_displacement%>%
                                 map(tidy) %>%
                                 bind_rows(.id = "ethnic_violence_index_dum") %>%
                                 mutate(model=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                                                       ethnic_violence_index_dum=="Syriac/Assyrian (Christian).0",
                                                     "No displacement","Displacement"),
                                        ethnic=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                                                        ethnic_violence_index_dum=="Arab (Sunni Muslim).1",
                                                      "Arab (Sunni Muslim)","Syriac/Assyrian (Christian)")),
                               vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
                               dot_args = list(aes(shape = model),
                                               size=3,color="black"),
                               whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  facet_wrap(.~ethnic)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))

gridExtra::grid.arrange(figure3_reduced,figure3_displacement)

# Table A.6: Alternative violence index without displacement ------------------

obs_data <- obs_data%>%
  mutate(violence_reduced_index = rowMeans(
    .[,c("hh_deaths_conflict",
         "prewar_res_damaged_destroyed",
         "prewar_biz_damaged_destroyed")]),
    violence_reduced_std = 
      (violence_reduced_index-mean(violence_reduced_index,na.rm=TRUE))/
      mean(violence_reduced_index,na.rm=TRUE))

m10<- lm_robust(currently_hosting ~ 
                  violence_reduced_std,
                clusters = current_village,
                se_type ="stata",
                data = obs_data)

m11 <- lm_robust(currently_hosting ~ 
                   violence_reduced_std+
                   education_before_conflict+ 
                   work_before_conflict+
                   prev_res_land_type+ 
                   prev_res_who_own+ 
                   prev_res_struc_type+ 
                   ethnic+
                   respondent_gender+
                   age,
                 clusters = current_village,
                 se_type ="stata",
                 data = obs_data)


m12 <-lm_robust(hosted_number_duration ~ 
                  violence_reduced_std,
                clusters= current_village, 
                se_type="stata",
                data= obs_data)

m13 <- lm_robust(hosted_number_duration ~ 
                   violence_reduced_std+
                   education_before_conflict+
                   work_before_conflict+
                   prev_res_land_type+
                   prev_res_who_own+
                   prev_res_struc_type+
                   ethnic+
                   respondent_gender+
                   age,
                 clusters = current_village, 
                 se_type ="stata", 
                 data = obs_data)

m14 <- lm_robust(currently_hosting ~ 
                   displaced_during_crisis,
                 clusters = current_village,
                 se_type ="stata",
                 data = obs_data)                             

m15 <- lm_robust(currently_hosting ~ 
                   displaced_during_crisis+
                   education_before_conflict+
                   work_before_conflict+
                   prev_res_land_type+
                   prev_res_who_own+ 
                   prev_res_struc_type+
                   ethnic+ 
                   respondent_gender+
                   age,
                 clusters= current_village,
                 se_type="stata",
                 data= obs_data)


m16 <-lm_robust(hosted_number_duration ~ 
                  displaced_during_crisis,
                clusters = current_village,
                se_type ="stata",
                data = obs_data)

m17 <- lm_robust(hosted_number_duration ~ 
                   displaced_during_crisis+
                   education_before_conflict+
                   work_before_conflict+
                   prev_res_land_type+
                   prev_res_who_own+
                   prev_res_struc_type+
                   ethnic+
                   respondent_gender+ 
                   age,                        
                 clusters = current_village,
                 se_type ="stata",
                 data = obs_data)

coef_names_reduced_index <- c("Reduced violence index",
                              "Experience of displacement",
                              "Postgraduate degree","Primary school",
                              "Secondary school","University degree",
                              "Agricultural occupation","Domestic work",
                              "Informal work","Private sector employee", 
                              "Salaried occupation","Small business owner",
                              "Urban residence","In family hands",
                              "Other private owner","Own property",
                              "Standalone house","Tent","Syriac-Assyrian",
                              "Gender (Male)","Age")

screenreg(list(m10,m14,m11,m15,m12,m16,m13,m17),
          include.ci = FALSE,
          omit.coef="(Intercept)",
          custom.coef.names = coef_names_reduced_index,     
          scalebox = 0.6,
          float.pos="hptb!",
          caption = "Effect of violence (reduced index without displacement) and 
       displacement as separate predictor on binary hosting outcome (Models 1-4) 
       and on length x duration of IDP hosting (Models 5-8).")

# Figure A.5+6: Alternative split of violence index (binary) ---------------------------------
conjoint_data <- conjoint_data %>% 
  mutate(violence_index_dum_a = ifelse(displaced_during_crisis==1|
                                         hh_deaths_conflict==1|
                                         prewar_res_damaged_destroyed==1|
                                         prewar_biz_damaged_destroyed==1, 1,0))

m18 <- conjoint_data %>%
  split(.$violence_index_dum_a) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id,
                   se_type = "stata",
                   data = .)) # effect of violence 

m19 <- 
  conjoint_data %>%
  split(list(.$ethnic,.$violence_index_dum_a)) %>%
  map( ~ lm_robust(idp_selected ~
                     idp_hh_single+
                     idp_occ_farmer+
                     idp_child_sick+ 
                     idp_ethnic_outgroup+
                     idp_religion_christian,
                   clusters = submission_id, 
                   se_type = "stata",
                   data = .)) # effect of violence splitted by ethnic group 

dwplot(m18%>%
         map(tidy) %>%
         bind_rows(.id = "violence_index_dum")%>%
         mutate(model=ifelse(violence_index_dum==0,"Low violence","High violence")),
       vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
       dot_args = list(aes(shape = model),
                       size=3,color="black"),
       whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        text = element_text(size=15))

dwplot(m19%>%
         map(tidy) %>%
         bind_rows(.id = "ethnic_violence_index_dum") %>%
         mutate(model=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                               ethnic_violence_index_dum=="Syriac/Assyrian (Christian).0",
                             "Low violence","High violence"),
                ethnic=ifelse(ethnic_violence_index_dum=="Arab (Sunni Muslim).0"|
                                ethnic_violence_index_dum=="Arab (Sunni Muslim).1",
                              "Arab (Sunni Muslim)","Syriac/Assyrian (Christian)")),
       vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2),
       dot_args = list(aes(shape = model),
                       size=3,color="black"),
       whisker_args= list(color="black"))+
  scale_y_discrete(labels = label)+
  facet_wrap(.~ethnic)+
  ylab("")+
  xlab("Change in Pr(Host IDP)")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size=15),
        strip.text.x = element_text(size = 16))

