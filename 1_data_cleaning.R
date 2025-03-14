#########################################################################
# Replication file for 
# "Violence, Displacement, and Support for Internally Displaced Persons:
# Evidence from Syria"
# Authors: Alexandra Hartman, Benjamin Morse, Sigrid Weber
# Last update: 12/02/2021 
# Script: Data cleaning
#########################################################################

# This script cleans the data and constructs all variables for analysis

# load required packages -------------------------------------------------------
library(tidyverse)
library(lubridate)
# reshape2 is called in namespace

# read survey data -------------------------------------------------------------
data <- read_csv("data.csv")

# clean-up respondents socio-economic characteristics---------------------------

# Ethnicity (focusing on Arabs and Syriac-Assyrians)
data <- data %>% 
  mutate(ethnic= case_when(
    ethnic== "Arab (Except Palestinian)" ~ "Arab (Sunni Muslim)",
    ethnic== "Other" ~ "Arab (Sunni Muslim)", # contains only "Syrian" as reply
    ethnic== "Palestinian" ~ NA_character_,
    ethnic== "Turkmen"~ NA_character_,
    ethnic=="Syriac-Assyrian"~"Syriac/Assyrian (Christian)",
    TRUE ~ NA_character_))

# Respondents' gender
data <- data %>% 
  mutate(respondent_gender=
           ifelse(interviewee_type =="Most important female in household",
                  "Female",respondent_gender)) %>% 
  select(-interviewee_type)

# Age
data <- data %>%
  mutate(respondent_dob=
           replace(respondent_dob,respondent_dob==ymd(21890101),ymd(19890101)),
         age=floor(decimal_date(ymd(20171228)) - decimal_date(respondent_dob)),
         age=ifelse(age<=15,NA,age)) %>% 
  select(-respondent_dob)

# Job before the conflict
data <- data %>% 
  mutate(work_before_conflict=
           case_when(
             work_before_conflict== "Agriculture - farm worker" ~ "Agricultural occupation",
             work_before_conflict== "Agriculture - owned farm" ~ "Agricultural occupation",
             work_before_conflict== "Landlord" ~ "Small business owner",
             work_before_conflict== "Business owner" ~ "Small business owner",
             work_before_conflict== "Salaried work for civil society"~ "Salaried occupation",
             work_before_conflict== "Salaried work for government [including teachers, policemen, etc.)" ~ "Salaried occupation",
             work_before_conflict== "Salaried work in private company" ~ "Private sector employee",
             work_before_conflict== "Caring for family members" ~ "Domestic work",
             work_before_conflict== "Other" ~ NA_character_,
             TRUE ~ as.character(work_before_conflict)))

# Residence ownership, rural vs urban land, structure of housing
data <- data %>% 
  mutate(prev_res_who_own=
           case_when(
             prev_res_who_own=="Myself and someone else together"~"Own property",
             prev_res_who_own=="Someone I am unrelated to" ~ "Other private owner",
             prev_res_who_own=="Myself" ~ "Own property",
             prev_res_who_own=="Someone in my extended family" ~ "In family hands",
             prev_res_who_own=="Someone in my immediate family" ~ "In family hands",
             prev_res_who_own=="I don't know"~ NA_character_,
             prev_res_who_own=="An armed group"~ "Government/ conflict actor",
             prev_res_who_own=="The government"~ "Government/ conflict actor",
             TRUE ~as.character(prev_res_who_own)
           ))

# turn character variables into factors
data <-  data %>% mutate_if(sapply(., is.character), as.factor)

# Rural vs urban land and housing structure
levels(data$prev_res_land_type)<-c("Rural","Rural","Urban")
levels(data$prev_res_struc_type) <- c("Apartment",NA,"Standalone house","Tent")

# clean displacement and violence variables ------------------------------------
data <-data %>%
  mutate(displaced_from_prewar_residence=
           case_when(
             place_of_previous_residence == "No" ~ 1,
             place_of_previous_residence == "Yes" ~ 0,
             TRUE ~ NA_real_),
         displace_num_crisis=ifelse(displace_num_crisis > 55,1,displace_num_crisis),
         hh_deaths_6yrs=ifelse(hh_deaths_6yrs==2013,1,hh_deaths_6yrs),
         hh_deaths_6yrs_dum=ifelse(hh_deaths_6yrs>0,1,0),
         prewar_res_damaged_destroyed=case_when(
           prewar_residence_status == "Damaged but can be repaired" ~ 1,
           prewar_residence_status == "Destroyed" ~ 1,
           TRUE ~ 0),
         prewar_biz_damaged_destroyed=case_when(
           prewar_business_status == "Damaged but can be repaired" ~ 1,
           prewar_business_status == "Destroyed" ~ 1,
           TRUE ~ 0),
         # only count "non-natural deaths"
         hh_deaths_conflict = ifelse(hh_deaths_6yrs_dum==1&
                                       hh_deaths_age>=6&
                                       hh_deaths_age<=60,1,0),
         displaced_during_crisis=ifelse(displace_num_crisis>=1,1,0)
  )

#  clean- up hosting behaviour variables ---------------------------------------
data <- data %>% mutate(
  host_how_many_crisis=ifelse(host_how_many_crisis==150,NA,host_how_many_crisis),
  host_how_long=ifelse(host_how_long<0|host_how_long>360,0,host_how_long),
  hosted_number_duration = ifelse(host_how_long!=0,host_how_long*host_how_many_crisis,0)
)

# construct final violence index -----------------------------------------------

# index averages over 4 dummy variables:
# 1. Whether your houshold has experienced the death of a household member (age 6 to 60) in the last 6 years 
# 2. Whether your pre-war residence was damaged or destroyed
# 3. Whether your pre-war business was damaged or destroyed 
# 4. Whether you were displaced during the crisis

data <- data %>%
  mutate(violence_index= 
           rowMeans(data[,c("displaced_during_crisis",
                            "hh_deaths_conflict", 
                            "prewar_res_damaged_destroyed",
                            "prewar_biz_damaged_destroyed")]),
         violence_index_dum = ifelse(violence_index >= 
                                       mean(violence_index,na.rm = TRUE), 1,0))

# remove variables and store data for observational analysis -------------------
data <- data %>% select(-X1,-place_of_previous_residence, -hh_deaths_6yrs,
                        -hh_deaths_age, -displace_num_crisis,-prewar_residence_status,
                        -prewar_business_status, -displaced_from_prewar_residence,
                        -hh_deaths_6yrs_dum) 
obs_data <- data %>% select(-starts_with("photo")) # remove photo variables

# reshape photo decisions for conjoint analysis --------------------------------
data <- data%>%
  mutate(
    selected1=ifelse(photoa1==photoa_decision,1,0),
    selected2=ifelse(photoa2==photoa_decision,1,0),
    selected3=ifelse(photob1==photob_decision,1,0),
    selected4=ifelse(photob2==photob_decision,1,0),
    selected5=ifelse(photoc1==photoc_decision,1,0),
    selected6=ifelse(photoc2==photoc_decision,1,0)
  ) %>% select(-c(photoa_decision,photob_decision,photoc_decision))

photo <- reshape2::melt(data,
                        id.var="submission_id",   
                        measure.vars=c("photoa1",
                                       "photoa2",
                                       "photob1",
                                       "photob2",
                                       "photoc1",
                                       "photoc2"),
                        variable.name="profile",
                        value.name="photo")

data <- reshape2::melt(data,
                       measure.vars=c("selected1",  
                                      "selected2", 
                                      "selected3",
                                      "selected4",
                                      "selected5",
                                      "selected6"),
                       variable.name="profile",
                       value.name="idp_selected")

data <- data%>%
  mutate(photo=photo$photo)%>% filter(photo!="-")%>%
  mutate(photo=as.numeric(photo),
         idp_presented_first=ifelse(profile=="selected1"|
                                      profile=="selected3"|
                                      profile=="selected5",1,0))

# create dummies for IDP characteristics----------------------------------------
data <- data%>%
  mutate(idp_hh_single=ifelse(photo==1 |photo==3 |photo==5 |photo==7 |photo==9|
                                photo==11|photo==13|photo==15|photo==17|
                                photo==19|photo==21|photo==23|photo==25|
                                photo==27|photo==29|photo==31,1,0),
         idp_occ_farmer=ifelse(photo==3 |photo==4 |photo==7 |photo==8 | 
                                 photo==11|photo==12|photo==15|photo==16| 
                                 photo==19|photo==20|photo==23|photo==24| 
                                 photo==27|photo==28|photo==31|photo==32,1,0),
         idp_child_sick=ifelse(photo==5 |photo==6 |photo==7 |photo==8 |
                                 photo==13|photo==14|photo==15|photo==16|
                                 photo==21|photo==22|photo==23|photo==24| 
                                 photo==29|photo==30|photo==31|photo==32,1,0),
         idp_ethnic_outgroup=ifelse(photo==9 |photo==10|photo==11|photo==12|
                                      photo==13|photo==14|photo==15|photo==16| 
                                      photo==25|photo==26|photo==27|photo==28| 
                                      photo==29|photo==30|photo==31| photo==32,
                                    1,0),
         idp_religion_christian=ifelse(photo<17,1,0)  
  )

conjoint_data <- data %>% select(-starts_with("photo"),-profile)

# code variables for conjoint analysis -----------------------------------------

# remove anything but data for survey experiment and observational data --------
rm(list=setdiff(ls(), c("conjoint_data","obs_data")))

# write datafiles if necessary--------------------------------------------------
#write.csv(obs_data, "REPLICATION/obs_data.csv") 
#write.csv(conjoint_data,"REPLICATION/experimental_data.csv")
