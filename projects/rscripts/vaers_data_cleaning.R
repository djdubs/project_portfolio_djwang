# Daniel Wang
# ST 563 Project
# 4/28/26
# Importing/Cleaning data

## This is the only client-side code to change
# Set file path to location of VAERS data for 2018-2021
setwd("projects/data")

library(tidyverse)

# ______________________________________________________________________________
# 2018-2021 VAERS data
vaers_2018 <- read.csv("2018VAERSDATA.csv")
vaers_2019 <- read.csv("2019VAERSDATA.csv")
vaers_2020 <- read.csv("2020VAERSDATA.csv")
vaers_2021 <- read.csv("2021VAERSDATA.csv")

# Each year contains the same columns, so we can concatenate each dataframe at once
vaers.2018_2021 <- rbind(vaers_2018, vaers_2019, vaers_2020, vaers_2021)

vaers_clean <- vaers.2018_2021 %>%
  arrange(VAERS_ID, desc(RECVDATE)) %>%
  distinct(VAERS_ID, .keep_all = T) %>%
  mutate(RECVDATE = as.Date(RECVDATE, format = "%m/%d/%Y")) %>%
  mutate(STATE = ifelse(STATE == "", NA, STATE)) %>%
  mutate(SEX = ifelse(SEX == "", NA, SEX)) %>%
  mutate(RPT_DATE = as.Date(RPT_DATE, format = "%m/%d/%Y")) %>%
  mutate(DIED = ifelse(DIED == "", 0, 1)) %>%
  mutate(DATEDIED = ifelse(DATEDIED == "", NA, DATEDIED)) %>%
  mutate(L_THREAT = ifelse(L_THREAT == "", 0, 1)) %>%
  mutate(ER_VISIT = ifelse(ER_VISIT == "", 0, 1)) %>%
  mutate(HOSPITAL = ifelse(HOSPITAL == "", 0, 1)) %>%
  mutate(X_STAY = ifelse(X_STAY == "", 0, 1)) %>%
  mutate(DISABLE = ifelse(DISABLE == "", 0, 1)) %>%
  mutate(RECOVD = case_when(
    RECOVD =="Y" ~ 1,
    RECOVD =="N" ~ 0,
    RECOVD =="U" ~ NA,
    RECOVD ==""  ~ NA
  )) %>%
  mutate(VAX_DATE = as.Date(VAX_DATE, format = "%m/%d/%Y")) %>%
  mutate(ONSET_DATE = as.Date(ONSET_DATE, format = "%m/%d/%Y")) %>%
  mutate(V_FUNDBY = ifelse(V_FUNDBY == "", NA, V_FUNDBY)) %>%
  mutate(SEVERITY = case_when(
    DIED == 1 ~ 5,
    L_THREAT == 1 ~ 4,
    HOSPITAL == 1 ~ 3,
    ER_VISIT == 1 | ER_ED_VISIT == "Y" ~ 2,
    OFC_VISIT == "Y" ~ 1,
    TRUE ~ 0)) %>%
  mutate(high_severity = ifelse(SEVERITY >= 3, 1, 0))

# only kept the most recent record for each ID
# Dates were converted to R's default format for dates: YYYY-MM-DD
# Any blank strings or unknown values, U, for any categorical variable were converted to an NA value
# Any blank strings for date variables were converted to an NA value
# Binary variables (Y/N) were converted to dummy variables to equal 1 or 0. 
# columns with character strings that do not represent a variable (i.e. SYMPTOM_TEXT) are untouched but temporarily kept in the dataframe
# columns included in VAERS form 2 are untouched (FORM_VERS and all following columns)
# Severity and high_severity represent the categorical and binary responses respectively
# currently unsure of the discrepancy between ER_VISIT and ER_ED_VISIT, they do not have consistent matching values
# ______________________________________________________________________________


# ______________________________________________________________________________
# 2018-2021 VAERS SYMPTOMS data
sym_2018 <- read.csv("2018VAERSSYMPTOMS.csv")
sym_2019 <- read.csv("2019VAERSSYMPTOMS.csv")
sym_2020 <- read.csv("2020VAERSSYMPTOMS.csv")
sym_2021 <- read.csv("2021VAERSSYMPTOMS.csv")

sym.2018_2021 <- rbind(sym_2018, sym_2019, sym_2020, sym_2021)

# symptoms of interest: pyrexia, dyspnoea, erythema, asthenia
sym_clean <- sym.2018_2021 %>%
  mutate(pyrexia = case_when(
    SYMPTOM1 == "Pyrexia" ~ 1,
    SYMPTOM2 == "Pyrexia" ~ 1,
    SYMPTOM3 == "Pyrexia" ~ 1,
    SYMPTOM4 == "Pyrexia" ~ 1,
    SYMPTOM5 == "Pyrexia" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(dyspnoea = case_when(
    SYMPTOM1 == "Dyspnoea" ~ 1,
    SYMPTOM2 == "Dyspnoea" ~ 1,
    SYMPTOM3 == "Dyspnoea" ~ 1,
    SYMPTOM4 == "Dyspnoea" ~ 1,
    SYMPTOM5 == "Dyspnoea" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(erythema = case_when(
    SYMPTOM1 == "Erythema" ~ 1,
    SYMPTOM2 == "Erythema" ~ 1,
    SYMPTOM3 == "Erythema" ~ 1,
    SYMPTOM4 == "Erythema" ~ 1,
    SYMPTOM5 == "Erythema" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(asthenia = case_when(
    SYMPTOM1 == "Asthenia" ~ 1,
    SYMPTOM2 == "Asthenia" ~ 1,
    SYMPTOM3 == "Asthenia" ~ 1,
    SYMPTOM4 == "Asthenia" ~ 1,
    SYMPTOM5 == "Asthenia" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(covid19 = case_when(
    SYMPTOM1 == "COVID-19" ~ 1,
    SYMPTOM2 == "COVID-19" ~ 1,
    SYMPTOM3 == "COVID-19" ~ 1,
    SYMPTOM4 == "COVID-19" ~ 1,
    SYMPTOM5 == "COVID-19" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(chest_pain = case_when(
    SYMPTOM1 == "Chest pain" ~ 1,
    SYMPTOM2 == "Chest pain" ~ 1,
    SYMPTOM3 == "Chest pain" ~ 1,
    SYMPTOM4 == "Chest pain" ~ 1,
    SYMPTOM5 == "Chest pain" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(dec_app = case_when(
    SYMPTOM1 == "Decreased appetite" ~ 1,
    SYMPTOM2 == "Decreased appetite" ~ 1,
    SYMPTOM3 == "Decreased appetite" ~ 1,
    SYMPTOM4 == "Decreased appetite" ~ 1,
    SYMPTOM5 == "Decreased appetite" ~ 1,
    TRUE ~ 0
  )) %>%
  group_by(VAERS_ID) %>%
  summarise(
    pyrexia = max(pyrexia),
    dyspnoea = max(dyspnoea),
    erythema = max(erythema),
    asthenia = max(asthenia),
    covid19 = max(covid19),
    chest_pain = max(chest_pain),
    dec_app = max(dec_app),
    .groups = "drop"
  )

# seven indicator variables for various symptoms
# ______________________________________________________________________________


# ______________________________________________________________________________
# 2018-2021 VAERS VAX data
vax_2018 <- read.csv("2018VAERSVAX.csv")
vax_2019 <- read.csv("2019VAERSVAX.csv")
vax_2020 <- read.csv("2020VAERSVAX.csv")
vax_2021 <- read.csv("2021VAERSVAX.csv")

vax.2018_2021 <- rbind(vax_2018, vax_2019, vax_2020, vax_2021)

vax_clean <- vax.2018_2021 %>%
  mutate(cov19 = ifelse(VAX_TYPE == "COVID19", 1, 0)) %>%
  mutate(varzos = ifelse(VAX_TYPE == "VARZOS", 1, 0)) %>%
  mutate(flu4 = ifelse(VAX_TYPE == "FLU4", 1, 0)) %>%
  mutate(ppv = ifelse(VAX_TYPE == "PPV", 1, 0)) %>%
  mutate(hpv9 = ifelse(VAX_TYPE == "HPV9", 1, 0)) %>%
  mutate(sp = ifelse(VAX_TYPE == "SMALL", 1, 0)) %>%
  mutate(yf = ifelse(VAX_TYPE == "YF", 1, 0)) %>%
  group_by(VAERS_ID) %>%
  summarise(
    cov19 = max(cov19),
    varzos = max(varzos),
    flu4 = max(flu4),
    ppv = max(ppv),
    hpv9 = max(hpv9),
    sp = max(sp),
    yf = max(yf),
    .groups = "drop"
  )
# ______________________________________________________________________________


# ______________________________________________________________________________
# Joining VAERS, symptom, and vaccine data
# some observations will be lost due to IDs present in VAERS data not being present in SYMPTOMS data and vice versa
vaers_dat <- vaers_clean %>%
  inner_join(sym_clean, by="VAERS_ID") %>%
  inner_join(vax_clean, by="VAERS_ID") %>%
  select(c(1, 37, 38, 4, 7, 21, 39:45, 46:52))
# ______________________________________________________________________________


# ______________________________________________________________________________
# Saving data for future use

save(vaers_dat, file = "vaers.2018-2021.RData")
# ______________________________________________________________________________
