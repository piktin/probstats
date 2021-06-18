##############################################
##    END OF TERM ASSIGNMENT            ######
#   Probability and Statistics          ######
##                                      ######
##############################################


#install necessary packages
install.packages("openxlsx")
install.packages('caret')
install.packages("data.table")
install.packages("ggpubr")

#calling the needed libraries:
library('openxlsx')
library('data.table')
library("ggpubr")
library("dplyr")
library(tidyverse)
library(caret)


#set file path to the excel file biomarkers
file_path <- "/Probability and Statistics/Assesment/biomarkers.xlsx"

##Read the file into dataframe
bio <- read.xlsx(file_path)

typeof(bio)
#View dataframe
View(bio)

#Find more info about dataframe, only biomarkers has characters.
str(bio)

summary(bio)

#check for invalid values
table(is.na(bio)) #there are 9 NANs in the set

# set filepath for covariates
file_path1 <-"/Probability and Statistics/Assesment/covariates.xlsx"


###############################################
##Question 1: Statistical hypothesis testing ##
###############################################


bio1 <- as.data.table(read.xlsx(file_path)) #open xlsx file with data.table

#Split the biomarkers into 2 columns with Px Id and Time point
bio1[, c("PatientID",
         "Timepoint") :=
       tstrsplit(Biomarker, "-",
                 fixed = TRUE)]

View(bio1)# Check that new columns are made

#Convert the Patient ID to numeric to sort the table
bio1[, PatientID :=
       as.numeric(PatientID)]

#Sort the table based on Patient Id (numeric)
(bio_sort <- bio1[order(PatientID),])

#Here you see that row 348 has only NANs values in all columns

#Check how many NANs values in bio_sort df.
table(is.na(bio_sort))

#Create new df bio1_ that removes all rows that have NANs value (row 348)
bio1_ <- na.omit(bio_sort)


#Creates new df bio2 that convert long data into wide data: move the timepoint of different px + timepoint to different columns 
bio2 <- dcast(bio1_,
              PatientID ~ Timepoint,
              value.var = c("IL-8", "VEGF-A","OPG","TGF-beta-1", "IL-6", "CXCL9", "CXCL1","IL-18","CSF-1"))
View(bio2)


#### Read the covariates into dataframe cova_df and find more info about it

cova_df <- as.data.table(read.xlsx(file_path1))

# Find out about NANs values in cova_df: there are 2 NANs value in the set
table(is.na(cova_df))


#Left merge cova_df with bio2 (wide) df based on PatientID
bio_cova <- merge(bio2, cova_df, 
      all.x = TRUE,
      by = "PatientID")
View(bio_cova)


############################
##    Males only          ##
############################

#Create df of males patients only from covariates file
cova_male <- cova_df[`Sex.(1=male,.2=female)`
                    == 1,]

#Semijoin the males bio_cova to the bio2 df, creating new df showing only those Px who are males
setkey(bio2, PatientID)
male_ <- bio2[bio_cova[cova_male,
                   which = TRUE]]
View(male_)

#Find means of all columns in males patients
mal <- colMeans(male_, na.rm = TRUE)
mal

##################################
###     FeMales only         ####
##################################

#Create df of females patients only from covariates file

cova_female <- cova_df[`Sex.(1=male,.2=female)`
                     == 2,]


##Semijoin the female bio_cova to the bio2 df, showing only those Px who are females

setkey(bio2, PatientID)

female_ <- bio2[bio_cova[cova_female,
                       which = TRUE]]
View(female_)

#Find means of all colums in females patients
fem <- colMeans(female_, na.rm = TRUE)
fem

##create df of males px of 9 biomarkers at inclusion only:
vars <- grepl("0", names(male_))
c_male <- male_[, ..vars]
View(c_male)

#find means of all columns in males
male_mean <- colMeans(c_male, na.rm = TRUE)
male_mean

#Create df of females px of 9 biomarkers at inclusion
vars_f <- grepl("0", names(female_))
c_female <- female_[, ..vars_f]

#find means of all columns in females
female_mean <- colMeans(c_female, na.rm = TRUE)
female_mean

##################
###   IL-8     ###
##################

#Create df of males px of IL-8 at inclusion (0 months)

vars_il8 <-  grepl("IL-8_0", names(male_))
c_male_il8 <- male_[,..vars_il]


#Create df of females px of IL-8 at inclusion (0 months)

vars_il8_f <- grepl("IL-8_0", names(female_))
c_female_il8 <- female_[, ..vars_il8_f] 


# Check if it is normally distributed
shapiro.test(c_female_il8$'IL-8_0weeks') # it is not, p is 0.001489
shapiro.test(c_male_il8$'IL-8_0weeks') # it is n.distributed, p = 0.5545


#Since it is NOT normally distributed, use fligner.test to test if variance is homogeneneos
fligner.test(c_male_il8$'IL-8_0weeks', c_female_il8$'IL-8_0weeks')
# p = 0.3351 > 0.05 so varequal = true

#Check if the variance is equal between the 2 samples using F.test, assuming that data of biomarkers are normally-dis
var.test(c_male_il8$'IL-8_0weeks', c_female_il8$'IL-8_0weeks', alternative = "two.sided")
#p-value = 0.5645 so not reject H0 that var = equal

#Find t distribution information about
#the null hypothesis:no different in IL-8 level at inclusion bwt males and females, removed all NANs value
#Assuming they are normally distributed

t.test(c_male_il8, c_female_il8,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE) #p=0.3294


#Since the data is not normally distributed, use unpaired 2 samples Wilcox test would be more correct:
wilcox.test(c_male_il8$`IL-8_0weeks`, c_female_il8$`IL-8_0weeks`,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p=0.4372



##################################################
# Vascular endothelial growth factor A (VEGF-A) ##
##################################################

#Create df of males px of VEGF-A at inclusion (0 months)

vars_vegf <-  grepl("VEGF-A_0", names(male_))
c_male_vegf <- male_[,..vars_vegf]


#Create df of females px of VEGF-A at inclusion (0 months)

vars_vegf_f <- grepl("VEGF-A_0", names(female_))
c_female_vegf <- female_[, ..vars_vegf_f] 


# Check if is is normally distributed
shapiro.test(c_female_vegf$'VEGF-A_0weeks') # it is n.d, p = 0.6825
shapiro.test(c_male_vegf$'VEGF-A_0weeks') # it is not n.d, p = 0.038

#Since there is no normally distributed, use fligner.test to test if var homogeneneos
fligner.test(c_male_vegf$'VEGF-A_0weeks', c_female_vegf$'VEGF-A_0weeks')
# p = 0.2998, varequal = true


#Find t distribution information about
#the null hypothesis:no different in VEGF-A level at inclusion bwt males and females, removed all NANs value

t.test(c_male_vegf,c_female_vegf,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
# p-value = 0.04139

#Since males populations is not n.d, use wilcox test to be more correct:

wilcox.test(c_male_vegf$'VEGF-A_0weeks', c_female_vegf$'VEGF-A_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p = 0.04146


#######################
#####     OPG     ####
#######################


#Create df of males px of OPG at inclusion (0 months)

vars_opg <-  grepl("OPG_0", names(male_))
c_male_opg <- male_[,..vars_opg]


#Create df of females px of OPG at inclusion (0 months)

vars_opg_f <- grepl("OPG_0", names(female_))
c_female_opg <- female_[, ..vars_opg_f] 


#Check if the variance is equal between the 2 samples using F.test, assuming that data of biomarkers are normally-dis
var.test(c_male_opg$'OPG_0weeks', c_female_opg$'OPG_0weeks', alternative = "two.sided")
#p-value = 0.4678 so not reject H0 that var1 = var2


# Check if is is normally distributed
shapiro.test(c_female_opg$'OPG_0weeks') # it is n.d, p = 0.305
shapiro.test(c_male_opg$'OPG_0weeks') # it is not n.d, p = 0.000186

#Since there is no normally distributed, use fligner.test to test if variance is homogeneneos
fligner.test(c_male_opg$'OPG_0weeks', c_female_opg$'OPG_0weeks')
# p is 0.3691, varequal is true



#Find t distribution information about
#the null hypothesis:no different in OPG level at inclusion bwt males and females, removed all NANs value

t.test(c_male_opg,c_female_opg,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)


#Since males populations is not n.d, use wilcox test to be more correct:
wilcox.test(c_male_opg$'OPG_0weeks', c_female_opg$'OPG_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p = 0.06536

################
#TGF-beta-1  ###
################

#Create df of males px of TGF-beta-1 at inclusion (0 months)

vars_tgf <-  grepl("TGF-beta-1_0", names(male_))
c_male_tgf <- male_[,..vars_tgf]


#Create df of females px of TGF-beta-1 at inclusion (0 months)

vars_tgf_f <- grepl("TGF-beta-1_0", names(female_))
c_female_tgf <- female_[, ..vars_tgf_f] 


#Check if the variance is equal between the 2 samples using F.test, assuming that data of biomarkers are normally-dis
var.test(c_male_tgf$'TGF-beta-1_0weeks', c_female_tgf$'TGF-beta-1_0weeks', alternative = "two.sided")
#p-value = 0.9278 so not reject H0 that var1 = var2

# Check if is is normally distributed
shapiro.test(c_female_tgf$'TGF-beta-1_0weeks') # it is not n.d, p = 0.00499
shapiro.test(c_male_tgf$'TGF-beta-1_0weeks') # it is not n.d p= 0.001692

#Since there is no normally distributed, use fligner.test to test if var homogeneneos
fligner.test(c_male_tgf$'TGF-beta-1_0weeks', c_female_tgf$'TGF-beta-1_0weeks')
#p is 0.3368, varequal is true 

#Find t distribution information about
#the null hypothesis:no different in TGF-beta-1 level at inclusion bwt males and females, removed all NANs value

t.test(c_male_tgf,c_female_tgf,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)

#Since both populations are not n.d, use wilcox test to be more correct:
wilcox.test(c_male_tgf$'TGF-beta-1_0weeks', c_female_tgf$'TGF-beta-1_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p = 0.05163

##################
##     IL-6   ####
##################

#Create df of males px of IL-6 at inclusion (0 months)

vars_il6 <-  grepl("IL-6_0", names(male_))
c_male_il6 <- male_[,..vars_il6]


#Create df of females px of IL-6 at inclusion (0 months)

vars_il6_f <- grepl("IL-6_0", names(female_))
c_female_il6 <- female_[, ..vars_il6_f] 

#Check if the variance is equal between the 2 samples using F.test, assuming that data of biomarkers are normally-dis
var.test(c_male_il6$'IL-6_0weeks', c_female_il6$'IL-6_0weeks', alternative = "two.sided")
#p-value = 0.03585 less than 0.05 alpha level so Reject H0 that var1 = var2, so var equal= false

# test for normal distributed. which both data is so use of var.test was appropriate
shapiro.test(c_male_il6$'IL-6_0weeks') # it it not n.d, p = 0.0004733
shapiro.test(c_female_il6$'IL-6_0weeks') # is is not n.d, p = 2.00xe^-5

#Since there is no normally distributed, use fligner.test to test if var homogeneneos
fligner.test(c_male_il6$'IL-6_0weeks', c_female_il6$'IL-6_0weeks')
# p = 0.4006, so varequal = true



#Find t distribution information about
#the null hypothesis:no different in IL-6 level at inclusion bwt males and females, removed all NANs value

t.test(c_male_il6, c_female_il6,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)

#Since both populations are not n.d, use wilcox test to be more correct:
wilcox.test(c_male_il6$'IL-6_0weeks', c_female_il6$'IL-6_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p = 0.3597


####################
###    CXCL-9   ####
####################

#Create df of males px of CXCL9 at inclusion (0 months)

vars_cxcl9 <-  grepl("CXCL9_0", names(male_))
c_male_cxcl9 <- male_[,..vars_cxcl]


#Create df of females px of CXCL9 at inclusion (0 months)

vars_cxcl9_f <- grepl("CXCL9_0", names(female_))
c_female_cxcl9 <- female_[, ..vars_cxcl9_f] 


#Check if the variance is equal between the 2 samples using F.test, assuming that data of biomarkers are normally-dis
var.test(c_male_cxcl9$'CXCL9_0weeks',c_female_cxcl9$'CXCL9_0weeks', alternative = "two.sided")
#p-value = 0.5463 so NOT Reject H0 that var1 = var2, so var equal= true

# test for normal distributed. which both data is so use of var.test was appropriate
shapiro.test(c_male_cxcl9$'CXCL9_0weeks') # it it not n.d, p = 1.613e-05
shapiro.test(c_female_cxcl9$'CXCL9_0weeks') # it is not n.d, p = 2.108e-08

#Since both populations are not normally distributed, best to use fligner.test to check for if var are equal
fligner.test(c_male_cxcl9$'CXCL9_0weeks', c_female_cxcl9$'CXCL9_0week')
# p is 0.15 so varequal is true

#Find t distribution information about
#the null hypothesis:no different in CXCL9 level at inclusion bwt males and females, removed all NANs value

t.test(c_male_cxcl9, c_female_cxcl9,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)

#Since both populations are not n.d, use wilcox test to be more correct:
wilcox.test(c_male_cxcl9$'CXCL9_0weeks',c_female_cxcl9$'CXCL9_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p = 0.8615


###############
##  CXCL1   ###
###############


#Create df of males px of CXCL1 at inclusion (0 months)

vars_cxcl1 <-  grepl("CXCL1_0", names(male_))
c_male_cxcl1 <- male_[,..vars_cxcl1]


#Create df of females px of CXCL1 at inclusion (0 months)

vars_cxcl1_f <- grepl("CXCL1_0", names(female_))
c_female_cxcl1 <- female_[, ..vars_cxcl1_f] 


#Check if the variance is equal between the 2 samples using F.test, assuming that data of biomarkers are normally-dis
var.test(c_male_cxcl1$'CXCL1_0weeks',c_female_cxcl1$'CXCL1_0weeks', alternative = "two.sided")
#p-value = 0.9943 so NOT Reject H0 that var1 = var2, so var equal= true


# test for normal distributed. which both data is so use of var.test was appropriate
shapiro.test(c_male_cxcl1$'CXCL1_0weeks') # it it not n.d, p = 7.515e-05
shapiro.test(c_female_cxcl1$'CXCL1_0weeks') # it is not n.d, p = 0.002602

#Since both populations are not normally distributed, best to use fligner.test to check for if var are equal
fligner.test(c_male_cxcl1$'CXCL1_0weeks', c_female_cxcl1$'CXCL1_0weeks')
#p is 0.41 so varequal is true

#Find t distribution information about
#the null hypothesis:no different in CXCL1 level at inclusion bwt males and females, removed all NANs value

t.test(c_male_cxcl1, c_female_cxcl1,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)

#Since both populations are not n.d, use wilcox test to be more correct:
wilcox.test(c_male_cxcl1$'CXCL1_0weeks', c_female_cxcl1$'CXCL1_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p = 0.003242


##################
## IL-18      ####
##################

#Create df of males px of IL-18 at inclusion (0 months)

vars_il18 <-  grepl("IL-18_0", names(male_))
c_male_il18 <- male_[,..vars_il18]


#Create df of females px of IL-18 at inclusion (0 months)

vars_il18_f <- grepl("IL-18_0", names(female_))
c_female_il18 <- female_[, ..vars_il18_f] 


# test for normal distributed. which both data is so use of var.test was appropriate
shapiro.test(c_male_il18$'IL-18_0weeks') # it it n.d, p = 0.612
shapiro.test(c_female_il18$'IL-18_0weeks') # it is n.d, p = 0.9625

#Since data is normally distributed, using F-test is suitable:
var.test(c_male_il18$'IL-18_0weeks', c_female_il18$'IL-18_0weeks', alternative = "two.sided")
#p-value = 0.7083 bigger than 0.05 alpha level so NOT Reject H0 that var1 = var2, so var equal= true


#Use t-test as data is normally distributed and variance is equal
#Find t distribution information about
#the null hypothesis: no different in IL-18 level at inclusion bwt males and females, removed all NANs value

t.test(c_male_il18, c_female_il18,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)

#wilcox.test:
wilcox.test(c_male_il18$'IL-18_0weeks', c_female_il18$'IL-18_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p-value = 0.2853


##################
###    CSF-1  ####
##################

#Create df of males px of CSF-1 at inclusion (0 months)

vars_csf <-  grepl("CSF-1_0", names(male_))
c_male_csf <- male_[,..vars_csf]


#Create df of females px of CSF-1 at inclusion (0 months)

vars_csf_f <- grepl("CSF-1_0", names(female_))
c_female_csf <- female_[, ..vars_csf_f] 

#Check if the variance is equal between the 2 samples using F.test, assuming that data of biomarkers are normally-dis
var.test(c_male_csf$'CSF-1_0weeks', c_female_csf$'CSF-1_0weeks', alternative = "two.sided")

#p-value = 0.005419 < 0.05 alpha level so Reject H0 that var1 = var2, so var equal= False
shapiro.test(c_male_csf$'CSF-1_0weeks')# it is normally distributed, p = 0.9218
shapiro.test(c_female_csf$'CSF-1_0weeks') # it is not n.d, p = 0.04536

#Since one of populations is not n.d, best to use Fligner test to check for variance equality
fligner.test(c_male_csf$'CSF-1_0weeks', c_female_csf$'CSF-1_0weeks')
#p = 0.4084, thus varqequal = true


#Find t distribution information about
#the null hypothesis:no different in CSF-1 level at inclusion bwt males and females, removed all NANs value

t.test(c_male_csf, c_female_csf,var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)

#Since one of populations is not n.d, use wilcox test to be more correct:
wilcox.test(c_male_csf$'CSF-1_0weeks', c_female_csf$'CSF-1_0weeks',var.equal=TRUE, conf.level= 0.95, na.rm= TRUE)
#p = 0.01114

##################################
#Part 1c: calculate type 1 error: #
##################################

type1_error <- 1- (1-0.05)^9
type1_error

# With Bonferroni corection: alpha/n = 0.05/9
alfa_bonf <- 0.05/9
alfa_bonf
type1_error_bonf <- 1 - (1-(0.05/9))^9
type1_error_bonf

#make dataframe with biomarkers and raw-p-values found in t-tests in previous section        
biomarkers <- c('IL8','VEGF-A', 'OPG','TGF-beta-1','IL-6','CXCL9','CXCL1','IL-18','CSF-1')
p_raw_t <- c(0.3294,0.04139,0.1326,0.0457,0.2515,0.9849,0.006167, 0.2496 ,0.005807 )

biomarkers_p_raw <- data.frame(biomarkers, p_raw_t)

#find adjusted p-values of t-tests:
bonfe_t <- p.adjust(biomarkers_p_raw$p_raw_t,
         method = "bonferroni")
bonfe_t


#make dataframe with biomarkers and raw-p-values found in wilcox-tests in previous section 
p_raw_wil <- c(0.4372 , 0.04146 , 0.06536, 0.05163 , 0.3597 , 0.8615 , 0.003242 , 0.2853, 0.01114  )

biomarker_p_raw_wil <- data.frame(biomarkers, p_raw_wil)

#find adjusted p-values of wilcox tests
bonfe_wil <- p.adjust(biomarker_p_raw_wil$p_raw_wil,
                    method = "bonferroni")
bonfe_wil


####################################
###   PART 2: Regression Model  ####
####################################


##########
#Part a ##
##########

#remove NANs in the bio_cova df
bio_cova_om <- na.omit(bio_cova)


# Create df that only includes the variables of interest, biomarkers at inclusion and covariates
BIO_COVA <- select(bio_cova_om,'Vas-12months', 'IL-8_0weeks' ,'VEGF-A_0weeks' , 'OPG_0weeks', 'TGF-beta-1_0weeks' , 'IL-6_0weeks' ,
                            'CXCL9_0weeks' , 'CXCL1_0weeks' , 'IL-18_0weeks' , 'CSF-1_0weeks' , 'Age', 'Sex.(1=male,.2=female)' ,
                            'Smoker.(1=yes,.2=no)', 'VAS-at-inclusion')
View(BIO_COVA)

#set random variables:
set.seed(123)

#rename the columns of the BIO-COVA set to shorter names, easier to use later for producing training model:
BIO_COVA <- BIO_COVA %>% 
  rename (VAS12 = 'Vas-12months', IL8 = 'IL-8_0weeks',VEGFA = 'VEGF-A_0weeks', OPG = 'OPG_0weeks',
         TGFB1 = 'TGF-beta-1_0weeks', IL6 = 'IL-6_0weeks', CXCL9 = 'CXCL9_0weeks' ,CXCL1 = 'CXCL1_0weeks',
         IL18 = 'IL-18_0weeks' ,CSF1 =  'CSF-1_0weeks' , Age = 'Age', Sex = 'Sex.(1=male,.2=female)' ,
         Smoker = 'Smoker.(1=yes,.2=no)', VASinclusion = 'VAS-at-inclusion')

################################
##Process to split the data:  ##
################################

#set n1 as length of the BIO_COVA df
n1 <- nrow(BIO_COVA)

#set number of testing set (20% of total data)
ntest <- round(0.2*n1)

#set number of training  set (80% of total data)
ntrain <- n1 - ntest

# Split the data into two sets: training and test sets
train_rows <- sample(1:n1, ntrain)
BIO_COVA_train <- BIO_COVA[train_rows,]
BIO_COVA_test <- BIO_COVA[-train_rows,]


# Build the model using the training set (80% of the patients)
model <- lm(VAS12 ~ IL8 + VEGFA + OPG + TGFB1 + IL6 + CXCL9 + CXCL1 + IL18 + CSF1 + Age + Sex + Smoker + VASinclusion,
            data = BIO_COVA_train, na.action=na.omit)

# Summarize the model
summary(model)

#view the training set
View(BIO_COVA_train)


#Predict using training set:
predict(model)

#Create a explo_var df that comprises of biomarkers at inclusion and other covariates, for ease of use later
explo_var <- BIO_COVA_train$IL8 + BIO_COVA_train$VEGFA +BIO_COVA_train$OPG+ BIO_COVA_train$TGFB1 + BIO_COVA_train$IL6+
  BIO_COVA_train$CXCL9 + BIO_COVA_train$CXCL1 + BIO_COVA_train$IL18 + BIO_COVA_train$CSF1 + BIO_COVA_train$Age+
  BIO_COVA_train$Sex +BIO_COVA_train$Smoker + BIO_COVA_train$VASinclusion


#Measurements of model prediction error using training set:

rmse_train <- sqrt(mean((predict(model, BIO_COVA_train) - BIO_COVA_train$VAS12)^2))
rmse_train
#2.345238

mae_train <- mean(abs(predict(model, BIO_COVA_train) - BIO_COVA_train$VAS12))
mae_train
#1.938147


###############################
##Plot using training set ####
###############################

#Simple plot using training set, explanatory variables vs response:
plot(explo_var, BIO_COVA_train$VAS12,
     xlab="Various biomarkers at inclusion and other covariates",
     ylab="VAS at 12 months",main="VAS in 12 months and different covariates",
     col=2)
abline(lm(BIO_COVA_train$VAS12 ~ explo_var), col = 'blue')

#Find correlation of the predicted(fitted) and the actual values of VAS at 12 months in training set:
cor(BIO_COVA_train$VAS12,fitted)
# cor = 0.6379601


#Try plot with logs of covariates/explanatory variables
explo_log <- log(explo_var)

plot(explo_log, BIO_COVA_train$VAS12,
     xlab="Various biomarkers at inclusion and other covariates",
     ylab="VAS at 12 months",main="VAS in 12 months and different covariates",
     col=2)
abline(lm(BIO_COVA_train$VAS12 ~ explo_log), col = 'blue')

#####Set sub-data sets for easier visualisation later:
fitted <- predict(model)
resid <- residuals(model)

#save to csv file the fitted values from training set:
write.csv(fitted, "C:/Users/Wind/Documents/Edinburgh/Probability and Statistics/Assesment/Bio_Cova_Training_Predict.csv")

#histogram of residuals from training model
hist(resid, xlab= 'Residuals', ylab='Frequency', main = ' Histogram of redisuals of training model')


#Using ggplot for better visuals
#create df2 that is made up from predicted values of model and residuals
df2 <- data.frame(
  resid_train = residuals(model),
  pred_train = predict(model))

#Use ggplot: predicted values vs abs residuals
plot_train <- ggplot(df2, aes(pred_train, abs(resid_train))) +
  geom_point() +
  geom_smooth()

plot_train +  ggtitle("Predicted Values from training set vs absolute values of residuals") +
  xlab("Predicted Values using train set") + ylab("Absolute values of residuals")
  
  
#create predicted df that made of predicted values and covariates
predicted_df <- data.frame(vas12_pred = predict(model, BIO_COVA_train), covariates=BIO_COVA_train$IL8 + BIO_COVA_train$VEGFA +BIO_COVA_train$OPG+ BIO_COVA_train$TGFB1 + BIO_COVA_train$IL6+
                             BIO_COVA_train$CXCL9 + BIO_COVA_train$CXCL1 + BIO_COVA_train$IL18 + BIO_COVA_train$CSF1 + BIO_COVA_train$Age+
                             BIO_COVA_train$Sex +BIO_COVA_train$Smoker + BIO_COVA_train$VASinclusion)


# this is the predicted line of multiple linear regression
plot_mul <- ggplot(data = predicted_df, aes(x = covariates, y = vas12_pred)) + 
  geom_point(color='red') +
  geom_line(color='blue',data = predicted_df, aes(x=covariates, y=vas12_pred))

plot_mul +  ggtitle("Plot of Covariates vs VAS at 12 months using training set") +
  xlab("Biomarkers level at inclusion and other covariates") + ylab("VAS at 12 months")


#set train set with prediction interval
train_pred_int <- predict(model, BIO_COVA_train, interval = "prediction")

#bind it together in 1 df with prediction intervals
mydata_train <- cbind(BIO_COVA_train, train_pred_int)
View(mydata_train)

#set new set of covariates from training set_
covariates=BIO_COVA_train$IL8 + BIO_COVA_train$VEGFA +BIO_COVA_train$OPG+ BIO_COVA_train$TGFB1 + BIO_COVA_train$IL6+
  BIO_COVA_train$CXCL9 + BIO_COVA_train$CXCL1 + BIO_COVA_train$IL18 + BIO_COVA_train$CSF1 + BIO_COVA_train$Age+
  BIO_COVA_train$Sex +BIO_COVA_train$Smoker + BIO_COVA_train$VASinclusion

#create plot using training set: covariates vs VAS12
ptrain <- ggplot(mydata_train, aes(covariates, BIO_COVA_train$VAS12)) +
  geom_point() +
  stat_smooth(method = lm)

# Add prediction intervals:
ptrain + geom_line(aes(y = lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y = upr), color = "red", linetype = "dashed") +
  ggtitle("Plot of Covariates vs VAS at 12 months using training set with Prediction Interval in red") +
  xlab("Biomarkers level at inclusion and other covariates") + ylab("VAS at 12 months")

#Make plot of predicted values vs actual values using training set
ptrain1 <- ggplot(mydata_train, aes(x = fit, y= VAS12)) +
  geom_point(color='red') +
  stat_smooth(method = 'lm', color = 'blue') +
  ggtitle("Plot of Predicted values vs Actual Values using training set") +
  xlab("Predicted Values for VAS at 12 months") + ylab("Actual Values")

ptrain1

#Another way to make plot. simple without prediction intervals:
plot_sing <- ggplot(data = predicted_df, aes(x = covariates, y = vas12_pred)) + 
  geom_point(color='red') +
  geom_smooth(method = "lm", se = FALSE, color = 'blue')

plot_sing +  ggtitle("Plot of Covariates vs VAS at 12 months using training set") +
  xlab("Biomarkers level at inclusion and other covariates") + ylab("VAS at 12 months")


#check if linearity can be assumed in our model. As it looks like horizontal line, residual vs fitted
#plot shows no pattern, we can assume linearity
#residual plot and others:
autoplot(model)



#################################################
#####           Use test set              #######
#################################################

#rmse, mae for test set using model from training data but with test set (out of sample evaluation)
rmse <- sqrt(mean((predict(model, BIO_COVA_test) - BIO_COVA_test$VAS12)^2))
rmse
#3.208264

mae <- mean(abs(predict(model, BIO_COVA_test) - BIO_COVA_test$VAS12))
mae
#2.697855


###########################
##    Make predictions  ###
###########################

#set random variables:
set.seed(123)

predict_test <- predict(model, newdata = BIO_COVA_test, type = 'response')
predict_test


##################################################
## Evaluating prediction error with test set  ####
##################################################

###Another way to measure the regression model performance: find R2, RMSE and MAE
predictions <- model %>% predict(BIO_COVA_test)

data.frame( R2 = R2(predictions, BIO_COVA_test$VAS12),
            RMSE = RMSE(predictions, BIO_COVA_test$VAS12),
            MAE = MAE(predictions, BIO_COVA_test$VAS12))
#R2=0.0976122, RMSE = 3.208264, MAE= 2.697855

#find prediction error rate from test set:
prediction_error_rate <- RMSE(predictions, BIO_COVA_test$VAS12)/mean(BIO_COVA_test$VAS12)
prediction_error_rate
#0.9083888

#check if correlation exists between predicted and actual values in test set
cor(BIO_COVA_test$VAS12, predict_test)
#cor = 0.312495, quite low.



########################################
####     PLOT  with test set       #####
########################################


#find residuals between observed and predicted in test set:
residuals_test <- BIO_COVA_test$VAS12 - predict_test 

#create subset of predict values using test set with prediction interval:
pred_test_int <- predict(model, newdata = BIO_COVA_test, interval = "prediction")

#bind all together in 1 dataframe to make plot
mydata <- cbind(BIO_COVA_test, pred_test_int, residuals_test)
View(mydata)

#subset of covariates in test set:
covariates_test=BIO_COVA_test$IL8 + BIO_COVA_test$VEGFA +BIO_COVA_test$OPG+ BIO_COVA_test$TGFB1 + BIO_COVA_test$IL6+
  BIO_COVA_test$CXCL9 + BIO_COVA_test$CXCL1 + BIO_COVA_test$IL18 + BIO_COVA_test$CSF1 + BIO_COVA_test$Age+
  BIO_COVA_test$Sex +BIO_COVA_test$Smoker + BIO_COVA_test$VASinclusion

#create plot of covariates vs VAS at 12 months from test set:
p <- ggplot(mydata, aes(covariates_test, BIO_COVA_test$VAS12)) +
  geom_point() +
  stat_smooth(method = lm)

# Add prediction intervals
p + geom_line(aes(y = lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y = upr), color = "red", linetype = "dashed") +
  ggtitle("Plot of Covariates vs VAS at 12 months using testing set with Prediction Interval in red") +
  xlab("Biomarkers level at inclusion and other covariates") + ylab("VAS at 12 months") + theme_classic()


#create plot of predicted values vs observed values:
plottest <- ggplot(mydata, aes(x = fit, y= VAS12)) +
  geom_point(color='red') +
  stat_smooth(method = 'lm', color = 'blue') +
  ggtitle("Plot of Predicted values vs Actual VAS at 12 months using testing set") +
  xlab("Predicted Values for VAS at 12 months") + ylab("Actual VAS at 12 months") + theme_classic()

plottest


#Create plot of predicted values vs residuals of test set:
plot_test_resid <- ggplot(mydata, aes(x = fit, y= residuals_test)) +
  geom_point(color='red') +
  stat_smooth(method = 'lm', color = 'blue') +
  ggtitle("Plot of Predicted values vs Residuals using testing set") +
  xlab("Predicted Values for VAS at 12 months") + ylab("Residuals") + theme_classic()

plot_test_resid

################################################################################
##                          End of Assignment                               ####
################################################################################
