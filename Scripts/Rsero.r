install.packages("devtools")
library(devtools)
library(tidyverse)
devtools::install_github("nathoze/Rsero") 
library(Rsero)
#example of a serological survey 1
 p.infection=0.2
 N.samples=500
 age = round(runif(n=N.samples, min =1, max = 70))
 seropositive = runif(n=N.samples)<p.infection 
simulated.survey  = SeroData(age_at_sampling = age, Y = seropositive) 


sex= c(rep('males',250), rep('females', 250))
simulated.survey  = SeroData(age_at_sampling =  age,
                             Y = seropositive,
                             sex = sex,
                             location  = 'Paris',
                             sampling_year = 2015) 
seroprevalence(simulated.survey)
#> [1] "Mean: 0.23    2.5%: 0.19    97.5%: 0.27"
seroprevalence.plot(simulated.survey,YLIM=0.3)


sex= c(rep('males',250), rep('females', 250))
simulated.survey  = SeroData(age_at_sampling =  age, Y = seropositive, sex = sex, location  = 'Paris', sampling_year = 2015, category = sex) 


#example of a serological survey 2
data('one_peak_simulation')


# seroprevalence estimation
seroprevalence(one_peak_simulation)
#> [1] "Mean: 0.11    2.5%: 0.09    97.5%: 0.15"

#plot of seroprevalence

seroprevalence.plot(one_peak_simulation)
#> [[1]]

ConstantModel = FOImodel(type = 'constant')

FOIfit.constant = fit(data = one_peak_simulation,  model = ConstantModel, chains=1)

seroprevalence.fit(FOIfit.constant, YLIM=0.5)
#> [[1]]


#model of one outbreak
OutbreakModel = FOImodel( type='outbreak', K=1)

FOIfit.outbreak = fit( data = one_peak_simulation,  model = OutbreakModel, chains=1)
seroprevalence.fit(FOIfit.outbreak)          

#compare the two models
DIC.constant = compute_information_criteria(FOIfit.constant)
DIC.outbreak = compute_information_criteria(FOIfit.outbreak)


#vietnam chik data 
#load the data and organise them in a SeroData Format

data('chikv_vietnam')
chikv_vietnam$location =  as.character(chikv_vietnam$location)
chik.sero = SeroData(age_at_sampling = chikv_vietnam$age,
                     Y = chikv_vietnam$Y,
                     category = chikv_vietnam$location,
                     reference.category = 'Hue',
                     sampling_year = 2017)



data('chikv_vietnam')
chikv_vietnam$location =  as.character(chikv_vietnam$location)
chik.sero.aggregated = SeroData(age_at_sampling = chikv_vietnam$age,
                                Y = chikv_vietnam$Y,
                                sampling_year = 2017)

seroprevalence(chik.sero.aggregated) # Value of the seroprevalence 

seroprevalence.plot(chik.sero.aggregated) # plots of the seroprevalence vs age




#real data analysis
library(readxl)
library(here)
#import data
  
chikdata <- read_excel(here("Data", "chikungunya_data_Uganda.xlsx"))



chikdata<- chikdata %>% select(Age_Yrs, IgM_CHIK, Year, Sex)

#omit na.omit in chikdata
chikdata <- na.omit(chikdata)



#recode Nengative to Negative
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "Nengative"] <- "Negative"

# Recode "NA" string to actual NA
chikdata$IgM_CHIK[chikdata$IgM_CHIK == "NA"] <- NA

# Drop all real NA values
chikdata <- chikdata[!is.na(chikdata$IgM_CHIK), ]



#recode Nengative to Negative
chikdata$Sex[chikdata$Sex == "FEAMALE"] <- "FEMALE"

# Recode "NA" string to actual NA
chikdata$Sex[chikdata$Sex == "NA"] <- NA

# Drop all real NA values
chikdata <- chikdata[!is.na(chikdata$Sex), ]

#recode the IgM_CHIK variable, Positive to TRUE , Negative to FALSE
chikdata$IgM_CHIK <- ifelse(chikdata$IgM_CHIK == "Positive", TRUE, FALSE)


chikdata$Sex=  as.character(chikdata$Sex)
chikdata$Age_Yrs=  as.integer(chikdata$Age_Yrs)
#convert chikdata to SeroData format

chik.sero.aggregated = SeroData(age_at_sampling = chikdata$Age_Yrs,
                                Y = chikdata$IgM_CHIK,
                                category = chikdata$Sex,
                                reference.category = "MALE",
                                sampling_year = 2020)
#seroprevalence estimation
seroprevalence(chik.sero.aggregated) # Value of the seroprevalence
seroprevalence.plot(chik.sero.aggregated) # plots of the seroprevalence vs age
#fit the model constant model
ConstantModel = FOImodel(type = 'constant')
FOIfit.constant = fit(data = chik.sero.aggregated,  model = ConstantModel, chains=1)
seroprevalence.fit(FOIfit.constant, YLIM=0.5)

#fit the model outbreak model
OutbreakModel = FOImodel( type='outbreak', K=1)
FOIfit.outbreak = fit( data = chik.sero.aggregated,  model = OutbreakModel, chains=1)
seroprevalence.fit(FOIfit.outbreak)
#compare the two models
  
DIC.constant = compute_information_criteria(FOIfit.constant)
DIC.outbreak = compute_information_criteria(FOIfit.outbreak)
# DIC values
DIC.constant
DIC.outbreak

#Model of the force of infection

#Some models have additional parameters, for instance the outbreak models require that we specify the number of peaks K.


model = list()
model[[1]] = FOImodel(type='constant')
model[[2]] = FOImodel(type='piecewise', K=2)
model[[3]] = FOImodel(type='piecewise', K=3)
model[[4]] = FOImodel(type='outbreak', K=1)
model[[5]] = FOImodel(type='outbreak', K=2)
model[[6]] = FOImodel(type='constantoutbreak', K=1)
model[[7]] = FOImodel(type='constantoutbreak', K=2)
model[[8]] = FOImodel(type='independent')
model[[9]] = FOImodel(type='constant', sp=0.95, se=0.99)
model[[10]] = FOImodel(type='outbreak', K=1, sp=0.95, se=0.99)
model[[11]] = FOImodel(type='constant', seroreversion = 1)



#Fitting a model to the serological data

#Our objective is to infer the parameters for the differents models, assess whether the models are adequate at reproducing the data, and select the best model(s)
Fit_9 = fit(data = chik.sero.aggregated,
            model = model[[9]],
            chain=4,
            cores=4, 
            iter=5000)

Fit_10 = fit(data = chik.sero.aggregated,
             model = model[[10]],
             chain=4,
             cores=4, 
             iter=5000)


#Analyzing the results
#The Rsero packages includes a certain number of functions to analyze the MCMC chains and the posterior distributions of the parameters. We first compare the different models and evaluate the goodness of fit
#The compute_information_criteria function implemented in Rsero evaluates the AIC, DIC, WAIC, and the PSIS-LOO
#The AIC, which requires an estimation of the maximum likelihood, is obtained for the ensemble of parameters in the chain that maximizes the likelihood

C_9 = compute_information_criteria(Fit_9)
C_10 = compute_information_criteria(Fit_10)

print(C_9)
print(C_10)


p = seroprevalence.fit(Fit_9)
print(p)
p = seroprevalence.fit(Fit_10)
print(p)



plot_posterior(Fit_10)


parameters_credible_intervals(Fit_10)


plot(Fit_10, YLIM =0.3) 



Fit_10$data
Fit_10$model
Fit_10$fit


traceplot_Rsero(Fit_10)




