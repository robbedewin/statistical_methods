library(readr)
titanic_missing <- read_csv("~/KULeuven/Bioinformatics_2024-2025/Statistical_Methods/Part I/Data and R-code/Lecture5_titanic.txt", 
                             col_types = cols(age = col_double()))
View(Lecture5_titanic)

titanic_missing<-subset(titanic_missing,select=c('survived','pclass','sex','age'))




# Load the library
library(mice)

# Tell mice to create 100 imputed datasets
imp <- mice(depression, m = 100)

# 'with' takes the imputed object and the model you want to run
fit <- with(data = imp, 
            exp = glm(survived ~ pclass + sex + age, family = binomial))

# 'pool' combines the results
final_estimates <- pool(fit)
summary(final_estimates)
summary(fit)


# How i would solve this
# Import Data

depression <- data...

# MAR because it is dependent on other observations

library(mice) 
depression_imputed <- mice(depression, m = 50)

# gaussian is the standard family, but to underline the difference with survival were we used binomial because data was 0 or 1?
fit <- with(data = depression_imputed, 
            exp = glm(BDI_final ~ BDI_baseline + income + treatment, family = gaussian))

final_estimates <- pool(fit)
summary(fit)

# For the Inverse Probability weighting

# Create a new variable that is 1 is complete and 0 otherwise
depression$r<-as.numeric(!is.na(depression$income))


## Fitting the logistic regression model to  calculate the probabilities of being complete

depression.ipw.glm<-glm(r ~ BDI_final + BDI_baseline + treatment, data=depression,family=binomial)
summary(depression.ipw.glm)

## Calculating the weights: Inverse Probabilities

depression.missing$w<-1/fitted(depression.ipw.glm)

# perform standard analysis but use w in weights argument


depression.results.ipw<- glm(BDI_final ~ BDI_baseline + income + treatment, 
                          data=depression, 
                          weights=depression$w, 
                          family=gaussian)
summary(depression.results.ipw)