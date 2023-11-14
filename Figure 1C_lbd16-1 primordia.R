library(tidyverse)
library(emmeans)

setwd("/Users/Dorota/Dropbox/PhD/LBD manuscript/")

E2 <-read_csv("20221101-Primordia.csv")
head(E2)

#Generalized linear model for data with binominal distribution
model1 <- glm(cbind(`non-emerged`,emerged) ~ genotype*condition , family = binomial(link = "logit"), 
                  data = E2)
summary(model1)
plot(model1)

#lsmeans with manual contrasts
emm1 <- emmeans(model1, ~ genotype*condition)
summary(emm1)

#contrast for comparison of responses to 75/125mM NaCl for each genotype
cont1 <- list(col75resp=c(0,0,1,0,-1,0),
             col125resp=c(1,0,0,0,-1,0),
             lbd75resp=c(0,0,0,1,0,-1),
             lbd125resp=c(0,1,0,0,0,-1))
contrast1  <- contrast(emm1, method = cont1, type ="response")
contrast1

#contrast for comparison of differences btw genotypes at 0,75 and 125 mM NaCl

cont2 <- list(lbdvsCol0=c(0,0,0,0,-1,1),
              lbdvsCol75=c(0,0,-1,1,0,0),
              lbdvsCol125=c(-1,1,0,0,0,0))
              
contrast2  <- contrast(emm1, method = cont2, type ="response")
contrast2


sink("Stats-LBD-primordia.txt")
print(summary(model1))
print(contrast1)
print(contrast2)
sink()
