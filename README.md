# Missing-Data
Simulated Missing Data is stored in the folder. We are using it for testing purposes only.
# This example estimates the  effect of participation in career academy programs on income
# This example used data from the Education Longitudinal Study of 2002 (ELS)
# More details about dataset can be found at http://nces.ed.gov/surveys/els2002/

#load dataset 
dataP<- read.csv(file="AERA_data1.csv")

#===============================================
#CODE FOR TREATMENT OF MISSING DATA BY MULTIPLE IMPUTATION
library(mice) #load package that performs multiple imputation by chained equations


#examine the number of missing cases
missing.indicator <- data.frame(is.na(dataP))
propMissing <- apply(missing.indicator,2,mean)

#create dummy missing value indicators
names(missing.indicator)[propMissing>0] <- paste(names(dataP)[propMissing>0],"NA",sep="")#convert dummy missing indicators from logical to numeric variables
for (var in 1:ncol(missing.indicator)) {
    missing.indicator[,var] <- as.numeric(missing.indicator[,var]) }


#merge covariate names with missing indicator names
dataP <- cbind(dataP,missing.indicator[,propMissing>0])


#show percentage missing
print(round(propMissing,3))


#create a treatment indicator
dataP$treat <- factor(dataP$trt)

#======================
#impute separately treated and untreated groups

#create storage
long.imputation <- c()

#loop through the treatment groups
for (group in 0:1) {


#creates a list of predictors of missing data with a mininum correlation of 0.1
#and at least 50% of useful data
predictor.selection <- quickpred(subset(dataP,treat==group), mincor=0.1, minpuc=0.5,method='pearson',
                                 exclude = "ID_group")


#impute variables by from least missing to most missing
#Using multiple imputation by chained equations
#with predictive mean matching as the univariate imputation method
imputation <- mice(subset(dataP,treat==group), m=5, method="pmm", visitSequence="monotone",
                  predictorMatrix = predictor.selection)


#extract stacked data files
long.imputation = rbind(long.imputation,complete(imputation, action="long"))

} #finish loop

#extract a single imputation dataset
imputation1 <- subset(long.imputation, subset=.imp==1)

imputation1

#create a list of all imputed datasets
library(mitools) #load package to combine imputed datasets 
#this object was specifically designed to be analyzed with the survey package
allImputations <- imputationList(list(
    subset(long.imputation, subset=.imp==1),
    subset(long.imputation, subset=.imp==2),
    subset(long.imputation, subset=.imp==3),
    subset(long.imputation, subset=.imp==4), 
    subset(long.imputation, subset=.imp==5)))


#------------------------------------
#ESTIMATE PROPENSITY SCORES WITH LOGISTIC REGRESSION    



#Create a vector of the covariates that will be balanced by using 
#propensity score methods
#define formula for propensity score model
covariateNames <-  c("X1ij",
      "X2ij", 
      "X3ij",
      "X4ij", 
      "X5ij", 
      "X6ij",
      "X7ij",
      "X8ij",
      "X9ij",
      "X10ij", 
      "X11ij", 
      "X12ij",     
      "X13ij", 
      "X14ij", 
      "X15ij", 
      "X16ij", 
      "X17ij", 
      "X18ij", 
      "X19ij", 
      "X20ij", 
      "X21ij", 
      "W1i_rep", 
      "W2i_rep", 
      "W3i_rep", 
      "W4i_rep")
      

#check whether any dummy missing indicators are redundant because
#variables have missing values for the same cases
missingCorrelations <- cor(missing.indicator[,propMissing>0.05])
diag(missingCorrelations) <- 0
maxCorrelations <- apply(missingCorrelations,2,max) 
dummyNAnames <- names(maxCorrelations)[maxCorrelations < 0.8]
maxCorrelationsHigh <- maxCorrelations[!duplicated(maxCorrelations)]
dummyNAnames <- c(dummyNAnames,names(maxCorrelationsHigh)[maxCorrelationsHigh>=0.8])



#add missing value indicators for variables with more than 5%
#of missing data to covariateNames
#merge covariate names with missing indicator names
covariateNames <- c(covariateNames, dummyNAnames)

#create formula based on covariate list
#this includes both individual level and school level covariates

psFormula <- paste(covariateNames, collapse="+")
# psFormula<- formula(imputation1, "treat~X1ij + X2ij + X3ij + X4ij + X5ij + X6ij + X7ij
#                           + X8ij + X9ij + X10ij + X11ij + X12ij
#                           + X13ij + X14ij + X15ij + X16ij + X17ij + X18ij + X19ij
#                           + X20ij + X21ij + W1i_rep + W2i_rep + W3i_rep + W4i_rep +
#                        (1|ID_group)")

psFormula <- formula(paste("treat~",psFormula, sep=""))
print(psFormula)

#define design object that describes the sample characteristics
#the variable psu identifies the primary sampling units (cluster ids)
#the variable STRATA_ID identifies the strata ids
#the variable bystuwt identifies the base-year sampling weights for
#respondents of the 2002 and 2004 waves (Base year and 1st follow-up)


#fit a logistic regression model with adjustments for clustering and stratification
#to the first imputed dataset

# Start the clock!
{ptm <- proc.time()

psModel1 <- glmer(treat~X1ij + X2ij + X3ij + X4ij + X5ij + X6ij + X7ij
                  + X8ij + X9ij + X10ij + X11ij + X12ij
                  + X13ij + X14ij + X15ij + X16ij + X17ij + X18ij + X19ij
                  + X20ij + X21ij + W1i_rep + W2i_rep + W3i_rep + W4i_rep + 1|ID_group,
                 data = imputation1, family=binomial(link = "logit"))

# Stop the clock
proc.time() - ptm}

# #fit a logistic regression model to all imputed datasets
# psModelAll <- with(imputation1, glmer(psFormula))

#extract propensity scores from first dataset
pScores <-  fitted(psModel1)
imputation1$pScores <-  pScores

# #extract propensity scores from all imputed datasets
# pScoresAll <- sapply(psModelAll, fitted)

#combine propensity scores across imputed datasets by taking the mean
# pScoresMean <- apply(pScoresAll,1,mean)

#add propensity score mean to imputed datasets
# allImputations <- update(allImputations, pScores = pScoresMean)






# #========================================================



#obtain summary statistics of common support for
#propensity scores estimated with the first imputed dataset

#for logistic regression 
with(imputation1, by(pScores,treat,summary))

#for all propensity scores
# by(imputation1[,63:65],dataP.imputed$treat,summary)

#create a table
tableCommonSupport = rbind(
    summary(imputation1[imputation1$treat==1,7:31]),
    summary(imputation1[imputation1$treat==0,7:31]))
rownames(tableCommonSupport) = c(rep("Treated",6),rep("Control",6))
write.csv(tableCommonSupport, file="Table_common_support.csv")   

#obtain proportion of treated cases above maximum control cases 
#and proportion of control cases below minum treated cases
#for logistic regression
with(imputation1, 100*c(
   mean(as.numeric(pScores[treat==1] > max(pScores[treat==0]))),
   mean(as.numeric(pScores[treat==0] < min(pScores[treat==1])))))

#obtain proportions of treated cases above maximum control cases
percentageAbove = with(imputation1, 100*c(
     mean(as.numeric(pScores[treat==1] > max(pScores[treat==0])))))
     # mean(as.numeric(pScoresRf[treat==1] > max(pScoresRf[treat==0]))),
     # mean(as.numeric(pScoresGBM[treat==1,] > max(pScoresGBM[treat==0,])))))

#obtain proportions of control cases below minimum treated cases
percentageBelow = with(imputation1, 100*c(
    mean(as.numeric(pScores[treat==0] < min(pScores[treat==1])))))
    # mean(as.numeric(pScoresRf[treat==0] < min(pScoresRf[treat==1]))),
    # mean(as.numeric(pScoresGBM[treat==0,] < min(pScoresGBM[treat==1,])))))


#evaluate common support with box and whiskers plot
imputation1$treat = factor(imputation1$treat)
library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
tiff("Chapter2_figure2-4.tif", res=600, compression = "lzw", height=6, width=15, units="in")
bwplot( pScores~treat, data = imputation1,  ylab = "Propensity Scores", xlab = "Treatment",auto.key = TRUE)
dev.off()

# 
# tiff("Chapter2_figure2-5.tif", res=600, compression = "lzw", height=6, width=15, units="in")
# bwplot( pScoresRf~treat ,  data = imputation1,  ylab = "Propensity Scores", xlab = "Treatment",auto.key = TRUE)
# dev.off()
# 
# 
# tiff("Chapter2_figure2-6.tif", res=600, compression = "lzw", height=6, width=15, units="in")
# bwplot( pScoresGBM~treat,  data = dataP.imputed,  ylab = "Propensity Scores", xlab = "Treatment",auto.key = TRUE)
# dev.off()

#evaluate common support with kernel density plots
require(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
tiff("Chapter2_figure2-4b.tif", res=600, compression = "lzw", height=6, width=15, units="in")
densityplot( ~pScores, groups=treat, plot.points=F, xlim=c(0,1), lwd=2,
             data = imputation1,  ylab = "Propensity Scores", xlab = "Treatment",auto.key = TRUE)
dev.off()


# tiff("Chapter2_figure2-5b.tif", res=600, compression = "lzw", height=6, width=15, units="in")
# densityplot( ~pScoresRf, groups=treat, plot.points=F, xlim=c(0,1), lwd=2,
#              data = dataP.imputed,  ylab = "Propensity Scores", xlab = "Treatment",auto.key = TRUE)
# dev.off()
# 
# 
# tiff("Chapter2_figure2-6b.tif", res=600, compression = "lzw", height=6, width=15, units="in")
# densityplot( ~pScoresGBM, groups=treat, plot.points=F, xlim=c(0,1), lwd=2,
#              data = dataP.imputed,  ylab = "Propensity Scores", xlab = "Treatment",auto.key = TRUE)
# dev.off()

#obtain histograms of propensity scores estimated with logistic regression
#for all inputed datasets
#first stack the imputed datasets
# allImputationsStacked <- data.frame() 
# for (imp in 1:5) { temp <- cbind(allImputations$imputations[[imp]],imputation = imp)
#                    allImputationsStacked = rbind(allImputationsStacked, temp) }
# allImputationsStacked$treat <- factor(allImputationsStacked$treat,
#                                       levels=c(0,1),labels=c("Untreated","Treated"))
# allImputationsStacked$imputation <- factor(allImputationsStacked$imputation,
#                                            labels=paste("Imputation",1:5))

# library(lattice) #used for the density plots
# lattice.options(default.theme = standard.theme(color = FALSE))
# tiff("Chapter2_figure2-7.tif", res=600, compression = "lzw", height=6, width=15, units="in")
# densityplot( ~pScores | imputation, data = allImputationsStacked, plot.points=F, lwd=2,
#              groups = treat, xlab = "Propensity Scores", auto.key = TRUE)
# dev.off()

# 
# tiff("Chapter2_figure2-8.tif", res=600, compression = "lzw", height=6, width=15, units="in")
# bwplot( pScores~treat | imputation, data = allImputationsStacked, lwd=2,
#               ylab = "Propensity Scores", auto.key = TRUE)
# dev.off()

#save results
#save(list=c("ELS.data.imputed", "allImputations","covariateNames","psFormula"), 
 # file="Chapter3_ELS_data_imputed_example_career_academy.Rdata", compress=T)

#save final results
#save(list=ls(), file="Chapter2_All_results.Rdata", compress=T)
------------------------------------------------------------------------------


#perfom full matching
library(MatchIt)  
library(optmatch)
fullmatching <- matchit(psFormula, distance= imputation1$pScores, 
                        data = imputation1, method = "nearest")
  #diagnose covariate balance
balance.fullmatching<- summary(fullmatching, standardize = T)

  #extract the table of balance
table.balance<- balance.fullmatching$sum.matched

#Estimate ATT with optimal full matched data using regression

#obtained match data
  data.fullmatching<- match.data(fullmatching)
#check number of subclasses created
  table(data.fullmatching$subclass)
  #check the weights
  table(data.fullmatching$weights)
  sum(data.fullmatching$weights)
  
#estimate the treatment effect  
  library(survey)
  design.fullmatching<- svydesign(ids=~1, wieghts =~weights,
                                 data=data.fullmatching)
#fit regression model
  model.fullmatching<- svyglm(Yij~treat, design.fullmatching, family
                                = gaussian())
  summary(model.fullmatching)
#estimate treatment effects adjusting for possible violation of independence
#due to matching of obseravtions (cluster effects)
  design.fullmatching2<- svydesign(ids=~subclass, wieghts =~weights,
                                  data=data.fullmatching)
  design.fullmatching2<- as.svrepdesign(design.fullmatching2, type = "bootstrap",
                                        replicates = 1000)
  model.fullmatching2<- svyglm(Yij~treat, design.fullmatching2, family
                              = gaussian())
  summary(model.fullmatching2)
  
  #the ENd!
