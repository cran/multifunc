## ----set-options, echo=FALSE, cache=FALSE, include=FALSE----------------------
library(knitr)
opts_chunk$set(tidy=FALSE,
               prompt=FALSE,
               warning=FALSE,
               message = FALSE)
opts_chunk$set(comment="#")

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("multifunc", username="jebyrnes")

## ---- warning=FALSE, results="hide", message=FALSE, error=FALSE---------------
library(multifunc)

#for plotting
library(ggplot2)
library(patchwork)

#for data
library(tidyr)
library(dplyr)
library(purrr)
library(forcats)

#for analysus
library(car)

## -----------------------------------------------------------------------------
#Read in data  to run sample analyses on the biodepth data
data(all_biodepth)

allVars<-c("biomassY3", "root3", "N.g.m2",  "light3", "N.Soil", "wood3", "cotton3")

varIdx<-which(names(all_biodepth) %in% allVars)

## -----------------------------------------------------------------------------
#######
# Now, specify what data we're working with - Germany
# and the relevant variables we'll be working with
#######

germany<- all_biodepth %>%
  filter(location=="Germany")


# which variables have > 2/3 of the values not NA?
german_vars <- whichVars(germany, allVars)

# What are the names of species in this dataset
# that have at least some values > 0?
species <- relevantSp(germany,26:ncol(germany))

spIDX <- which(names(germany) %in% species) #in case we need these

## -----------------------------------------------------------------------------
germanyForPlotting <- germany %>%
  select(Diversity, german_vars) %>%
  pivot_longer(cols = german_vars,
               names_to = "variable",
               values_to = "value") %>%
  mutate(variable = factor(variable))


#make the levels of the functions into something nice for plots
levels(germanyForPlotting$variable) <- c('Aboveground Biomass', 'Root Biomass', 'Cotton Decomposition', 'Soil Nitrogen', 'Plant Nitrogen')

## -----------------------------------------------------------------------------
germanyLabels <- germanyForPlotting %>%
  group_by(variable) %>%
  nest() %>%
  mutate(fits = map(data, ~lm(value ~ Diversity, data=.x)),
         stats = map(fits, broom::glance)) %>%
  unnest(stats) %>%
  ungroup() %>%
  mutate(labels = paste0("p = ", round(p.value, 3), ", ", 
                        expression(R^2), " = ", round(r.squared, 2)),
         labels =  gsub("p = 0 ", "p < 0.001 ", labels),
         Diversity = 7, 
         value=c(2000, 1200, 1.5, 6, 20))

## ---- fig.width=12, fig.height=8----------------------------------------------
ggplot(aes(x=Diversity, y=value),data=germanyForPlotting) +
  geom_point(size=3)+
  facet_wrap(~variable, scales="free") +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) + 
  xlab("Species Richness") +
  ylab("Value of Function") +
  geom_text(data=germanyLabels,
            aes(label=labels)) +
 theme(panel.grid = element_blank())

## -----------------------------------------------------------------------------
#pull out only those plots planted in monoculture
monoGermany<- germany %>%
  filter(rowSums(select(., ACHMIL1:VICTET1))==1) %>%
  #get mono names
  pivot_longer(ACHMIL1:VICTET1,
               names_to = "mono") %>%
    filter(value !=0) %>%
  select(-value) %>%
#pivot by function
  pivot_longer(german_vars,
               names_to = "variable") %>%
  
#get the mean for each monoculture
group_by(variable, mono) %>%
summarize(value = mean(value))
  
#now, who is the best performer?
monoGermany %>%
  group_by(variable) %>%
  filter(value == max(value)) 

## -----------------------------------------------------------------------------
germany <- germany %>%
  mutate(N.Soil = -1*N.Soil +max(N.Soil, na.rm=T))

## -----------------------------------------------------------------------------
spList<-sAICfun("biomassY3", species, germany)
spList

## -----------------------------------------------------------------------------

redund<-getRedundancy(german_vars, species, germany)
coefs<-getRedundancy(german_vars, species, germany, output="coef")
stdCoefs<-stdEffects(coefs, germany, german_vars, species)

#for example
redund

## ----fig.height=8, fig.width=12-----------------------------------------------

#plot the num. functions by fraction of the species pool needed
posCurve<-divNeeded(redund, type="positive")

posCurve$div<-posCurve$div/ncol(redund)

pc<-qplot(nfunc, div, data=posCurve, group=nfunc, geom=c("boxplot"))+
  geom_jitter(size=4, position = position_jitter(height=0.001, width = .04))+
  ylab("Fraction of Species Pool\nwith Positive Effect\n")+ 
  xlab("\nNumber of Functions")+theme_bw(base_size=24)+ylim(c(0,1))


negCurve<-divNeeded(redund, type="negative")

negCurve$div<-negCurve$div/ncol(redund)

nc<-qplot(nfunc, div, data=negCurve, group=nfunc, geom=c("boxplot"))+
  geom_jitter(size=4, position = position_jitter(height=0.001, width = .04))+
  ylab("Fraction of Species Pool\nwith Negative Effect\n")+ 
  xlab("\nNumber of Functions")+theme_bw(base_size=24)+ylim(c(0,1))

#combine these into one plot
pc + nc +
   plot_annotation(tag_levels = 'A')

## -----------------------------------------------------------------------------
allOverlapPos <- map_df(2:nrow(redund),
                        ~getOverlapSummary(redund, m = .x, denom = "set"),
                        .id = "m")

allOverlapPos

## -----------------------------------------------------------------------------
posNeg<-data.frame(Species = colnames(redund),
                   Pos = colSums(filterOverData(redund)),
                   Neg = colSums(filterOverData(redund, type="negative")))

#plot it
ggplot(aes(x=Pos,y= Neg), data=posNeg) +
  geom_jitter(position = position_jitter(width = 0.05, height = 0.05), size=5, shape=1) +
  theme_bw(base_size=18) +
  xlab("\n# of Positive Contributions per species\n") + 
  ylab("# of Negative Contributions per species\n") +
  annotate("text", x=0, y=3, label="a)") +
  stat_smooth(method="glm", colour="black", size=2, 
              method.args = list(family=poisson(link="identity"))) +
  geom_abline(slope = 1, intercept = 0, lty = 2, color = "red")

## -----------------------------------------------------------------------------
#now get the standardized effect sizes
posNeg<-within(posNeg, {
  
  stdPosMean <- colSums(filterCoefData(stdCoefs))/Pos
  stdPosMean[which(is.nan(stdPosMean))] <-0
  
  stdNegMean <- colSums(filterCoefData(stdCoefs, type="negative"))/Neg
  stdNegMean[which(is.nan(stdNegMean))] <-0  
})


ggplot(aes(x=stdPosMean,y= stdNegMean), data=posNeg) +
  #  geom_point(size=3) +
  geom_jitter(position = position_jitter(width = .02, height = .02), size=5, shape=1) +
  theme_bw(base_size=18) +
  xlab("\nAverage Standardized Size of\nPositive Contributions") + 
  ylab("Average Standardized Size of\nNegative Contributions\n") +
  stat_smooth(method="lm", colour="black", size=2)+
  annotate("text", x=0, y=0.25, label="b)") +
  geom_abline(slope=-1, intercept=0, size=1, colour="black", lty=3)

## -----------------------------------------------------------------------------

#add on the new functions along with the averaged multifunctional index
germany<-cbind(germany, getStdAndMeanFunctions(germany, german_vars))
#germany<-cbind(germany, getStdAndMeanFunctions(germany, vars, standardizeZScore))

## -----------------------------------------------------------------------------

#plot it
ggplot(aes(x=Diversity, y=meanFunction),data=germany)+geom_point(size=3)+
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) + 
  xlab("\nSpecies Richness") +
  ylab("Average Value of Standardized Functions\n")

## ----fig.width=16, fig.height=8-----------------------------------------------
#reshape for plotting everything with ggplot2
germanyMeanForPlotting <- germany %>%
  select(Diversity, biomassY3.std:meanFunction) %>%
  pivot_longer(cols = c(biomassY3.std:meanFunction),
               names_to = "variable") %>%
    mutate(variable = fct_inorder(variable))
#nice names for plotting
levels(germanyMeanForPlotting$variable) <- c('Aboveground Biomass', 'Root Biomass', 'Cotton Decomposition', 'Soil Nitrogen', 'Plant Nitrogen', 'Mean Multifuncion Index')

#plot it
ggplot(aes(x=Diversity, y=value),data=germanyMeanForPlotting)+geom_point(size=3)+
  facet_grid(~variable) +
  theme_bw(base_size=15)+
  stat_smooth(method="lm", colour="black", size=2) + 
  xlab("\nSpecies Richness") +
  ylab("Standardized Value of Function\n")

## -----------------------------------------------------------------------------
#statistical fit
aveFit<-lm(meanFunction ~ Diversity, data=germany)
Anova(aveFit)
summary(aveFit)

## -----------------------------------------------------------------------------
germanyThresh<-getFuncsMaxed(germany, german_vars, threshmin=0.05, threshmax=0.99, prepend=c("plot","Diversity"), maxN=7)

## -----------------------------------------------------------------------------
mfuncGermanyLinear08<-glm(funcMaxed ~ Diversity, data=subset(germanyThresh, germanyThresh$thresholds=="0.8"), family=quasipoisson(link="identity"))

Anova(mfuncGermanyLinear08, test.statistic="F")
summary(mfuncGermanyLinear08)

## -----------------------------------------------------------------------------

gcPlot<-subset(germanyThresh, germanyThresh$thresholds %in% qw(0.2, 0.4, 0.6, 0.8)) #note, using qw as %in% is a string comparison operator

gcPlot$percent<-paste(100*gcPlot$thresholds, "%", sep="")

qplot(Diversity, funcMaxed, data=gcPlot, facets=~percent) +
  stat_smooth(method="glm", 
              method.args = list(family=quasipoisson(link="identity")),
              colour="red", lwd=1.2) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("Species Richness") +
  theme_bw(base_size=14) +
  geom_text(data=data.frame(percent = unique(gcPlot$percent),
                       lab = paste(letters[1:4], ")", sep=""),
                       Diversity=2,
                       funcMaxed=6
                       ), mapping=aes(x=Diversity, y=funcMaxed, label=lab))

## -----------------------------------------------------------------------------

germanyThresh$percent <- 100*germanyThresh$thresholds
ggplot(data=germanyThresh, aes(x=Diversity, y=funcMaxed, group=percent)) +
    ylab(expression("Number of Functions" >= Threshold)) +
    xlab("Species Richness") +
    stat_smooth(method="glm", 
                method.args = list(family=quasipoisson(link="identity")), 
                lwd=0.8, fill=NA, aes(color=percent)) +
    theme_bw(base_size=14) +
    scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red")

## ---- fig.height=8, fig.width=12----------------------------------------------
germanyLinearSlopes<-getCoefTab(funcMaxed ~ Diversity,
                                data = germanyThresh, 
                                coefVar = "Diversity",
                                family = quasipoisson(link="identity"))


######
# Plot the values of the diversity slope at
# different levels of the threshold
######
germanSlopes <- ggplot(germanyLinearSlopes, aes(x=thresholds*100,
                                                y = estimate,
                                                ymax = estimate + 1.96 * std.error,
                                                ymin = estimate - 1.96 * std.error)) +
  geom_ribbon(fill="grey50") +
  geom_point() +
  ylab("Change in Number of Functions per Addition of 1 Species\n") +
  xlab("\nThreshold (%)") +
  geom_abline(intercept=0, slope=0, lwd=1, linetype=2) +
  theme_bw(base_size=14)

germanSlopes

## -----------------------------------------------------------------------------
germanIDX <- getIndices(germanyLinearSlopes, germanyThresh, funcMaxed ~ Diversity)
germanIDX

## ---- fig.height=8, fig.width=12----------------------------------------------

germanyLinearSlopes$estimate[which(germanyLinearSlopes$thresholds==germanIDX$Tmde)]

germanyThresh$IDX <- 0
germanyThresh$IDX [which(germanyThresh$thresholds %in% 
                           c(germanIDX$Tmin, germanIDX$Tmax, germanIDX$Tmde))] <- 1



ggplot(data=germanyThresh, aes(x=Diversity, y=funcMaxed, group=percent)) +
  ylab(expression("Number of Functions" >= Threshold)) +
  xlab("Species Richness") +
  geom_smooth(method="glm", 
              method.args = list(family=quasipoisson(link="identity")), 
              fill=NA, aes(color=percent, lwd=IDX)) +
  theme_bw(base_size=14) +
  scale_color_gradient(name="Percent of \nMaximum", low="blue", high="red") +
  scale_size(range=c(0.3,5), guide="none") +
  annotate(geom="text", x=0, y=c(0.2,2,4.6), label=c("Tmax", "Tmde", "Tmin")) +
  annotate(geom="text", x=16.7, y=c(germanIDX$Mmin, germanIDX$Mmax, germanIDX$Mmde), label=c("Mmin", "Mmax", "Mmde"))


## ---- fig.height=8, fig.width=12----------------------------------------------
germanSlopes + annotate(geom="text", y=c(-0.01, -0.01, -0.01, germanIDX$Rmde.linear+0.02), x=c(germanIDX$Tmin*100, germanIDX$Tmde*100, germanIDX$Tmax*100, germanIDX$Tmde*100),  label=c("Tmin", "Tmde", "Tmax", "Rmde"), color="black") 
  
  

## ---- error=FALSE, warning=FALSE, message=FALSE-------------------------------
####
# Now we will look at the entire BIODEPTH dataset
#
# Note the use of dplyr to generate the new data frame using all locations from biodepth.
# If you're not using dplyr for your data aggregation, you should be!
# It will save you a lot of time and headaches in the future.
####

#Read in data  to run sample analyses on the biodepth data
data(all_biodepth)
sub_biodepth<-subset(all_biodepth, all_biodepth$location %in% c("Sheffield", "Portugal", "Sweden"))
sub_biodepth$location<-factor(sub_biodepth$location, levels=c("Portugal", "Sweden", "Sheffield"))

allVars<-qw(biomassY3, root3, N.g.m2,  light3, N.Soil, wood3, cotton3)
varIdx<-which(names(sub_biodepth) %in% allVars)



#re-normalize so that everything is on the same 
#sign-scale (e.g. the maximum level of a function is 
#the "best" function)
sub_biodepth <- sub_biodepth %>%
  group_by(location) %>%
  mutate(
    light3 <- ifelse(sum(is.na(light3)==0),
      -1*light3+max(light3, na.rm=T),
      NA),
    
    N.Soil <- ifelse(sum(is.na(N.Soil)==0),
      -1*N.Soil+max(N.Soil, na.rm=T),
      NA)
  )  %>%
  ungroup()


#get thresholds
bdThreshes<-sub_biodepth %>%
  group_by(location) %>%
  nest() %>%
  mutate(fmaxed = map(data, getFuncsMaxed,
                         vars=allVars,
                         maxN=8)) %>%
  ungroup() %>%
  unnest(fmaxed)



####look at slopes

#note, maxIT argument is due to some fits needing more iterations to converge
bdLinearSlopes<-getCoefTab(funcMaxed ~ Diversity,
                           data=bdThreshes,
                           groupVar=c("location", "thresholds"), 
                           coefVar="Diversity",
                           family=quasipoisson(link="identity"),
                           control=list(maxit=800))


## ---- error=FALSE, warning=FALSE, message=FALSE-------------------------------

indexTable <- lapply(levels(bdLinearSlopes$location), function(x){
  slopedata <- subset(bdLinearSlopes, bdLinearSlopes$location==x) 
  threshdata <- subset(bdThreshes, bdThreshes$location==x) 
  ret <- getIndices(slopedata, threshdata, funcMaxed ~ Diversity)
  ret<-cbind(location=x, ret)
  ret
})
indexTable <- do.call(rbind, indexTable)

indexTable <- rbind(data.frame(location="Germany", germanIDX), indexTable)


indexTable

## ---- fig.height=8, fig.width=10----------------------------------------------
library(gridExtra)

bdCurvesAll<-qplot(Diversity, funcMaxed, data=bdThreshes, group=thresholds, alpha=I(0)) +
  facet_wrap(~location)+#, scales="free") +
  scale_color_gradient(low="blue", high="red", name="Proportion of \nMaximum", guide=FALSE) +
  stat_smooth(method="glm", lwd=0.8, fill=NA,
              method.args = list(family=gaussian(link="identity"), 
                                 control=list(maxit=200)), aes(color=thresholds)) +
  ylab("\nNumber of Functions ≥ Threshold\n\n") +
  xlab("Species Richness") +
  theme_bw(base_size=15)  + 
  geom_text(data=data.frame(location = levels(bdThreshes$location), lab = paste(letters[1:3], ")", sep=""), thresholds=c(NA, NA, NA)), x=1, y=7, mapping=aes(label=lab))



#Plot it!
slopePlot<-ggplot(bdLinearSlopes, aes(x=thresholds*100, 
                                      y=estimate,
                                      ymin = estimate - 2*std.error,
                                      ymax = estimate + 2*std.error)) +
  geom_ribbon(fill="grey50") +
  geom_point() +
  ylab("Change in Number of Functions \nper Addition of 1 Species\n") +
  xlab("Threshold (%)") +
  facet_wrap(~location)+#, scale="free") +
  geom_abline(intercept=0, slope=0, lwd=0.6, linetype=2) +
  theme_bw(base_size=15)+ 
  geom_text(data=data.frame(location = levels(bdThreshes$location), 
                            lab = paste(letters[4:6], ")", 
                                        sep=""), 
                            thresholds=c(NA, NA, NA),
                            std.error=c(NA, NA, NA),
                            estimate = c(NA, NA, NA)
                            ), 
            x=0.05, y=0.27, mapping=aes(label=lab))


###Plot them in a single figure
grid.arrange(bdCurvesAll,
             slopePlot) 

