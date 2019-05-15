# All code for cost-effectiveness paper



### Example 1: BRCA

##
# Plot MRS and INB for WAS data
# HKatki 10 7 16
##

##
# Read excel sheet of all BRCA information. Each row is a risk-threshold: sens, spec, PPV, NPV, MRS, etc. 
##
setwd("/Users/katkih/hkCurrent/Dropbox/RiskStratification/MRS/MRSmethodology/OptimalThreshold1TimeScreen")
library(readxl)
brca <- as.data.frame(read_excel("data/BRCA.xlsx", "Thresholds"))
brcaGenPop <- as.data.frame(read_excel("data/BRCA_GenPop.xlsx"))
# brca <- as.data.frame(read_excel("E://Staff/Ionut Bebu/METHODS/Optimal_Screening/Screening_Hormuzd/DATA/BRCA.xlsx", "Thresholds"))
# brcaGenPop <- as.data.frame(read_excel("E://Staff/Ionut Bebu/METHODS/Optimal_Screening/Screening_Hormuzd/DATA/BRCA_GenPop.xlsx"))

##
# Set path to the LaTeX file, from which the figures/ subdirectoty is under.
# For weird security reasons, LaTeX will not let you go back from the path of the LaTeX file.  
# This is why I can't have a figures/ subdirectory outside of /docs that I would prefer.
##
setwd(paste(getwd(),"/doc/StatMed/ama/figures",sep=""))


##
# Helper Function to calculate deltas for Cost and Effectiveness, and Incremental Net Benefit (INB)
##
INB <- function(E,C,probDplusMplus,probMplus,probDplus) {
  delE <- (E[2]-E[3])*probDplusMplus
  delC <- C[1] + C[2]*probMplus + (C[3]-C[4]-C[2])*probDplusMplus
  INB <- delE - delC
  e <- (E[2]-E[3])-(C[3]-C[4])
  INBall <- INB - ((e+C[2])*probDplus-C[2])
  cbind(delE,delC,INB,INBall)
}    



##
# Fix effectiveness and costs
##
dollarsperLYG <- 51207 # dollars per year of life
E <- c(80,78,73)*dollarsperLYG

# Manchanda costs
# C_AJ <- c(100,250,14558,19615) # BRCAPRO test cost, BRCA test cost, cost of prevention & early treatment, cost of late treatment
# C_GenPop <- c(100,1000,14558,19615*0.85) # BRCAPRO test cost, BRCA test cost, cost of prevention & early treatment, cost of late treatment

# Long and Ganz costs: initial + during + last year costs
breastcancercosts <- 86013+8*7547+63790 # average 10 years of surviving before death by breast cancer
ovariancancercosts <- 124838+3*13724+87218 # average 5 years of surviving
# RRM+RRSO+cancer treatment costs.  Assume yearly breast cancer surveillance costs the same for C_1 and C_2, so can ignore
C_1 <- 12286+7393 + 0.027*breastcancercosts + 0.012*ovariancancercosts 
# breast treatment avg cost, ovarian treatment avg cost, zero cost if no cancer
C_2 <- 0.53*breastcancercosts+0.299*ovariancancercosts+ 0.171*0 
C_AJ <- c(100,249,C_1,C_2) # BRCAPRO test cost, BRCA test cost, cost of prevention & early treatment, cost of late treatment
C_GenPop <- c(100,2200,C_1,C_2) # BRCAPRO test cost, BRCA test cost, cost of prevention & early treatment, cost of late treatment

# Net Effectiveness and risk thresholds
e_AJ <- (E[2]-E[3])-(C_AJ[3]-C_AJ[4]) # Net Effectiveness of early intervention
e_GenPop <- (E[2]-E[3])-(C_GenPop[3]-C_GenPop[4]) # Net Effectiveness of early intervention
# Optimal risk thresholds
R_AJ <- C_AJ[2]/(e_AJ+C_AJ[2])
R_GenPop <- C_GenPop[2]/(e_GenPop+C_GenPop[2])



##
# Calculate costs and effectivenesses for Ashkenazi-Jews
##
probDplusMplus <- brca$MRS/2 + brca$OverallPtest*brca$BRCAprev
probDplus <- brca$BRCAprev[2]
allAJ <- data.frame(brca$Threshold,brca$MRS,INB(E,C_AJ,probDplusMplus,brca$OverallPtest,probDplus))
colnames(allAJ)[1:2] <- c("Threshold","MRS")
INBallAJ <- (e_AJ+C_AJ[2])*probDplus-C_AJ[2]
allAJ$INBgain <- allAJ$INB - max(0,INBallAJ)
probDplusMminus <- probDplus - probDplusMplus
probDminusMplus <- brca$OverallFP*(1-probDplus)
probDminusMminus <- (1-brca$OverallFP)*(1-probDplus)

########
# FOR AJs
# This shows that that the INB is very sensitive to the cost of the screening BRCAPRO test
# If it were costless, there are a range of risk thresholds for screening, but at $100, screening is never cost-worthy
# When screening is never cost-worthy, which of the all-or-nothing actions are best?
# It always flips at optimal risk threshold equals prevalence: below prevalence, test everyone, above prevalence, test no one.
#########
# Numerically calculate the bounds on dollars per life-year gained so that screening at optimal threshold is best decision
# This means calculate INBgain as dollarsperLYG increases

allAJ$perLYG <- ((C_AJ[2]/allAJ$Threshold - C_AJ[2]) + (C_AJ[3]-C_AJ[4])) / ((E[2]-E[3])/dollarsperLYG)
allAJ$e <- (E[2]-E[3])/dollarsperLYG * allAJ$perLYG - (C_AJ[3]-C_AJ[4])
# The reason I call these "new" to distinguish from the $INB as calculated in OLD WORK below, which is irrelevant now 
#allAJ$newINB <- (allAJ$e+C_AJ[2])*brca$OverallTP*brca$BRCAprev[2] - C_AJ[1] - C_AJ[2]*brca$OverallPtest
allAJ$newINB <- allAJ$e*probDplusMplus - C_AJ[1] - C_AJ[2]*probDminusMplus
allAJ$newINBall <- allAJ$e*probDplus - C_AJ[2]*(1-probDplus)
allAJ$newINBgain <- allAJ$newINB - ifelse(allAJ$newINBall<0,0,allAJ$newINBall)
allAJ$NB_M <- allAJ$delE - allAJ$delC # Net Benefit for non-optimal threshold NB_M, lambda is fixed dollarsperLYG

# Calculate the Net Benefits if each threshold were the optimal threshold:use the perLYG associated with that threshold
e <- E/dollarsperLYG # get rid of willingness-to-pay lambda
E_M <- e[1]*(1-probDplus) + e[2]*probDplusMplus + e[3]*probDplusMminus
C_M <- C_AJ[1] + C_AJ[2]*(1-probDminusMminus) + C_AJ[3]*probDplusMplus + C_AJ[4]*probDplusMminus
allAJ$NB_MofR <- allAJ$perLYG*E_M - C_M # Net Benefit for optimal threshold NB_M(R), where perLYG=lambda(R)
allAJ$NB_none <- allAJ$perLYG*e[1]*(1-probDplus) + (allAJ$perLYG*e[3]-C_AJ[2]-C_AJ[4])*probDplus
allAJ$NB_all  <- (allAJ$perLYG*e[1]-C_AJ[2])*(1-probDplus) + (allAJ$perLYG*e[2]-C_AJ[2]-C_AJ[3])*probDplus
allAJ$INB_gain <- allAJ$NB_MofR - pmax(allAJ$NB_none,allAJ$NB_all)

# Plot the NBs, but they are so close to each other, you can't see the difference.  
# Instead, show this as a table below.
# par(mfrow=c(1,1),pty="m")
# plot(allAJ$Threshold,allAJ$NB_MofR,type="l",log="x")
# lines(allAJ$Threshold,allAJ$NB_none,lty=2,col="blue")
# lines(allAJ$Threshold,allAJ$NB_all,lty=3,col="red")




##
# Make all plots for Ashkenazi Jews
##

##
# Figure 1: Plot of Net Benefit vs. (non-optimal) risk thresholds.  This is the simple view.
##
par(mfrow=c(1,2),pty="s")

plot(allAJ$Threshold, allAJ$delE - allAJ$delC,xlim=c(1.6e-04,1.3),log="x",type="l",lty=1,
     xlab="Risk Threshold",ylab="Net Benefit for screening Ashkenazi-Jewish women",
     main=paste("lambda=$",dollarsperLYG,"per life-year gained"))
# Net Benefit for testing everyone
points(1.6e-04,allAJ$delE[2]-allAJ$delC[2]+C_AJ[1],pch=8)
text(1.7e-04,allAJ$delE[2]-allAJ$delC[2],"test\nall",cex=0.9,pos=1)
# Net Benefit for testing no one
points(1,allAJ$delE[4970]-allAJ$delC[4970]+C_AJ[1],pch=8) 
text(1,allAJ$delE[4970]-allAJ$delC[4970]+250,"test\nno\none",cex=0.9,adj=c(0.25,0))
# Optimal threshold
points(R_AJ,allAJ$delE[4]-allAJ$delC[4],pch=8) 
text(R_AJ+0.0003,allAJ$delE[4]-allAJ$delC[4]-100,"optimal\nthreshold",cex=0.9,pos=1)

# Now plot INBgain to compare to NBgain
plot(brca$Threshold,allAJ$newINBgain,log="x",xlim=c(0.005,0.5),ylim=c(-150,0),type="l",lty=1,
     xlab="Optimal Risk Threshold (R)",ylab="Incremental Net Benefit gain",main="INB gain")
abline(h=0,lty=3)
abline(v=brca$BRCAprev[2],lty=3)
text(brca$BRCAprev[2],-150,"If INB_Gain<0, then test none",cex=0.9,pos=4)
text(brca$BRCAprev[2],-150,"If INB_Gain<0,\nthen test all",cex=0.9,pos=2)

# Save this figure
filename <- "Fig1AJ.pdf" ; dev.copy2pdf(file=filename)
embedFonts(file=filename,outfile=filename) #Needed for TexStudio to read fonts (looks fine in pdf file w/o needed to embed)


##
# Figure 2: Make a table of this output to explain INB_Gain
##
allAJ[2:10,c(1,8,14,15,16,17)]


##
# WebAppendix Figure: Plot INB for screening vs. INB_none and vs. INB_all in separate plots
##
par(mfrow=c(1,3),pty="m")
plot(brca$Threshold,allAJ$NB_MofR-allAJ$NB_none,log="x",xlim=c(0.005,0.5),ylim=c(-150,500),type="l",lty=1,
     xlab="Optimal Risk Threshold (R)",ylab="INB_none",main="INB_none")
abline(h=0) ; abline(v=brca$BRCAprev[2],lty=3)

plot(brca$Threshold,allAJ$NB_MofR-allAJ$NB_all,log="x",xlim=c(0.005,0.5),ylim=c(-150,500),type="l",lty=1,
     xlab="Optimal Risk Threshold (R)",ylab="INB_all",main="INB_all",col="blue")
abline(h=0); abline(v=brca$BRCAprev[2],lty=3)

# Now plot INBnone and INBall in the same plot and highlight which part is the INBgain 
plot(brca$Threshold,allAJ$newINB,type="l",lty=1,log="x",xlim=c(0.005,0.5),ylim=c(-150,500),
     xlab="Optimal Risk Threshold (R)",ylab="INB_gain",main="INB_gain")
lines(brca$Threshold,allAJ$newINB-allAJ$newINBall,lty=1,col="blue")
lines(brca$Threshold,allAJ$newINBgain,lty=1,col="red")
abline(h=0); abline(v=brca$BRCAprev[2],lty=3)
legend(0.03,250,lty=1,col=c("black","blue","red"),bty="n",legend=c("INB_none","INB_all","INB_gain"),cex=1)

# Save this figure
filename <- "INBnonevsINBallvsINBgain.pdf" ; dev.copy2pdf(file=filename)
embedFonts(file=filename,outfile=filename) #Needed for TexStudio to read fonts (looks fine in pdf file w/o needed to embed)



##
# Figure 3: NB gain from DCA versus risk threshold and second axis of willingness to pay lambda
##
par(mar=c(7,4,4,2)+0.1) #mgp=c(2.2,1,0)) # Make space at bottom for 2nd x-axis for lambda
par(mfrow=c(1,1),pty="m")

brca$NBgaincostless <- brca$NetBenefit-pmax(0,brca$NBalwaysact)
plot(brca$Threshold,brca$NBgaincostless,log="x",xaxt="n",xlim=c(0.0002,0.5),type="l",ylim=c(-0.01,0.01),lty=2,
     xlab="",ylab="Net Benefit gain from DCA",main="NB gain from DCA")
abline(h=0)
abline(v=brca$BRCAprev[2],lty=3)
text(0.05,-0.01,"If NB_Gain<0,\nthen test none",cex=0.9,pos=4)
text(brca$BRCAprev[2],-0.01,"If NB_Gain<0,\nthen test all",cex=0.9,pos=2)
# add 1st x-axis for optimal risk thresholds R
ticks <- c(0.0002,0.0004,0.001,0.003,0.01,0.02,0.1,0.5)
axis(1, ticks, labels=ticks,line=0) ; mtext("Optimal Risk Threshold (R)",1,line=2)

# Define the NBgain for a costing screening test
brca$NBgainwithcost <- brca$NBgaincostless-C_AJ[1]/allAJ$e
lines(brca$Threshold,brca$NBgainwithcost,lty=1)
legend(2e-04,0.011,lty=c(2,1),bty="n",
       legend=c("Ignoring screening test cost","With screening test cost"),cex=1)

# add 2nd x-axis for willingness to pay lambda
lambdas <- c("+$217207","+$92706","+$18007","-$15193","-$26813","-$29303","-$31295","-$31694")
axis(1, ticks, labels=lambdas,line=4,cex.axis=1) ; mtext("Willingness to pay (lambda)",1,line=6)

# Save this figure
filename <- "testcostsAJ.pdf" ; dev.copy2pdf(file=filename)
embedFonts(file=filename,outfile=filename) #Needed for TexStudio to read fonts (looks fine in pdf file w/o needed to embed)

##
# Ignore all below stuff that is commented out
##

# # Alternate Figure 1 plots on separate figures:
# # Plot NBgains, w/ and w/o screening test costs, first vs. optimal risk threshold, second vs. dollars per LYG
# # Define the usual NBgain for a costless screening test
# brca$NBgaincostless <- brca$NetBenefit-pmax(0,brca$NBalwaysact)
# par(mfrow=c(1,2),pty="s")
# plot(brca$Threshold,brca$NBgaincostless,log="x",xlim=c(0.005,0.5),type="l",ylim=c(-0.01,0.01),lty=2,
#      xlab="Optimal Risk Threshold (R)",ylab="Net Benefit gain from DCA",main="NB gain from DCA")
# abline(h=0)
# abline(v=brca$BRCAprev[2],lty=3)
# text(0.05,-0.01,"If NB_Gain<0,\nthen test none",cex=0.9,pos=4)
# text(brca$BRCAprev[2],-0.01,"If NB_Gain<0,\nthen test all",cex=0.9,pos=2)
# 
# # Define the NBgain for a costing screening test
# brca$NBgainwithcost <- brca$NBgaincostless-C_AJ[1]/allAJ$e
# lines(brca$Threshold,brca$NBgainwithcost,lty=1)
# legend(0.04,0.011,lty=c(2,1),bty="n",
#        legend=c("Ignoring screening test cost","With screening test cost"),cex=0.7)
# 
# # Plot usual NBgain for a costless screening test versus dollars per LYG
# # This version is for Manchanda costs
# # plot(allAJ$perLYG,brca$NBgaincostless,xlim=c(100,1e5),type="l",ylim=c(-0.001,0.01),lty=2,log="x",
# #      xlab="Dollars per life-year gained",ylab="Net Benefit gain",main="NB gain vs. dollars\n per life-year gained")
# # abline(h=0)
# # legend(5e2,0.01,lty=c(2),bty="n",
# #        legend=c("Ignoring screening test cost"),cex=0.8)
# 
# # Plot usual NBgain for a costless screening test versus dollars per LYG
# # This version is for Long and Ganz costs
# plot(allAJ$perLYG,brca$NBgaincostless,xlim=c(-3.1e4,1e4),type="l",ylim=c(-0.01,0.01),lty=2,
#      xlab="Dollars per life-year gained (lambda)",ylab="Net Benefit gain from DCA",main="NB gain from DCA")
# abline(h=0)
# legend(-2e4,0.01,lty=c(2),bty="n",legend=c("Ignoring screening test cost"),cex=0.8)
# abline(v=allAJ$perLYG[115],lty=3) # the dollars per LYG at prevalence: -29609
# text(allAJ$perLYG[115],-0.01,"If NB_Gain<0,\nthen test all",cex=0.9,pos=4)

# # Plot INBgain versus dollars per LYG.  INB gain is negative everywhere, including at negative dollars per LYG.
# # This version is for Long and Ganz costs
# plot(allAJ$perLYG,allAJ$newINBgain,xlim=c(-3.5e4,1e3),type="l",ylim=c(-200,0),lty=1,
#      xlab="Dollars per life-year gained",ylab="INB gain",main="INB gain vs. dollars\n per life-year gained")
# abline(h=0)

# Save this figure
# filename <- "figures/testcostsAJ.pdf" ; dev.copy2pdf(file=filename)
# embedFonts(file=filename,outfile=filename) #Needed for TexStudio to read fonts (looks fine in pdf file w/o needed to embed)

# # Define the usual NBgain for a costless screening test - I don't plot this since I don't see a use for it
# brca$NBgaincostless <- brca$NetBenefit-pmax(0,brca$NBalwaysact)
# plot(allAJ$perLYG,brca$NBgaincostless,xlim=c(100,1e5),type="l",ylim=c(-0.001,0.01),lty=2,log="x",
#      xlab="Dollars per life-year gained",ylab="Net Benefit gain",main="Ashkenazi Jews")
# abline(h=0)
# legend(5e2,0.01,lty=c(2),bty="n",
#       legend=c("Ignoring screening test cost"),cex=0.8)

# Define the NBgain for a costing screening test - I don't plot this because there's nothing interesting to say here.
# brca$NBgainwithcost <- brca$NBgaincostless-C_AJ[1]/allAJ$e
#lines(allAJ$perLYG,brca$NBgainwithcost)
# legend(5e2,-0.005,lty=c(1,2),bty="n",
#        legend=c("With screening test cost","Ignoring screening test cost"),cex=0.8)




#####
## ALL WORK FOR THE GENERAL POPULATION
###

##
# Calculate costs and effectivenesses for General-Population
##
probDplusMplus <- brcaGenPop$MRS/2 + brcaGenPop$OverallPtest*brcaGenPop$BRCAprev
probDplus <- brcaGenPop$BRCAprev[2]
probDplusMminus <- probDplus - probDplusMplus
probDminusMplus <- brcaGenPop$OverallFP*(1-probDplus)
probDminusMminus <- (1-brcaGenPop$OverallFP)*(1-probDplus)
allGenPop <- data.frame(brcaGenPop$Threshold,brcaGenPop$MRS,INB(E,C_GenPop,probDplusMplus,brcaGenPop$OverallPtest,probDplus))
colnames(allGenPop)[1:2] <- c("Threshold","MRS")
INBallGenPop <- (e_GenPop+C_GenPop[2])*probDplus-C_GenPop[2]
allGenPop$INBgain <- allGenPop$INB - max(0,INBallGenPop)

# Numerically calculate the bounds on dollars per life-year gained so that screening at optimal threshold is best decision
# This means calculate INBgain as dollarsperLYG increases
allGenPop$perLYG <- ((C_GenPop[2]/allGenPop$Threshold - C_GenPop[2]) + (C_GenPop[3]-C_GenPop[4])) / ((E[2]-E[3])/dollarsperLYG)
allGenPop$e <- (E[2]-E[3])/dollarsperLYG * allGenPop$perLYG - (C_GenPop[3]-C_GenPop[4])
# The reason I call these "new" is to distinguish from the $INB as calculated in OLD WORK below, which is irrelevant now 
allGenPop$newINB <- (allGenPop$e+C_GenPop[2])*brcaGenPop$OverallTP*brcaGenPop$BRCAprev[2] - C_GenPop[1] - C_GenPop[2]*brcaGenPop$OverallPtest
allGenPop$newINBall <- (allGenPop$e+C_GenPop[2])*brcaGenPop$BRCAprev[2] - C_GenPop[2]
allGenPop$newINBgain <- allGenPop$newINB - ifelse(allGenPop$newINBall<0,0,allGenPop$newINBall)

##
# Sensitivity Analysis Section 4.2: Sensitivity to cost of the definitive test (genetic test cost)
# This calculates the genetic test cost c_0 such that testing everyone in Gen Pop is always best. 
# THis key c_0 cost is done for various various dollars per LYG (lambda)
# Need to hardcode the indices representing each dollars per LYG threshold
##
fixeddollarsperLYG <- c(0,2e4,5e4,1e5) # dollars per LYG thresholds
e_GenPop <- (78-73)*fixeddollarsperLYG-(C_GenPop[3]-C_GenPop[4]) # Net Effectivenesses of early intervention
key <- c(138,85,54,34)
crit_c0 <- (C_GenPop[1]+e_GenPop*probDplusMminus[key]) / probDminusMminus[key]
crit_c0

##
# Various Plots for General Population
##

# Ignore the below that is commented out

# Note that a plot versus optimal risk threshold requires only NB gain with screening test cost and net effectiveness of early intervention
# It does not require the cost of the definitive test (in this case, the cost of BRCA genotyping).
# This could be an advantage if this cost is hard to know - like in lung cancer screening
# However, we never know what are acceptable risk thresholds, but we do know what are acceptable dollars per LYG.
# But to plot versus dollars per LYG, I have to calculate that from the optimal risk threshold, and this
# requires definitive test cost.  So that is the value of knowing the definitive test cost.

# # Plot 2 figures:  INBgain and usual NBgain (i.e. assume costless screening test) vs. Risk-Threshold
# par(mfrow=c(1,2))
# plot(brcaGenPop$Threshold,allGenPop$newINBgain,xlim=c(0.0005,0.1),type="l",log="x",
#      xlab="Optimal Risk Threshold",ylab="Incremental Net Benefit (INB) gain",main="General Population")
# abline(h=0)
# #abline(v=perLYG[optimal],lty=3)
# # Define the usual NBgain for a costless screening test
# brcaGenPop$NBgaincostless <- brcaGenPop$NetBenefit-pmax(0,brcaGenPop$NBalwaysact)
# plot(brcaGenPop$Threshold,brcaGenPop$NBgaincostless,log="x",xlim=c(0.0005,0.1),type="l",ylim=c(-0.002,0.001),lty=2,
#      xlab="Optimal Risk Threshold",ylab="Net Benefit gain",main="General Population")
# abline(h=0)
# # Define the NBgain for a costing screening test
# brcaGenPop$NBgainwithcost <- brcaGenPop$NBgaincostless-C_GenPop[1]/allGenPop$e
# lines(brcaGenPop$Threshold,brcaGenPop$NBgainwithcost,lty=1)
# legend(3.5e-04,-0.0017,lty=c(1,2),bty="n",
#        legend=c("With screening test cost","Ignoring screening test cost"),cex=0.8)


# # Plot NBgain for a costless screening test vs dollars per LYG
# brcaGenPop$NBgaincostless <- brcaGenPop$NetBenefit-pmax(0,brcaGenPop$NBalwaysact)
# plot(allGenPop$perLYG,brcaGenPop$NBgaincostless,xlim=c(0,2e5),type="l",ylim=c(-0.002,0.001),lty=2,
#      xlab="Dollars per life-year gained",ylab="Net Benefit gain",main="General Population")
# abline(h=0)
# # Define the NBgain for a costing screening test
# brcaGenPop$NBgainwithcost <- brcaGenPop$NBgaincostless-C_GenPop[1]/allGenPop$e
# lines(allGenPop$perLYG,brcaGenPop$NBgainwithcost)
# legend(1.5e4,-0.0015,lty=c(1,2),bty="n",
#        legend=c("With screening test cost","Ignoring screening test cost"))


##
# Current Figure 4:  NB gain from DCA vs risk threshold and INB gain vs. dollars per LYG lambda
##
par(mfrow=c(1,2),pty="s")

# Define the usual NBgain for a costless screening test
brcaGenPop$NBgaincostless <- brcaGenPop$NetBenefit-pmax(0,brcaGenPop$NBalwaysact)
plot(brcaGenPop$Threshold,brcaGenPop$NBgaincostless,log="x",xlim=c(0.0005,0.1),type="l",ylim=c(-0.002,0.001),lty=2,
     xlab="Optimal Risk Threshold (R)",ylab="Net Benefit gain from DCA",main="NB gain from DCA")
abline(h=0)
# Define the NBgain for a costing screening test
brcaGenPop$NBgainwithcost <- brcaGenPop$NBgaincostless-C_GenPop[1]/allGenPop$e
lines(brcaGenPop$Threshold,brcaGenPop$NBgainwithcost,lty=1)
legend(5e-03,0.0011,lty=c(1,2),bty="n",
       legend=c("With screening test cost","Ignoring screening test cost"),cex=0.8)
abline(v=brcaGenPop$BRCAprev[2],lty=3)
text(5e-04,-0.002,"If NB_Gain<0,\nthen test none",cex=0.9,pos=4)
text(1.2e-02,-0.002,"If NB_Gain<0,\nthen test all",cex=0.9,pos=2)


# Plot INBgain vs dollars per LYG
plot(allGenPop$perLYG,allGenPop$newINBgain,xlim=c(0,5e5),type="l",
     xlab="Dollars per life-year gained (lambda)",ylab="Incremental Net Benefit (INB) gain",main="INB gain")
abline(h=0)
abline(v=allGenPop$perLYG[27],lty=3) # the dollars per LYG at prevalence: -29609
text(-1e4,-600,"If NB_Gain<0,\nthen test all",cex=0.9,pos=4)
text(3e05,-600,"If NB_Gain<0,\nthen test none",cex=0.9,pos=2)

#abline(v=perLYG[optimal],lty=3)

# If you can't figure this out, delete this
# plot(allGenPop$perLYG,allGenPop$newINB,xlim=c(-1e5,5e5),ylim=c(-1000,1000),type="l",
#      xlab="Dollars per life-year gained",ylab="Incremental Net Benefit (INB) gain",main="INB gain vs. dollars per life-year gained")
# abline(h=0)
# lines(allGenPop$perLYG,allGenPop$newINB-allGenPop$newINBall,lty=3)


# Save this figure
filename <- "testcostsGenPop.pdf" ; dev.copy2pdf(file=filename)
embedFonts(file=filename,outfile=filename) #Needed for TexStudio to read fonts (looks fine in pdf file w/o needed to embed)


##
# Create Table showing key stats (PPV, cNPV, Sens, Spec, AUC, etc) for various willingness to pay lambda
# Choose the willingness to pay lambdas to be the intervals that supports screening, based on 3 criteria:
# 1. INB gain=0
# 2. NB gain=0 (assume costless screening test)
# 3. NB gain=0 (assume screening test cost)
# Note that criteria 1 and 3 should be equivalent, so this is a test to be sure.
##

# Find the key thresholds where INBgain=0
left <- which.max(allGenPop$newINBgain>=0) # Find first point where INBgain switches from neg to pos
right<- max(which(allGenPop$newINBgain>=0)) # Find second point where INBgain switches from pos to neg
highest <- which(allGenPop$newINBgain==max(allGenPop$newINBgain,na.rm=TRUE))
key <- c(left,highest,right) # The 3 indices indicating where INB>0 and the optimal one in the middle
info <- matrix(c(allGenPop$perLYG[key],brcaGenPop$Threshold[key],allGenPop$newINBgain[key],brcaGenPop$OverallPtest[key],
                 brcaGenPop$PPV[key],brcaGenPop$cNPV[key],brcaGenPop$OverallTP[key],1-brcaGenPop$OverallFP[key]),nrow=3)
info <- data.frame(info)
colnames(info) <- c("Dollars/LYG","Threshold","INBgain","Positivity","PPV","cNPV","Sens","Spec")
rownames(info) <- c("Left INBgain>0","Best INBgain", "Right INBgain>0")
info$AUC <- (info$Sens+info$Spec)/2
info.INBgain <- info

# Find the key thresholds where NBgain=0 (costless screening test)
left <- 10+which.max(brcaGenPop$NBgaincostless[-(1:10)]>=0) # Find first point where switches from neg to pos, ignoring crap at begining
right<- max(which(brcaGenPop$NBgaincostless[1:1000]>=0)) # Find second point where INBgain switches from pos to neg, ignoring crap at end
highest <- which(brcaGenPop$NBgaincostless==max(brcaGenPop$NBgaincostless,na.rm=TRUE))
key <- c(left,highest,right) # The 3 indices indicating where INB>0 and the optimal one in the middle
info <- matrix(c(allGenPop$perLYG[key],brcaGenPop$Threshold[key],allGenPop$newINBgain[key],brcaGenPop$OverallPtest[key],
                 brcaGenPop$PPV[key],brcaGenPop$cNPV[key],brcaGenPop$OverallTP[key],1-brcaGenPop$OverallFP[key]),nrow=3)
info <- data.frame(info)
colnames(info) <- c("Dollars/LYG","Threshold","INBgain","Positivity","PPV","cNPV","Sens","Spec")
rownames(info) <- c("Left NBgaincostless>0","Best NBgaincostless", "Right NBgaincostless>0")
info$AUC <- (info$Sens+info$Spec)/2
info.NBgaincostless <- info

# Find the key thresholds where NBgain=0 (with screening test cost)
left <- which.max(brcaGenPop$NBgainwithcost>=0) # Find first point where switches from neg to pos
right<- max(which(brcaGenPop$NBgainwithcost>=0)) # Find second point where switches from pos to neg
highest <- which(brcaGenPop$NBgainwithcost==max(brcaGenPop$NBgainwithcost,na.rm=TRUE))
key <- c(left,highest,right) # The 3 indices indicating where INB>0 and the optimal one in the middle
info <- matrix(c(allGenPop$perLYG[key],brcaGenPop$Threshold[key],allGenPop$newINBgain[key],brcaGenPop$OverallPtest[key],
                 brcaGenPop$PPV[key],brcaGenPop$cNPV[key],brcaGenPop$OverallTP[key],1-brcaGenPop$OverallFP[key]),nrow=3)
info <- data.frame(info)
colnames(info) <- c("Dollars/LYG","Threshold","INBgain","Positivity","PPV","cNPV","Sens","Spec")
rownames(info) <- c("Left NBgainwithcost>0","Best NBgainwithcost", "Right NBgainwithcost>0")
info$AUC <- (info$Sens+info$Spec)/2
info.NBgaincost <- info

info <- rbind(info.INBgain,info.NBgaincostless,info.NBgaincost)
info


##
# Figure 5: Table showing key stats (PPV, cNPV, Sens, Spec, AUC, etc) for various willingness to pay lambda
# lambda are: $20k, $50k, $100k, $137k (optimal INB), $200k, $300k (top end supporting screening)
# Hardcode the indices
##

key <- c(85,54,34,27,20,14)
info <- matrix(c(allGenPop$perLYG[key],brcaGenPop$Threshold[key],allGenPop$newINBgain[key],brcaGenPop$OverallPtest[key],
                 brcaGenPop$PPV[key],brcaGenPop$cNPV[key],brcaGenPop$OverallTP[key],1-brcaGenPop$OverallFP[key]),
               nrow=6)
info <- data.frame(info)
colnames(info) <- c("Dollars/LYG","Threshold","INBgain","Positivity","PPV","cNPV","Sens","Spec")
rownames(info) <- c("Left INBgain>0","Best INBgain", "Right INBgain>0")
info$AUC <- (info$Sens+info$Spec)/2
info



##
# IGNORE EVERYTHING BELOW THIS
##


##
# Comparing two tests for cost-effectiveness
##

# calculate slope of concentration curve (Sens1-Sens2)/(p1-p2), check if greater than critical value R/pi
slope <- (0.9-0.8)/(0.6-0.44)
critical <- R_GenPop/brcaGenPop$BRCAprev[2]
ifelse(slope>critical,"test 1 is better","test 2 is better")

# Calculate minimum dollars per LYG for test 1 to be better
d <- ((C_GenPop[2]/(slope*brcaGenPop$BRCAprev[2]) - C_GenPop[2]) + (C_GenPop[3]-C_GenPop[4])) / ((E[2]-E[3])/dollarsperLYG)
d

# # Trial and error for dollar cost to be less than sensitivity gain of 0.1 for a reduction in positivity of 0.4-0.44=-0.04
# critical*((0.4-0.44)+(204-100)/1000)

# Compare if extra cost of identifying distant relatives is worth it: sens is equal, but positivity decreases 0.44-0.35=0.09
critical <- (R_GenPop/brcaGenPop$BRCAprev[2])*(0.35-0.44+(200-100)/2200)
ifelse(critical<0,"test 1 is better","test 2 is better")
# Trial and error for cost of extra cost screening test to be worth it: $298
(R_GenPop/brcaGenPop$BRCAprev[2])*(0.35-0.44+(298-100)/2200)

# Compare if reduced cost of NCCN is worth it: sens is equal, but positivity increases 0.60-0.68=-0.08
critical <- (R_GenPop/brcaGenPop$BRCAprev[2])*(0.68-0.60+(50-100)/2200)
ifelse(critical<0,"test 1 is better","test 2 is better")
# Trial and error for dollar cost to be less than sensitivity gain of 0 for a reduction in positivity of 0.60-0.68=-0.08
(R_GenPop/brcaGenPop$BRCAprev[2])*(0.68-0.60+(0-100)/2200)



######
## OLD WORK
# Plot INB and INBall versus threshold
# But the problem with this is that it fixes all costs and effectiveness, so there is only 1 threshld to use: the optimal.
# We fix optimal threshold, but then M+ is defined using each possible threshold, which is inconsistent.
#####
par(mfrow=c(1,2))
# First plot for AJs
plot(allAJ$Threshold,allAJ$INBgain,type="l",col="black",log="x",xlab="Threshold",ylab="INBgain",main="Ashkenazi Jews")#,xlim=c(0.002,0.95))#,ylim=c(0,0.01))
abline(h=0) # thresholds for INBgain>0 are permissible for screening
# where is INB maximized at: Risk threshold that maximizes INB 
abline(v=R_AJ,lty=3)
text(R_AJ+.000012,min(allAJ$INBgain,na.rm=TRUE),"max\n INB",cex=0.6,adj=c(0,0))
abline(v=brca$BRCAprev[2],lty=3)
text(brca$BRCAprev[2]-.000012,min(allAJ$INBgain,na.rm=TRUE),"max\n AUC",cex=0.6,adj=c(1,0))
text(1,0.5*(min(allAJ$INBgain,na.rm=TRUE)+max(allAJ$INBgain,na.rm=TRUE)),
     ifelse(R_AJ<brca$BRCAprev[2],"test everyone better\n than test no one","test everyone worse\n than test no one")
     ,cex=0.6,adj=c(1,0))
# Second plot for Gen Pop
plot(allGenPop$Threshold,allGenPop$INBgain,type="l",col="black",log="x",xlab="Threshold",ylab="INBgain",main="General Population")#,xlim=c(0.002,0.95))#,ylim=c(0,0.01))
abline(h=0) # thresholds for INBgain>0 are permissible for screening
# where is INB maximized at: Risk threshold that maximizes INB 
abline(v=R_GenPop,lty=3)
text(R_GenPop+.000012,min(allGenPop$INBgain,na.rm=TRUE),"max\n INB",cex=0.6,adj=c(0,0))
abline(v=brcaGenPop$BRCAprev[2],lty=3)
text(brcaGenPop$BRCAprev[2]-.000012,min(allGenPop$INBgain,na.rm=TRUE),"max\n AUC",cex=0.6,adj=c(1,0))
text(1,0.5*(min(allGenPop$INBgain,na.rm=TRUE)+max(allGenPop$INBgain,na.rm=TRUE)),
     ifelse(R_GenPop<brcaGenPop$BRCAprev[2],"test everyone better\n than test no one","test everyone worse\n than test no one")
     ,cex=0.6,adj=c(1,0))



# OLD WORK
# # Numerically calculate the bounds on dollars per life-year gained so that screening at optimal threshold is best decision
# # This means calculate INBgain as dollarsperLYG increases
######
# perLYG <- seq(1e3,2e5,1e3)
# Es <- (E[2]-E[3])/dollarsperLYG * perLYG # Each possible E[2]-E[3] (e_1-e_2 in the paper) as perLYG value increases
# es <- Es-(C_GenPop[3]-C_GenPop[4]) # Each possible e as perLYG value increases
# Rs <- C_GenPop[2]/(es+C_GenPop[2]) # optimal thresholds as perLYG value increases
# # Figure out the index pointing to each threshold R to get the 2x2 table stats in the whole brca data.frame
# count<-1 ; indices<-NA
# for (R in Rs) {
#   indices[count] <- which(brcaGenPop$Threshold==signif(R,2))
#   count <- count+1
# }
# probDplus <- brcaGenPop$BRCAprev[2]
# probTPs <- brcaGenPop$OverallTP[indices]*probDplus
# probPs <- brcaGenPop$OverallPtest[indices]
# INBs <- (es+C_GenPop[2])*probTPs - C_GenPop[1] - C_GenPop[2]*probPs
# INBalls <- (es+C_GenPop[2])*probDplus-C_GenPop[2]
# INBgains <- INBs - ifelse(INBalls<0,0,INBalls)
# screennone <- (Rs>probDplus) # 1=>screen none is better than test all, 0=> screen none is worse than test all
# optimal <- which.max(INBgains)
# optindex <- indices[optimal]
# c(Rs[optimal],INBgains[optimal],probPs[optimal],brcaGenPop$PPV[optindex],brcaGenPop$cNPV[optindex],
#   brcaGenPop$OverallTP[optindex],1-brcaGenPop$OverallFP[optindex])
# 
# # Find the key thresholds where INB=0
# left <- which.max(INBgains>=0) # Find first point where INBgain switches from neg to pos
# right<- max(which(INBgains>=0)) # Find second point where INBgain switches from pos to neg
# prev <- max(which(Rs>=probDplus))
# key <- c(left,optimal,right) # The 3 indices indicating where INB>0 and the optimal one in the middle
# keyind <- indices[key]
# info <- matrix(c(perLYG[key],Rs[key],INBgains[key],probPs[key],brcaGenPop$PPV[keyind],brcaGenPop$cNPV[keyind],
#            brcaGenPop$OverallTP[keyind],1-brcaGenPop$OverallFP[keyind]),nrow=3)
# info <- data.frame(info)
# colnames(info) <- c("Dollars/LYG","Threshold","INBgain","Positivity","PPV","cNPV","Sens","Spec")
# rownames(info) <- c("Left INBgain>0","Best INBgain", "Right INBgain>0")
# info$AUC <- (info$Sens+info$Spec)/2
# # Plot
# plot(perLYG,INBgains,type="l")
# abline(h=0)
# abline(v=perLYG[optimal],lty=3)





## OLD PLOTS  
# Make 4 plots, 2 for Ashkenazi-Jews (first row), 2 for General Population (second row)
# In each row: Plot deltas of cost and effectiveness vs. threshold, then plot MRS and INB vs. threshold
###########
# par(mfrow=c(2,2))
# par(mar=c(3,2,2,4)+0.1,mgp=c(2.2,1,0))
# 
# ##
# # First row of plots is for Ashkenazi-Jews
# ##
# 
# ##
# # Put delta Cost and delta Effectiveness on one plot
# ##
# plot(allAJ$Threshold,allAJ$delE,type="l",col="red",log="x",yaxt="n",xaxt="n",xlab="Threshold",ylab="",main="Ashkenazi-Jews")#,xlim=c(0.002,0.95))#,ylim=c(0,0.018))
# axis(4,col="red")
# mtext("delta Effectiveness",side=4,line=2,col="red")
# abline(v=brca$BRCAprev[10],lty=2)
# # text(brca$BRCAprev[10]-.0012,0.0002+.00003,"BRCA1/2 mutation",cex=0.6,adj=c(1,0))
# # text(brca$BRCAprev[10]-.0012,0.0002-.00003,"prevalence",cex=0.6,adj=c(1,1))
# # delta Cost 
# par(new=TRUE)
# plot(allAJ$Threshold,allAJ$delC,type="l",col="blue",log="x",yaxt="n",xlab="",ylab="")#,xlim=c(0.002,0.95))#,ylim=c(0,0.01))
# axis(2,col="blue")
# mtext("delta Cost",side=2,line=2,col="blue")
# legend("topright",col=rev(c("red","blue")),lty=1,legend=rev(c("delta Effective","delta Cost")),bty="n",cex=0.8)
# 
# ##
# # Put MRS and INB on one plot
# ##
# plot(allAJ$Threshold,allAJ$MRS,type="l",col="red",log="x",yaxt="n",xaxt="n",xlab="Threshold",ylab="",main="Ashkenazi-Jews")#,xlim=c(0.002,0.95),ylim=c(0,0.018))
# axis(4,col="red")
# mtext("Mean Risk Stratification (MRS)",side=4,line=2,col="red")
# abline(v=brca$BRCAprev[10],lty=2)
# # text(brca$BRCAprev[10]-.0012,0.0002+.00003,"BRCA1/2 mutation",cex=0.6,adj=c(1,0))
# # text(brca$BRCAprev[10]-.0012,0.0002-.00003,"prevalence",cex=0.6,adj=c(1,1))
# # where is INB maximized at
# abline(v=INBmaxat,lty=3)
# text(INBmaxat+.000012,0.0002+.00003,"max INB",cex=0.6,adj=c(0,0))
# # INB 
# par(new=TRUE)
# plot(allAJ$Threshold,allAJ$INB,type="l",col="blue",log="x",yaxt="n",xlab="",ylab="")#,xlim=c(0.002,0.95))#,ylim=c(0,0.01))
# axis(2,col="blue")
# mtext("INB",side=2,line=2,col="blue")
# legend("topright",col=rev(c("red","blue")),lty=1,legend=rev(c("MRS","INB")),bty="n",cex=0.8)
# 
# 
# 
# ##
# # Second row of plots is for General-Population
# ##
# 
# ##
# # Put delta Cost and delta Effectiveness on one plot
# ##
# plot(allGenPop$Threshold,allGenPop$delE,type="l",col="red",log="x",yaxt="n",xaxt="n",xlab="Threshold",ylab="",main="General Population")#,xlim=c(0.0002,0.95))#,ylim=c(0,0.018))
# axis(4,col="red")
# mtext("delta Effectiveness",side=4,line=2,col="red")
# abline(v=brcaGenPop$BRCAprev[10],lty=2)
# # text(brcaGenPop$BRCAprev[10]+.00012,0.00002+.000003,"BRCA1/2 mutation",cex=0.6,adj=c(0,0))
# # text(brcaGenPop$BRCAprev[10]+.00012,0.00002-.000003,"prevalence",cex=0.6,adj=c(0,1))
# # delta Cost 
# par(new=TRUE)
# plot(allGenPop$Threshold,allGenPop$delC,type="l",col="blue",log="x",yaxt="n",xlab="",ylab="")#,xlim=c(0.0002,0.95))#,ylim=c(0,0.01))
# axis(2,col="blue")
# mtext("delta Cost",side=2,line=2,col="blue")
# legend("topright",col=rev(c("red","blue")),lty=1,legend=rev(c("delta Effective","delta Cost")),bty="n",cex=0.8)
# 
# ##
# # Put MRS and INB on one plot
# ##
# plot(allGenPop$Threshold,allGenPop$MRS,type="l",col="red",log="x",yaxt="n",xaxt="n",xlab="Threshold",ylab="",main="General Population")#,xlim=c(0.0002,0.95))#,ylim=c(0,0.018))
# axis(4,col="red")
# mtext("Mean Risk Stratification (MRS)",side=4,line=2,col="red")
# abline(v=brcaGenPop$BRCAprev[10],lty=2)
# # text(brcaGenPop$BRCAprev[10]+.00012,0.00002+.000003,"BRCA1/2 mutation",cex=0.6,adj=c(0,0))
# # text(brcaGenPop$BRCAprev[10]+.00012,0.00002-.000003,"prevalence",cex=0.6,adj=c(0,1))
# # where is INB maximized at
# abline(v=INBmaxat,lty=3)
# text(INBmaxat+.000012,0.0002+.00003,"max INB",cex=0.6,adj=c(0,0))
# # INB
# par(new=TRUE)
# plot(allGenPop$Threshold,allGenPop$INB,type="l",col="blue",log="x",yaxt="n",xlab="",ylab="")#,xlim=c(0.0002,0.95))#,ylim=c(0,0.01))
# axis(2,col="blue")
# mtext("INB",side=2,line=2,col="blue")
# legend("topright",col=rev(c("red","blue")),lty=1,legend=rev(c("MRS","INB")),bty="n",cex=0.8)

  