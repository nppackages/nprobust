########################################################################
## REPLICATION FILES: Calonico, Cattaneo and Farrell (2018)
## nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference
## Last update: 21-Oct-2019
########################################################################
rm(list=ls(all=TRUE))
library(TeachingDemos); library(ggplot2);
library(nprobust)

#######################################################
## SET =TRUE TO GENERATE OUTPUTS
do.output = TRUE
#######################################################

## Install NPPACKAGE Package
if(do.output) txtStart("output/nprobust_1.txt", results = FALSE)
install.packages("nprobust")
if(do.output) txtStop()

## Help NPPACKAGE Package
if(do.output) txtStart("output/nprobust_2.txt", results = FALSE)
help(package="nprobust")
if(do.output) txtStop()

if(do.output) txtStart("output/nprobust_3.txt", results = FALSE)
library(nprobust)
if(do.output) txtStop()

if(do.output) txtStart("output/nprobust_4.txt")
chole = read.csv("nprobust_data.csv")
summary(chole)
if(do.output) txtStop()
attach(chole)

if(do.output) txtStart("output/nprobust_5.txt")
f0_chol1 = kdrobust(chol1, subset=(t==0), neval=7)
summary(f0_chol1)
if(do.output) txtStop()

summary(kdrobust(chol1, subset=(t==0), neval=7, bwselect="MSE-DPI"))
summary(kdrobust(chol1, subset=(t==0), neval=30, bwselect="IMSE-DPI"))

if(do.output) txtStart("output/nprobust_6.txt", results = FALSE)
ev1=seq(250,350,length.out=30)
ev2=seq(0,100,  length.out=30)
f0_chol1 = kdrobust(chol1, subset=(t==0), eval=ev1)
f1_chol1 = kdrobust(chol1, subset=(t==1), eval=ev1)
f0_chol2 = kdrobust(chol2, subset=(t==0), eval=ev1)
f1_chol2 = kdrobust(chol2, subset=(t==1), eval=ev1)
f0_cholf = kdrobust(cholf, subset=(t==0), eval=ev1)
f1_cholf = kdrobust(cholf, subset=(t==1), eval=ev1)
f0_comp  = kdrobust(comp,  subset=(t==0), eval=ev2)
f1_comp  = kdrobust(comp,  subset=(t==1), eval=ev2)
if(do.output) txtStop()

if(do.output) txtStart("output/nprobust_7.txt")
nprobust.plot(f0_chol1,f1_chol1, legendGroups=c("Control Group", "Treatment Group"), 
                                 xlabel="Cholesterol at Baseline 1", ylabel="Density")+theme(legend.position=c(.4,.2))
ggsave("output/kd-chol1.pdf", width=4, height=4)
nprobust.plot(f0_chol2,f1_chol2, xlabel="Cholesterol at Baseline 2",
                                 ylabel="Density")+theme(legend.position="none")
ggsave("output/kd-chol2.pdf", width=4, height=4)
nprobust.plot(f0_cholf,f1_cholf, xlabel="Cholesterol after Treatment",
                                 ylabel="Density")+theme(legend.position="none")
ggsave("output/kd-cholf.pdf", width=4, height=4)
nprobust.plot(f0_comp,f1_comp, xlabel="Treatment Compliance",
                               ylabel="Density")+theme(legend.position="none")
ggsave("output/kd-comp.pdf", width=4, height=4)
if(do.output) txtStop()

t.test(cholf[t==0],cholf[t==1])
t.test(comp[t==0],comp[t==1])

if(do.output) txtStart("output/nprobust_8.txt", results = FALSE)
ev=seq(250,350,length.out=30)
m0_cholf_1 = lprobust(cholf, chol1, subset=(t==0), eval=ev)
m1_cholf_1 = lprobust(cholf, chol1, subset=(t==1), eval=ev)
m0_comp_1  = lprobust(comp,  chol1, subset=(t==0), eval=ev)
m1_comp_1  = lprobust(comp,  chol1, subset=(t==1), eval=ev)
m0_cholf_2 = lprobust(cholf, chol2, subset=(t==0), eval=ev)
m1_cholf_2 = lprobust(cholf, chol2, subset=(t==1), eval=ev)
m0_comp_2  = lprobust(comp,  chol2, subset=(t==0), eval=ev)
m1_comp_2  = lprobust(comp,  chol2, subset=(t==1), eval=ev)
if(do.output) txtStop()

if(do.output) txtStart("output/nprobust_9.txt")
summary(lprobust(cholf, chol1, subset=(t==0), eval=ev[1:7]))
if(do.output) txtStop()

if(do.output) txtStart("output/nprobust_10.txt")
nprobust.plot(m0_cholf_1,m1_cholf_1, legendGroups=c("Control Group", "Treatment Group"),
                                     xlabel="Cholesterol at Baseline 1",
                                     ylabel="Cholesterol after Treatment")+theme(legend.position=c(.3,.8))
ggsave("output/lp-cholf-1.pdf", width=4, height=4)
nprobust.plot(m0_cholf_2,m1_cholf_2, xlabel="Cholesterol at Baseline 2",
                                     ylabel="Cholesterol after Treatment")+theme(legend.position="none")
ggsave("output/lp-cholf-2.pdf", width=4, height=4)
nprobust.plot(m0_comp_1,m1_comp_1, xlabel="Cholesterol at Baseline 1",
                                   ylabel="Treatment Complaince")+theme(legend.position="none")
ggsave("output/lp-comp-1.pdf", width=4, height=4)
nprobust.plot(m0_comp_2,m1_comp_2, xlabel="Cholesterol at Baseline 2",
                                   ylabel="Treatment Complaince")+theme(legend.position="none")
ggsave("output/lp-comp-2.pdf", width=4, height=4)
if(do.output) txtStop()

if(do.output) txtStart("output/nprobust_11.txt")
  summary(lpbwselect(cholf, chol1, subset=(t==0), eval=ev[1:7]))
if(do.output) txtStop()

summary(lpbwselect(cholf, chol1, subset=(t==0), eval = seq(250,350,length.out=30)[1:7]))
summary(lpbwselect(cholf, chol1, subset=(t==0), eval = seq(250,350,length.out=30)[1:7], bwselect="CE-DPI"))

if(do.output) txtStart("output/nprobust_12.txt")
summary(lpbwselect(cholf, chol1, subset=(t==0), eval = ev[1:7], bwselect="ALL"))
if(do.output) txtStop()

#######################################################
## REPLICATION: Figure 1 of Efron and Feldman (1991)
#######################################################
y = 0.25*chol1 + 0.75*chol2 - cholf
x = comp
t = t

m0=lprobust(y,x,subset=(t==0),neval=100)
m1=lprobust(y,x,subset=(t==1),neval=100)
nprobust.plot(m1,m0, xlabel="Compliance", ylabel="Cholesterol Difference")

