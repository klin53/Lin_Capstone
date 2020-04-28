---
title: "Lin_Kun_Capstone"
author: "Kun Lin"
date: "4/27/2020"
output: html_document
---

```{r setup, include=FALSE}
library(ez)
library(viridis)
library(tidyverse)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.


##Background
**1)Myotonic dystrophy type 1 (DM1) is a microsatellite repeat expansion disease, characterized by expanded CTG repeats (hundreds to thousands of repeats long in affected patients), these CTG repeats form intranuclear foci that sequester a family of RNA binding proteins, MBNL, functionally inactivating MBNL in the cells. There are many Central Nervous System (CNS) symptoms of the disease including developmental deplay, hypersomnia, although the underlying mechanisms are not clear. Prelminary data show that in neurons expressing expanded CTG repeats (hSyn-480CTG), Snap25 mRNA transcripts are not properly localzied to distal neurites. SNAP25 is a protein important for proper synaptic functions, regulating pre-synaptic release and post-synaptic receptor density.**
##Knowledge gap

**2) It is not known whether expressing 480CTG repeats in neurons reduces the level of NMDA receptors in the post-synaptic terminal**

##Prediction
**3)If expressing 480 CTG repeats sufficiently reduces the local Snap25 mRNA transcript levels (therefore SNAP25 protein concentration at the post0synapse), then NMDA receptor densities will be lower in cells expressing 480CTG repeats**

##Variables
**4) The dependent variable will be the fluorescence intensity of NMDA receptors using fluorescence imaging (the dependent variable is a continuous measured variable). The predictor variable will be the virus that's expressed in these neurons (sorted, 3 levels: No Viral Transduction; Transduction control (expresses GFP); 480CTG, expresses 480 CTG repeats) **

##Hypothesis
**5) Null hypothesis: Expression of 480 CTG repeats does not change the level of NMDA receptors in the post-synaptic terminal. Alternate hypothesis: Expression of 480CTG repeats reduces the level of NMDA receptors in the post-synaptic terminal.**

##Stats test chocie 
**6)I will use one-way completely random ANOVA. There are 3 predictor levels: no transduction, No Viral Transduction; Transduction control (expresses GFP); 480CTG, expresses 480 CTG repeats), therefore, t-tests are not appropriate. Since there is only one predictor variable, it cannot be a two-way ANOVA. It is completely random because the unit of quantification is individual cells, and measurement is taken once only. The measurements of individual cells are not intrinsically related, and the heterogeneity of primary cortical neruons is sufficient that so they can be considered independent when few cells (5) are sampled per experiemnt.**

##experimental flow and defense of individual replicates for purpose of analysis
**7) The cultured cells will be from single-cell suspension of cortical neurons from multiple embryos. The neurons will be plated on glass coverslips in 6 well plates. Each well will be randomized to location to receive a virus treatment: No Viral Transduction; Transduction control (expresses GFP); 480CTG, expresses 480 CTG repeats). All coverslips from a single experiment will be PFA fixed and stained together, then imaged using identical microscope and processing parameters before analysis. The units of analysis (independent replicates) will be individual neurons on the coverslip, this analysis can only be done with high-resolution images given that only the dendritic spine will be analyzed, the fluorescence intensity of NMDA receptors will be normalized to the surface area of the dendritic spine.Therefore individual cells are used as units of analysis, as opposed to whole coverslip as it's impractical to image the entier coverslip at high resolution (60X). 5 neurons will be randomly sampled from each group from each experiment, and the experiment will be repeated 5 times to get 25 neurons per condition. Although 5 neurons will be from the same culture from each group, there exists inherent heterogeneity in cortical neuron cultures given the numerous neuron types that exist in different cortical layers. Addtionally, averaging the measurements from multiple neurons within each culture doesn't yield anything meaningful when the quantification is done at the single cell level. Repeating the same experiment dozens of times becomes prohibitively expensive. The experiment is repeats 5 times to ensure that anyconclusions are not due to anomalies in one experiment. However, the indepdent replicates remain indiviual neurons as averaging the measurements doesn't yield anything meaningful.**

##Data Simulation
**8)**

```{r}
b = 150 #expected basal outcome value
a = .8 #expected fold-to-basal effect of treatment
f = 0.95 #viral transduction control, not expected to change outcoem measurement
sd = 30 #expected standard deviation of Outcome variable
n = 25 # number of independent replicates per group
sims = 100 #number of Monte Carlo simulations to run. 
```

```{r}
##plot data

CRdataMaker <- function(n, b, a, f, sd) { 
  
  
  a1 <- rnorm(n, b, sd) #basal or negative ctrl
  a2 <- rnorm(n, (b*a), sd) #positive control or some other treatment
  a3 <- rnorm(n, (b*f), sd) #treatment effect
    
    Outcome <- c(a1, a2, a3)
    Predictor <- c(rep(c("Not Transduced", "hSyn-GFP", "hSyn-CTG480"), each = n))
    ID <- as.factor(c(1:length(Predictor)))
    df <-data.frame(ID, Predictor, Outcome)
    }

dat <- CRdataMaker(n,b,a,f,sd)

ggplot(dat, aes(Predictor, Outcome))+  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult=1), 
               geom="crossbar", 
               width=0.2, 
               color="red"
               ) + 
  geom_jitter(width=0.1,size = 4, alpha=0.5)+ggtitle("Effect of Expanded CTG Repeats Expression on the Level of NMDA Receptors")+ylab("Fluorescence Intensity Normalized to Surface Area of Dendritic Spine, A.U.")+xlab("Virus Transduction")
```
##Monte Carlo power analysis

```{r}

pval <- replicate(
  sims, {
 
    sample.df <- CRdataMaker(n, b, a, f, sd)
    
    sim.ezaov <- ezANOVA(
            data = sample.df, 
            wid = ID,
            dv = Outcome,
            between = Predictor,
            type = 2
            )
  
  pval <- sim.ezaov$ANOVA[1,5]
    
    }
  )

pwr.pct <- sum(pval<0.05)/sims*100
paste(pwr.pct, sep="", "% power. Change 'n' in your initializer for higher or lower power.")
ggplot(data.frame(pval))+
  geom_histogram(aes(pval), color="#d28e00")+
  labs(x="p-value")
```
**The proposed sample size of 25 neurons/group is sufficiently large to test the hypothesis with 90% power, sufficient to detect the estiamted difference of 25% from baseline. **
