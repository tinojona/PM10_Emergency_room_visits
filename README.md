# PM10 and emergency room visits


This repository presents the analysis on the association between emergency room visits and exposure to PM10 (particulate matter with an aerodynamic diameter of less than 10 micrometers) in Basel, Switzerland. The final report (FINAL_REPORT.Rmd .html) summarizes the findings of the analysis, sets them into perspective with the latest literature and discusses them. 

<br> 

## Main results

We found PM10 exposure to increase the risk of an emergency room visit, which however decreased when we accounted for temperature confounding.

<p align="center">
  <img src="plots/model1_model2_all.png" alt="Figure 1" width="600"/>
</p>

*Cummulative relative risk for all-cause emergency room visits from Model 1 (erv ~ PM10) and Model 2 (erv ~ PM10 + temperature).*


<br>

## Folders

### The analysis folder

The `analysis` folder contains the different analysis conducted to meet the research objectives, starting off with descriptive statistics (01) and then the effect of PM10 (02) and the modification of this effect from heat (03). Additionally we conducted sensitivity analysis (04) and created a map of Basel (05), where our hospital and environmental data sources are indicated. 

### The data folder

The `data`folder contains the raw and processed environmental and hospitalization data and the processing routine. 


### The plots folder

This folder contains all the output that was generated in the analysis, which were then used in the final report. 
