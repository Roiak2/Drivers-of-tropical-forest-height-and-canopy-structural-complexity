# Water, Wind, & Time - Drivers of tropical forest height and canopy structural complexity

## Overview
This repository contains code and data for chapter 2 of my dissertation, which is to be published in *Ecology* soon under Ankori-Karlinsky et al., 2022. 

This dissertation focuses on the relationship between tree architecture, canopy structure, and chronic wind-exposure.

When referring to canopy structure and structural complexity, there are two salient metrics related to ecosystem function:

<img src="/Figures/Canopy metrics.png" height="500">

<p align="center">Drawing by Nina Berinstein</p> 

These metrics have been connected to productivity [[1]](#1) and have been posited to be reduced by chronic wind-exposure [[2]](#2).
 
Therefore, this chapter explores the following question:
**Does chronic wind-exposure reduce stand-level canopy height and structural complexity?**

## Methods

To address this question, we randomly sampled ~20,000 30 m-radius forested sites across Puerto Rico stratified by forest age, and quantified height and rugosity with 2016 airborne LiDAR from the USGS [[3]](#3).  

We then ran random forest models to predict each metric based on chronic wind-exposure, forest age, mean annual precipitation, topography, soil properties, and exposure to previous hurricanes (most recent in 1998).


## References
<a id="1">[1]</a> 
Gough, C.M., Atkins, J.W., Fahey, R.T. & Hardiman, B.S. (2019). High rates of primary production in structurally complex forests. Ecology, 100, e02864.
<a id="2">[2]</a> 
Eloy, C., Fournier, M., Lacointe, A. & Moulia, B. (2017). Wind loads and competition for light sculpt trees into self-similar structures. Nature Communications, 8, 1014.  
<a id="3">[3]</a> 
Carswell Jr., W.J. (2016). The 3D Elevation Program: summary for Puerto Rico (USGS Numbered Series No. 2015–3088). The 3D Elevation Program: summary for Puerto Rico, Fact Sheet. U.S. Geological Survey, Reston, VA.