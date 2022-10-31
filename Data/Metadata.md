# Drivers of tropical forest height and canopy structural complexity across heterogeneous landscapes:
---

This dataset was used for analyses and figures in Ankori-Karlinsky et al., 2022 in *Ecology*. 

The dataset contains 20,660 30 m-radius forested sites in Puerto Rico stratified by forest age class (see Martinuzzi et al., 2020 for details). 

For each site, we quantified canopy height and rugosity using 2016 airborne USGS LiDAR data.

The dataset was used to predict both metrics based on forest age, exposure to chronic winds, precipitation, topography, soil propoerties, and exposure to previous hurricanes.

## Metadata 

Each row is a 30 m-radius site.
Below is information on each column:  
    - *x* - longitude in NAD83 (EPSG 4269)
    - *y* - latitude in NAD83  
    - *Forest_Age* - Forest age class (1-5):    
            - 5-16 years   
            - 17-25 years   
            - 26-39 years   
            - 40-65 years   
            - 66+ years        
    - *CHM* - Mean height (m) from canopy height model   
    - *Max_Height* - Maximum height (m)   
    - *Rugosity* - Stdev in height from 15m moving windows   
    - *MAP* - Mean annual precipitation (mm/year) from PRISM   
    - *Wind_Binary* - Protected (0) vs. exposed (1)    
    - *Elev_Mean* - Mean elevation (m) from USGS DEM   
    - *Slope* - Mean slope (in degrees)   
    - *Soil* - Soil type (e.g. volcanic vs. limestone) from USGS    
    - *Hugo/Georges exposure* - Mean exposure from 0-1 to hurricanes Hugo (1989) and Georges (1998) from Boose et al., (2004)    
    - *AWS_mean* - Available water storage in soils from gSSURGO   


## Sharing/access Information

Please cite Ankori-Karlinsky et al., 2022 *Ecology* if you use these data.