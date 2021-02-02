Data availability
Data associated with this study are available on the Open Science Framework https://doi.org/10.17605/OSF.IO/TJ6C5. Bird distribution and habitat data are available from http://datazone.birdlife.org/species/requestdis; amphibian, mammal and reptiles distribution and habitat data are available from https://www.iucnredlist.org/; Historical and future climate simulations are available from http://doi.org/10.17605/OSF.IO/5XMRJ .

Code availability
Code associated with this study is available on https://github.com/qiaohj/diversity_in_e

Detailed description
System requirements and Runtime

This code was developed, tested and run on on Ubuntu 20.04 (Intel(R) Xeon(R) Silver 4110 CPU @ 2.10GHz, 1.5 TB RAM), using R v4.0.3.

Installation guide

No installation is required. The preprocessing and main scripts are run in the appropriate order detailed below. All the necessary data can be downloaded via http://doi.org/10.17605/OSF.IO/5XMRJ. Bird distribution and habitat data are available from http://datazone.birdlife.org/species/requestdis; amphibian, mammal and reptiles distribution and habitat data are available from https://www.iucnredlist.org/;



Folder structure in OSF (http://doi.org/10.17605/OSF.IO/5XMRJ)

'environmental layers' contains the historical and future environmental layers used in this project. We used R object format to save the layers to make the analysis faster and simpler.
'IUCN_List' contains the species list for different taxonomic groups.
'IUCN_Distribution' contains the range map for all the species in this analysis.
'mask' contains the mask layers which acted as the index for each pixel in the analysis. 
'Niche_Models' contains the occurrences, niche estimation, potential distribution from 2020 to 2100 under different GCM and SSP combinations, and the year-by-year, dispersal-based distribution for each species in the project. You can also generate the data with the scripts in https://github.com/qiaohj/diversity_in_e.

Overview of scripts used in GitHub (https://github.com/qiaohj/diversity_in_e)
Preparing data: folder "Prepare"

*  the scripts in folder 'ENM' were used to prepare the environmental layers, including 1) extract the monthly values from netCDF files (extract_netcdf_monthly.r); 2) get the annual variables based on the monthly values (get_annual_vars.r); 3) reproject the layers to ECK4 projection (annual_vars_2_eck4.r); and 4) fix the variables to fulfill the requirements in this project (fix_env_layers_to_2020.r). The final dataset were stored with the compressed R object format which can be downloaded via http://doi.org/10.17605/OSF.IO/5XMRJ (in 'environmental layers' folder).

*  the scripts in folder 'IUCN_2_Raster' were used to extract the range map and store them into the compressed R objects species by species. The objects can be downloaded via http://doi.org/10.17605/OSF.IO/5XMRJ (in 'IUCN_Distribution' folder).

*  the script 'create_eck4_mask.r' contains the code to generate an index-based mask for the analysis.

Testing the methods used in the project: folder "Methods"

The files in this folder are used as the pre-experiment scripts, to test the feasibility and effectiveness.
*  niche_definition.r: used to test the different methods to estimate the niche, remove the outliers, and so on.
*  dispersal_comparation.r: used to compare the different dispersal abilities. Here we tested 0-5 pixels per year, and moved 0-2 times per year. Finally, we decided to use 0 and 1 as the only dispersal distances in our project. 
*  dispersal_pathway.r: Test the methods to calculate the dispersal pathway.

Making models: folder "Model"

This folder contains the scripts to make the ecological niche models, dispersal models, and calculate the relevant diversity metrics in this project. The scripts should be run in the appropriate order detailed below.

*  models.r: used to generate the ecological niche models for each species in the project. 
*  predict.r: predicts the year-by-year potential distribution for the species from 2021 to 2100.
*  dispersal_1.r: estimate the dispersal-based distribution under the no exposure scenario.
*  dispersal_5.r: estimate the dispersal-based distribution under the 5-year exposure scenario.
*  diversity.r: generate a distribution matrix based on the species distribution maps.
*  diversity_matrix.r: calculate the diversity metrics using the matrix above.

Making Figures: folder "Figures"
The scripts in this folder are used to generate the figures (most of them are unused).
