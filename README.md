# Zoopl_carcass
Data and R code for publication ""Sinking dead" - How zooplankton contribute to Particulate Organic Carbon flux in the subantarctic Southern Ocean". Currently under review at the journal Limnology and Oceanography. 

The repository contains data files and the R script to perform the calculations and produce the figures. See the script (Carcass_script.R) for further details on data accessibility. 

Data files include: 
- **CTD** data on temperature, salinity, PAR and other parameters: CTD_stationPS1.csv and CTD_stationPS2.csv.
- Remote-sensed data on **ocean colour** as a proxy for chlorophyll a: MODISA_chl.csv
- Measured **chlorophyll a** in the upper water column (0-200 m): Chlorophyll_measured.csv
- Sinking **particle flux** data from the sediment traps at SOTS (Southern Ocean Time Series) site: IMOS_ABOS-SOTS_KF_20180322_SAZ47_FV01_SAZ47-20-2018_PARFLUX-Mark78H-21_END-20180322_C-20200416.nc
- **POC/PON** content of the upper water column (0-200 m), measured by water fitration of CTD sampler water: POM_measured
- **Prosome lengt** of the collected copepod *Neocalanus tonsus*: Prosome_length_sample.csv
- The measurements of **copepod prosome length, dry weight and carbon content**: Copepod_carbon.csv
- The measured **sinking velocity** by copepod carcasses: Sinking_speed.csv
- **Oxygen measurements** every 4 hours over a time of 24 hours in two batches (replicates): batch_1 & batch_2

For questions or comments, please pull a request or contact me under Svenja.Halfter@utas.edu.au. 
