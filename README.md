# [Movement-assisted localization from acoustic telemetry data](https://doi.org/10.1186/s40462-020-00199-6)

### Hostetter, N. J., and J. A. Royle

### Movement Ecology

### Code/Data DOI: [doi.org/10.1186/s40462-020-00199-6](https://doi.org/10.1186/s40462-020-00199-6)

### Please contact the first author for questions about the code or data requests: Nathan Hostetter (njhostet@uw.edu)
__________________________________________________________________________________________________________________________________________

## Abstract:
Background: Acoustic telemetry technologies are being increasingly deployed to study a variety of aquatic taxa including fishes, reptiles, and marine mammals. Large cooperative telemetry networks produce vast quantities of data useful in the study of movement, resource selection and species distribution. Efficient use of acoustic telemetry data requires estimation of acoustic source locations from detections at receivers (i.e., “localization”). Multiple processes provide information for localization estimation including detection/non-detection data at receivers, information on signal rate, and an underlying movement model describing how individuals move and utilize space. Frequently, however, localization methods only integrate a subset of these processes and do not utilize the full spatial encounter history information available from receiver arrays.
Methods: In this paper we draw analogies between the challenges of acoustic telemetry localization and newly developed methods of spatial capture-recapture (SCR). We develop a framework for localization that integrates explicit sub-models for movement, signal (or cue) rate, and detection probability, based on acoustic telemetry spatial encounter history data. This method, which we call movement-assisted localization, makes efficient use of the full encounter history data available from acoustic receiver arrays, provides localizations with fewer than three detections, and even allows for predictions to be made of the position of an individual when it was not detected at all. We demonstrate these concepts by developing generalizable Bayesian formulations of the SCR movement-assisted localization model to address study-specific challenges common in acoustic telemetry studies.
Results: Simulation studies show that movement-assisted localization models improve point-wise RMSE of localization estimates by > 50% and greatly increased the precision of estimated trajectories compared to localization using only the detection history of a given signal. Additionally, integrating a signal rate sub-model reduced biases in the estimation of movement, signal rate, and detection parameters observed in independent localization models.
Conclusions: Movement-assisted localization provides a flexible framework to maximize the use of acoustic telemetry data. Conceptualizing localization within an SCR framework allows extensions to a variety of data collection protocols, improves the efficiency of studies interested in movement, resource selection, and space-use, and provides a unifying framework for modeling acoustic data.

## Code 
1. [simulation_analysis](./simulation_analysis/): This folder contains the code to simulate and analyze data. It also generates the JAGS file for the Bayesian implementation.

## Data
All data are reproducible from the simulation script. 
