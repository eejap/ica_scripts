## FastICA of subsidence InSAR time series

This code executes FastICA on LiCSBAS InSAR time series analysis outputs that focus on regions of subsidence.

Following ICA, inelastic components of subsidence are identified.

All identified components are reconstructed and saved in a netCDF file, along with several other parameters.

This release accompanies work by Payne et al., submitted to JGR Solid Earth.

Below demonstrates an exmample output using 5 components and temporal ICA on ascending LiCS frame 086A_05410_131313, 2016-2022, in north-east Iran.
Independent component (IC) 2 is identified by the skew of the gradients as the inelastic component of subsidence.

![086A_05410_131313_output_demonstration_labels](https://github.com/user-attachments/assets/83542d5b-a270-4e1f-babe-cb563b98f97b)

## Citations
Payne, J. A.; Watson, A.; Maghsoudi, Y.; Lazecký, Milan; Watson, C. S.; Elliot, J. R. Nationwide assessment of subsidence induced hazard in Iran using Sentinel-1 InSAR. *AGU Fall Meetings Abstracts*. **2023**, NS32A-04.

## Acknowledgements
This work has been accomplished during Jessica Payne’s PhD at the University of Leeds, funded by SENSE CDT studentship.

COMET is the UK Natural Environment Research Council's Centre for the Observation and Modelling of Earthquakes, Volcanoes and Tectonics. LiCSAR is developed as part of the NERC large grant, "Looking inside the continents from Space" (NE/K010867/1). LiCSAR contains modified Copernicus Sentinel data [2014-] analysed by the COMET. LiCSAR uses JASMIN, the UK’s collaborative data analysis environment.

The Scientific Colour Maps (Crameri, 2018) is used.

## References
Lazecký, M.; Spaans, K.; González, P.J.; Maghsoudi, Y.; Morishita, Y.; Albino, F.; Elliott, J.; Greenall, N.; Hatton, E.; Hooper, A.; Juncu, D.; McDougall, A.; Walters, R.J.; Watson, C.S.; Weiss, J.R.; Wright, T.J. LiCSAR: An Automatic InSAR Tool for Measuring and Monitoring Tectonic and Volcanic Activity. *Remote Sens*. **2020**, 12, 2430, https://doi.org/10.3390/rs12152430.

Morishita, Y.; Lazecky, M.; Wright, T.J.; Weiss, J.R.; Elliott, J.R.; Hooper, A. LiCSBAS: An Open-Source InSAR Time Series Analysis Package Integrated with the LiCSAR Automated Sentinel-1 InSAR Processor. *Remote Sens*. **2020**, 12, 424, https://doi.org/10.3390/RS12030424 (Remote Sensing 2022 Best Paper Award).

Hyvarinen, A. & Oja, E. Indepedent component analysis: algorithms and applications. *Neural Networks*. **2000**, 13, 4, https://doi.org/10.1016/S0893-6080(00)00026-5

Ebmeier, S. K. Application of independent component analysis to multitemporal InSAR data with volcanic case studies. *JGR: Solid Earth*. **2016**, 121, 12, https://doi.org/10.1002/2016JB013765.
