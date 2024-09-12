# ICA of subsidence InSAR time series

This code executes ICA on LiCSBAS InSAR time series analysis outputs that focus on regions of subsidence.

Following ICA, inelastic components of subsidence are identified.

All identified components are reconstructed and saved in a netCDF file, along with several other parameters.

This release accompanies work by Payne et al., submitted to JGR Solid Earth.

Below demonstrates an exmample output using 5 components and temporal ICA on ascending LiCS frame 086A_05410_131313, 2016-2022, in north-east Iran.
Independent component (IC) 2 is identified by the skew of the gradients as the inelastic component of subsidence.

![086A_05410_131313_output_demonstration_labels](https://github.com/user-attachments/assets/94d8b57f-e5f8-4a3e-b5c5-9bde9ddf06b5)
