### Instructions to run FastICA for subsidence on several LiCS frames ###

1. Run LiCSBAS for your frames of interest \
  a. Use at coherence mask threshold of at least 0.3

2. Organise your inputs for ICA based on outputs from LiCSBAS \
  a. Create a data_dir and specify path in ica_run.py. This is where you store your inputs for ICA. \
  b. Create an out_dir and specify path in ica_run.py. This is where your outputs from ICA will be saved. \
  c. Create several directories within the data_dir: \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;aa. \$data_dir/cumh5 \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bb. \$data_dir/mask \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cc. \$data_dir/EQA.dem\_par \
  d. Save your LiCSBAS outputs in these directories with the below structures \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;aa. cum.h5 --> ${data}\_dir/cumh5/${frame}\_GEOCml\*GACOS\*\_cum.h5 \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;n.b./ this is found in the TS_GEOCml* directory \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bb. mask.tif --> ${data}\_dir/mask/${frame}\_GEOCml\*GACOS\*mask\*\_coh\_03\_mask.geo.tif \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;n.b./ create this using LiCSBAS_flt2geotiff.py \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cc. EQA.dem_par --> ${data}\_dir/EQA.dem_par/${frame}GEOCml\*GACOS\*mask\*\_EQA.dem_par \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;n.b./ this is found in the relevant GEOCml* directory

4. Ensure you are running on slurm scheduler or similar that runs batch processing \
  a. Edit batch_execute_ica_run.sh to suit your batch processing system \
  b. Edit execute_ica_run.sh to choose the list of LiCS frames for which you want to perform ICA (FRAMES_LIST) \
  c. Edit VARIABLES TO CHANGE in ica_run.py to choose: \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;aa. Number of components to isolate (n_components) \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bb. Temporal or spatial ICA (refer to Ebmeier et al 2016) \
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cc. Again ensure directory paths are correct

5. Finally, run the batch script in command line: \
sbatch batch_execture_ica_run.sh
