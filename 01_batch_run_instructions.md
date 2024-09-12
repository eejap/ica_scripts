### Instructions to run FastICA for subsidence on several LiCS frames ###

1. Run LiCSBAS for your frames of interest \
  a. Use at coherence mask threshold of at least 0.3

2. Organise your inputs for ICA based on outputs from LiCSBAS \
  a. Create a data_dir and specify path in ica_run.py. This is where you store your inputs for ICA. \
  b. Create an out_dir and specify path in ica_run.py. This is where your outputs from ICA will be saved. \
  c. Create several directories within the data_dir: \
    aa. \$data_dir/cumh5 \
      bb. \$data_dir/mask \
      cc. \$data_dir/EQA.dem\_par \
  d. Save your LiCSBAS outputs in these directories with the below structures \
      aa. cum.h5 --> ${data}\_dir/cumh5/${frame}\_GEOCml\*GACOS\*\_cum.h5                       n.b./ this is found in the TS_GEOCml* directory \
      bb. mask.tif --> ${data}\_dir/mask/${frame}\_GEOCml\*GACOS\*mask\*\_coh\_03\_mask.geo.tif    n.b./ create this using LiCSBAS_flt2geotiff.py \
      cc. EQA.dem_par --> ${data}\_dir/EQA.dem_par/${frame}GEOCml\*GACOS\*mask\*\_EQA.dem_par   n.b./ this is found in the relevant GEOCml* directory \

3. Ensure you are running on slurm scheduler or similar that runs batch processing
  a. Edit batch_execute_ica_run.sh to suit your batch processing system
  b. Edit execute_ica_run.sh to choose the list of LiCS frames for which you want to perform ICA (FRAMES_LIST)
  c. Edit VARIABLES TO CHANGE in ica_run.py to choose:
      aa. Number of components to isolate (n_components)
      bb. Temporal or spatial ICA (refer to Ebmeier et al 2016)
      cc. Again ensure directory paths are correct

4. Finally, run the batch script in command line:
sbatch batch_execture_ica_run.sh
