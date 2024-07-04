#!/bin/bash

# python interpreter
PYTHON_EXECUTABLE=/gws/smf/j04/nceo_geohazards/software/mambalics/bin/python

# Path to your Python script
SCRIPT_PATH=ica_run.py

# Path to the frames list file
FRAMES_LIST="/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran_data/ica_frames_run.txt"

# Loop over each frame in the frames list
while IFS= read -r frame; do
    echo $frame
    # Run the Python script with the current frame
    $PYTHON_EXECUTABLE $SCRIPT_PATH --frame "$frame"
done < "$FRAMES_LIST"

