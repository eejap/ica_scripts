#!/bin/bash

# python interpreter
PYTHON_EXECUTABLE=/gws/smf/j04/nceo_geohazards/software/mambalics/bin/python

# Path to your Python script
SCRIPT_PATH=ica_run_region.py

# Path to the region list
REGIONS_LIST="/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/ne_iran_data/region_names.txt"

# Loop over each region in the regions list
while IFS= read -r region; do
    echo "Submitting job for region: $region"

    # Run your Python script with the current region
    $PYTHON_EXECUTABLE $SCRIPT_PATH --region "$region"
done < "$REGIONS_LIST"
