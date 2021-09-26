#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh

pyscript='/home/jovyan/repo1-jtstanley/counter.py'
ufofile='/home/jovyan/repo1-jtstanley/UFO-pop-data/ufo_sightings.csv'
colnum=3

outfile='/home/jovyan/repo1-jtstanley/state_counts.txt'

conda activate swefs
python $pyscript $ufofile $colnum > $outfile
