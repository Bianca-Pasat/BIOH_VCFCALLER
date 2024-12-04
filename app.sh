#!/bin/bash
python app.py;

wdir=$(pwd);
filename=${wdir}/csv_file.csv;
exp_design=${wdir}/exp_design.csv;

nextflow run '/home/bianca/Desktop/nf-core-qualitycontrol_align_longhack' -profile docker --input $filename --samplesInfo $exp_design --sra false -resume
