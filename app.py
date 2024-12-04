from pywebio.input import *
from pywebio.output import *
import os
import pandas as pd
import csv

def FileInput():
    method = radio("Choose one", options=['FASTQ', 'SRA'])
    
    if method == 'FASTQ':
        
        csv_file = file_upload("Please a list of FASTQ file in csv format", accept=".csv",  multiple = True,  placeholder='Choose a file')
        exp_design = file_upload("Please upload the experimental design in tab delimited format", accept=".txt",  multiple = True,  placeholder='Choose a file')
        
        put_text("Your files have been successfully inputed to the pipeline")
        
        list_control = []  # a list to hold all the individual pandas DataFrames
        for csvfile in csv_file:
            open(csvfile["filename"], "wb").write(csvfile["content"])
            df = pd.read_csv(csvfile["filename"])
            list_control.append(df)
            
        list_contrast = []  # a list to hold all the individual pandas DataFrames
        for csvfile in exp_design:
            open(csvfile["filename"], "wb").write(csvfile["content"])
            df = pd.read_csv(csvfile["filename"], sep = "\t")
            list_contrast.append(df)
            
        list_control[0].to_csv('./csv_file.csv', index = False)
        list_contrast[0].to_csv('./exp_design.txt',  sep='\t', index = False) 
        
    if method == 'SRA':
        img = file_upload("Select an image:", accept="image/*")
        
if __name__ == '__main__':
    FileInput()
