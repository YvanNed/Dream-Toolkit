#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
Created on Thu Sep 29 14:42:30 2022
@author: arthurlecoz

event_xml_PD.py

Goal : gather the different events from the .xml file from profusion :
    -> Loop over the files of a directory
    -> Concatenate the information per subject
        -> info = 
            - name of event
            - start
            - duration
            - channel
    -> Put it in dataframe (.csv) : name_of_the_file_event_xml.csv
    -> save it in same dir

=============================================================================
"""
# %% Paths and Directory

import os, pathlib, xmltodict, pandas as pd

dir_to_xml = '/Users/thandrillon/Data/ICEBERG/PSG/EDF PD'
path_to_dir = pathlib.Path(dir_to_xml)

# %% Script

for filename in os.listdir(dir_to_xml):
    if filename.endswith(".XML"):
        print(filename)
        
        namelist = list(); startlist = list(); durationlist = list(); #channellist = list();
        
        path = path_to_dir/filename
        
        with open(
                path, 'r'
                ) as f:
            xmltext = f.read()
            
        xml = xmltodict.parse(
            xmltext, process_namespaces=True
            )['CMPStudyConfig']['ScoredEvents']
        
        for i, elem in enumerate(xml['ScoredEvent']) :
                namelist.append(elem['Name'])
                startlist.append(elem['Start'])
                durationlist.append(elem['Duration'])
               # channellist.append(elem['Input'])
        
        df = pd.DataFrame(
            {
                'Name' : namelist,
                'Start': startlist,
                'Duration': durationlist,
                #'Channel': channellist
                }
            )
        
        savename = filename[:-8] + "_event_xml.csv"
        df.to_csv(path_to_dir/savename)