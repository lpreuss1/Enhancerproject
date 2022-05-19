#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 14:45:30 2022

@author: leonardopreuss
"""
"""
Script for turning values from enhancer program into
mean for every 50 bp window.
This script needs the csv file, as well as the bed file for the region 
as an input
"""

import pandas as pd
import argparse


#adding commandline Arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-s', '-scores', type=argparse.FileType('r'))
parser.add_argument('-b', '--bed', type=argparse.FileType('r'))
args = parser.parse_args()

#create a data frame
df = pd.read_csv(args.s, sep= ' ')

#Change Column names
df = df.rename(columns = {'The':'steps', 'max': 'score'})

#removing null values()
df = df.dropna(axis=1)

# turning tep-values into genome coordinates
bed = pd.read_csv(args.bed, sep = '\t', 
                  names = ['chr', 'start', 'end'])

df.steps = df.steps +   bed['start'].loc[0]
#defining step and window size as well as the value 
#fopr slicing
ws = 300 #windowsize
sts = 50 #stepsize

sv = int(ws / sts) +1 #slice_value

def score_for_point(df, a, sv):
    """Computes score at a certain point"""
    mean = df['score'].iloc[a-sv:a].mean()
    return mean

def all_scores(df,sv):
    values = []
    for i,v in enumerate(df.score):
        if i == 0:
            values.append(df['score'].iloc[0])
        
        elif i < sv and i > 0:
            values.append(df['score'].iloc[0:i].mean())
            
        else:
            values.append(score_for_point(df, i, sv))
    return values

values = all_scores(df, sv)

#creates filename
filename = f"{args.s}.with_means"

#returns a csv with the new values as a column
df['median_scores'] = values
df['chr'] = bed.chr.loc[0]
df[['steps', 'median_scores', 'chr']].to_csv(filename, index=False, sep='\t')
