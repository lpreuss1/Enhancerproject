#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for turning values from enhancer program into
mean for every 50 bp window.
This script needs the csv file, as well as the bed file for the region 
as an input
"""

import pandas as pd
import argparse

def get_slice_value(ws, sts):
    '''generates slice value from windowsize(ws) and stepsize(sts)'''
    slicevalue = int(ws / sts) +1 
    return slicevalue

def score_for_point(df, a, sv):
    """Computes score at a certain point"""
    mean = df['score'].iloc[a-sv:a].mean()
    return mean

def all_scores(df,sv):
    '''returns the mean value for each step between all overlapping windows'''
    values = []
    for i,v in enumerate(df.score):
        if i == 0:
            values.append(df['score'].iloc[0])
        
        elif i < sv and i > 0:
            values.append(df['score'].iloc[0:i].mean())
            
        else:
            values.append(score_for_point(df, i, sv))
    return values

if __name__ == "__main__":
    
    #adding commandline Arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-ws', metavar=('-ws N'),
                        action= 'store',
                        type=int,
                        nargs='?',
                        default=500,
                        help='value for window-size\n default value = 500')
    parser.add_argument('-sts', metavar=('-sts N'),
                        action= 'store',
                        type=int,
                        nargs='?',
                        default=500,
                        help='value for step-size\n default value = 50')
    parser.add_argument('-s', '--scores', metavar=('\b'),
                        type=argparse.FileType('r'),
                        required = True,
                        help='Output csv-file from enhancer-program')
    parser.add_argument('-b', '--bed', metavar='\b',
                        type=argparse.FileType('r'),
                        required = True, 
                        help= 'bedfile of the scanned region')
    args = parser.parse_args()
    
    #create a data frame
    df = pd.read_csv(args.scores, sep= ' ')
    
    #Change Column names
    df = df.rename(columns = {'The':'steps', 'max': 'score'})
    
    #removing null values
    df = df.dropna(axis=1)
    
    # turning step-values into genome coordinates
    bed = pd.read_csv(args.bed, sep = '\t', 
                      names = ['chr', 'start', 'end'])
    
    df.steps = df.steps +   bed['start'].loc[0]
    
    
    #for slicing
    slicevalue = get_slice_value(args.ws, args.sts)
    mean_values = all_scores(df, slicevalue) 
    
    #creates filename
    filename = f"{str(args.scores.name)}.with_means"
    
    #returns a csv with the new values as a column
    df['median_scores'] = mean_values
    df['chr'] = bed.chr.loc[0]
    df[['steps', 'median_scores', 'chr']].to_csv(filename, index=False,
                                                 sep='\t')
