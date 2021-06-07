#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np

def options():
    parser = argparse.ArgumentParser(description='Repalce geneID based on personalized file')
    parser.add_argument('--TE', type=str, help='TEs annotation gff file (three columns)')
    parser.add_argument('--MODs', type=str, help='MODs annotations')
    parser.add_argument('--Gene', type=str, help='Gene annotation gff file ')
    parser.add_argument('--window', type=int32, help='window size')
    parser.add_argument('--stepsize', type=int32, help='stepsize for each movement')

    args = parser.parse_args()
    return args

def main():
    # Get options
    args = options()

    #TE file
    TEs = pd.read_csv(args.TE,index_col=False, names=["Chr","str", "end", "Repeat_class"], sep = "\t")
    #MODs file
    MODs = pd.read_csv(args.MODs, index_col = False, sep = ",")
    #Gene file
    Genes = pd.read_table(args.Gene, index_col=False, names=["Chr","str","end"],sep = "\t")

    #re-sturcture MODs file 
    MODs_clean = MODs.loc[:,["Chr","bp","MODs_type","nonref","ref"]]
    MODs_clean['Total'] = MODs_clean['ref'] + MODs_clean['nonref']

    #Create group for each chromosome or genomic region of interests
    Group = pd.unique(MODs_clean['Chr'])

    #Create sliding window for each chromsomes
    dd = {}
    dd["info"] = pd.DataFrame({"Chr":Group, "TEs":[s + "_TEs"  for s in Group], 
                  "MODs":[s + "_MODs"  for s in Group], "window": [s + "_window" for s in Group]})

    for i in Group:
        dd[str(i) + "_window"] = pd.DataFrame({'Start':range(0,max(MODs_clean.loc[MODs_clean['Chr'] == i,'bp']),args.stepsize), 
                      "End": range(args.window,max(MODs_clean.loc[MODs_clean['Chr'] == i,'bp']) + args.window,args.stepsize)})

        dd[str(i) + "_MODs"] = MODs_clean.loc[MODs_clean['Chr'] == i,["Chr","bp","nonref","Total"]]
        dd[str(i) + "_TEs"] = TEs.loc[TEs['Chr'] == i,["str","end"]]
        dd[str(i) + "_Genes"] = Genes.loc[TEs['Chr'] == i,["str","end"]]

    #Calculate frequency for each type of data and bind to corresponded column 
    for a in Group:
    for b in ["window"]:
        for c in ["TEs"]:
            for d in ["MODs"]:
                for e in ["Genes"]:
                
                    window = a+"_"+b
                    TEs = a+"_"+c
                    MODs = a+"_"+d
                    Genes = a+"_"+e
            
                    for x in range(0,len(dd[window])):
                        dd[window].loc[x,'TEs_freq']= len(dd[TEs].loc[(dd[TEs]['str'] > dd[window].loc[x,'Start']) &  
                                                       (dd[TEs]['end'] <= dd[window].loc[x,'End']),:])
                
                        dd[window].loc[x,'MODs_freq']= len(dd[MODs].loc[(dd[MODs]['bp'] > dd[window].loc[x,'Start']) &  
                                                       (dd[MODs]['bp'] <= dd[window].loc[x,'End']),:])
                        
                        dd[window].loc[x,'Gene_density']= len(dd[Genes].loc[(dd[Genes]['str'] > dd[window].loc[x,'Start']) &  
                                                       (dd[Genes]['end'] <= dd[window].loc[x,'End']),:])

if __name__ == '__main__':
    main()
