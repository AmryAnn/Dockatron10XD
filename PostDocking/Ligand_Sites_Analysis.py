# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:27:22 2020

@author: aarid

Analysis blind docking results
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


""" Load separated results into individual dataframes and then concatenate
those dataframes into one dataframe. Added a column 'Index' that assigns a 
unique number to each record in the dataframe."""

ls1 = pd.read_csv("ligand_binding_sites5_5000.txt", sep='\t')
ls2 = pd.read_csv("ligand_binding_sites5005_10000.txt", sep='\t')
ls3 = pd.read_csv("ligand_binding_sites10005_15000.txt", sep='\t')
ls4 = pd.read_csv("ligand_binding_sites15005_20000.txt", sep='\t')
ls5 = pd.read_csv("ligand_binding_sites20005_25000.txt", sep='\t')
ls6 = pd.read_csv("ligand_binding_sites25005_30000.txt", sep='\t')
ls7 = pd.read_csv("ligand_binding_sites30005_35000.txt", sep='\t')
ls8 = pd.read_csv("ligand_binding_sites35005_40200.txt", sep='\t')

LS_all = pd.concat([ls1, ls2, ls3, ls4, ls5, ls6, ls7, ls8])

LS_all['Index'] = range(len(LS_all))


""" Various analysis functions: basic affinity score stats, records with the 
top ten affinty scores over all and for each ligand, frequency of sites bound 
over all and for each ligand. """ 

def ligand_Pose_energy_stats(ligand_sites):
    print("Basic affinity score statistics:")
    print("Number of records: ", len(ligand_sites))
    mean_affinity = ligand_sites.Energy.mean()
    print("Average affinity score: ", mean_affinity)
    median_affinity = ligand_sites.Energy.median()
    print("Median affinity score: ", median_affinity)
    affinity_std = ligand_sites.Energy.std()
    print("Standard deviation: ", affinity_std)
    max_affinity = ligand_sites.Energy.max()
    min_affinity = ligand_sites.Energy.min()
    print("Best affinity score: ", min_affinity)
    print("Worst affinity score: ", max_affinity)
    nunique_affinity = ligand_sites.Energy.nunique()
    print("Number of unique affinity scores: ", nunique_affinity)
    return mean_affinity, median_affinity, affinity_std, max_affinity, \
min_affinity, nunique_affinity

#ligand_Pose_energy_stats(LS_all)

def top_ten_over_all(ligand_sites):
    sortE = ligand_sites.sort_values('Energy')
    worst_value = sortE.head(10).Energy.max()
    top_ten_scores = ligand_sites[ligand_sites.Energy <= worst_value]
    top_ten_sorted = top_ten_scores.sort_values('Energy')
    print("Top ten affinity scores over all: ")
    print(top_ten_sorted)
    return top_ten_sorted

#top_ten_over_all(LS_all)

def top_ten_for_ligand(ligand, ligand_sites):
    ligand_data = ligand_sites[ligand_sites.Lig == ligand]
    ligand_data_sorted = ligand_data.sort_values('Energy')
    worst_value = ligand_data_sorted.head(10).Energy.max()
    top_ten = ligand_data[ligand_data.Energy <= worst_value]
    top_ten_sorted = top_ten.sort_values('Energy')
    print("Top ten affinity scores for ligand: ", ligand)
    print(top_ten_sorted)
    return top_ten_sorted

#top_ten_for_ligand('L23', LS_all)

def freq_sites_bound(ligand_sites):    
    sites = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    freq = ligand_sites.groupby('Site').Index.count()
    plt.bar(sites, freq)
    plt.title("Number of Poses at Each Site")
    plt.xlabel("Site")
    plt.ylabel("Number of Poses")
    plt.show()
    print(freq)
    return sites, freq

#freq_sites_bound(LS_all)

def freq_sites_bound_for_ligand(ligand, ligand_sites):
    ligand_data = ligand_sites[ligand_sites.Lig == ligand]
    sites = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    freq = ligand_data.groupby('Site').Index.count()
    plt.bar(sites, freq)
    plt.title("Number of Poses at Each Site for Ligand: "+str(ligand))
    plt.xlabel("Site")
    plt.ylabel("Number of Poses")
    plt.show()
    print(freq)
    return sites, freq

#freq_sites_bound_for_ligand('L23', LS_all)
    
def best_affinity_at_each_site(ligand_sites):
    sites_best = ligand_sites.groupby('Site').Energy.min()
    print("Best affinity score at each site: ")
    print(sites_best)
    return sites_best

#best_affinity_at_each_site(LS_all)
    
def best_affinity_at_each_site_for_ligand(ligand, ligand_sites):
    ligand_data = ligand_sites[ligand_sites.Lig == ligand]
    sites_best = ligand_data.groupby('Site').Energy.min()
    print("Best affinity score at each site for ligand -", ligand+":")
    print(sites_best)
    return sites_best

#best_affinity_at_each_site_for_ligand('L23', LS_all)

