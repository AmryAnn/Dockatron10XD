# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:27:22 2020

@author: aarid

Analysis blind docking results
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import argparse
import logging


#""" Load separated results into individual dataframes and then concatenate
#those dataframes into one dataframe. Added a column 'Index' that assigns a 
#unique number to each record in the dataframe."""
#
#ls1 = pd.read_csv("ligand_binding_sites5_5000.txt", sep='\t')
#ls2 = pd.read_csv("ligand_binding_sites5005_10000.txt", sep='\t')
#ls3 = pd.read_csv("ligand_binding_sites10005_15000.txt", sep='\t')
#ls4 = pd.read_csv("ligand_binding_sites15005_20000.txt", sep='\t')
#ls5 = pd.read_csv("ligand_binding_sites20005_25000.txt", sep='\t')
#ls6 = pd.read_csv("ligand_binding_sites25005_30000.txt", sep='\t')
#ls7 = pd.read_csv("ligand_binding_sites30005_35000.txt", sep='\t')
#ls8 = pd.read_csv("ligand_binding_sites35005_40200.txt", sep='\t')
#
#LS_all = pd.concat([ls1, ls2, ls3, ls4, ls5, ls6, ls7, ls8])
#
#LS_all['Index'] = range(len(LS_all))


""" Various analysis functions: basic affinity score stats, records with the 
top ten affinty scores over all and for each ligand, frequency of sites bound 
over all and for each ligand. """ 

def LigandSites_DataFrame(ligand_sites):
    df = pd.read_csv(ligand_sites, sep='\t')
    df['Index'] = range(len(df))
    print(df)
    return df

def ligand_Pose_energy_stats(ligand_sites):
    print("\n Basic affinity score statistics:")
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
    print("Number of unique affinity scores: ", nunique_affinity, "\n")
    logging.info("\n Basic affinity score statistics:\n"
                 f"Number of records: {len(ligand_sites)}\n"
                 f"Average affinity score: {mean_affinity}\n"
                 f"Median affinity score: {median_affinity}\n"
                 f"Standard deviation: {affinity_std}\n"
                 f"Best affinity score: {min_affinity}\n"
                 f"Worst affinity score: {max_affinity}\n"
                 f"Number of unique affinity scores: {nunique_affinity}\n")
    return mean_affinity, median_affinity, affinity_std, max_affinity, \
min_affinity, nunique_affinity

#ligand_Pose_energy_stats(LS_all)

def top_ten_over_all(ligand_sites):
    sortE = ligand_sites.sort_values('Energy')
    worst_value = sortE.head(10).Energy.max()
    top_ten_scores = ligand_sites[ligand_sites.Energy <= worst_value]
    top_ten_sorted = top_ten_scores.sort_values('Energy')
    top_ten_sorted['Rank'] = range(1, len(top_ten_sorted)+1)
    print("Top ten affinity scores over all: ")
    print(top_ten_sorted)
    logging.info(f"\n Top ten affinity scores over all: \n{top_ten_sorted} \n")
    
    return top_ten_sorted

#top_ten_over_all(LS_all)

def top_ten_for_ligand(ligand, ligand_sites):
    ligand_data = ligand_sites[ligand_sites.Lig == ligand]
    ligand_data_sorted = ligand_data.sort_values('Energy')
    worst_value = ligand_data_sorted.head(10).Energy.max()
    top_ten = ligand_data[ligand_data.Energy <= worst_value]
    top_ten_sorted = top_ten.sort_values('Energy')
    top_ten_sorted['Rank'] = range(1, len(top_ten_sorted)+1)
    print("\n Top ten affinity scores for ligand: ", ligand)
    print(top_ten_sorted)
    logging.info(f"\n Top ten affinity scores for ligand {ligand}: \n{top_ten_sorted} \n")
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
    print("\n Number of Poses at Each Site:")
    print(freq)
    logging.info(f"\n Number of Poses at Each Site: \n{freq} \n")
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
    print("\n Number of Poses at Each Site for Ligand: "+str(ligand))
    print(freq)
    logging.info(f"\n Number of poses at each site for ligand {ligand}: \n{freq} \n")
    return sites, freq

#freq_sites_bound_for_ligand('L23', LS_all)
    
def best_affinity_at_each_site(ligand_sites):
    sites_best = ligand_sites.groupby('Site').Energy.min()
    print("\n Best affinity score at each site: ")
    print(sites_best)
    logging.info(f"\n Best affinity score at each site: \n{sites_best} \n")
    return sites_best

#best_affinity_at_each_site(LS_all)
    
def best_affinity_at_each_site_for_ligand(ligand, ligand_sites):
    ligand_data = ligand_sites[ligand_sites.Lig == ligand]
    sites_best = ligand_data.groupby('Site').Energy.min()
    print("\n Best affinity score at each site for ligand", ligand+":")
    print(sites_best)
    logging.info(f"\n Best affinity score at each site for ligand {ligand}: \n{sites_best} \n")
    return sites_best

#best_affinity_at_each_site_for_ligand('L23', LS_all)


""" Command line arguments and logfile """

if __name__ == '__main__':
    # Add command line argument parser
    parser = argparse.ArgumentParser(description="Further analysis of ensemble docking results.", \
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ligand_sites', type=str, required=True, \
                        help='The path to ligand_binding_sites.txt')
    parser.add_argument('--ligand_ID', type=str, required=True, \
                        help='The identifier for the ligand you are analyzing')
    
    # Add logger to create a log file
    logging.basicConfig(filename="DockatronAnalysis2.log", format="%(asctime)s - %(levelname)s - %(message)s",
                        level=logging.NOTSET)
    logging.basicConfig(filename="Step2.log", format="%(asctime)s - %(levelname)s - %(message)s",
                        level=logging.NOTSET)
    logging.getLogger('matplotlib.font_manager').disabled = True
    logger = logging.getLogger(__name__)
    logging.info("*****Begin ensemble docking analysis part 2*****")
    
    # Assign command line arguments
    args = parser.parse_args()
    SITES_PATH = Path(args.ligand_sites)
    if not SITES_PATH.exists() or not SITES_PATH.is_file():
        logging.error(f"ligand_binding_sites.txt path {SITES_PATH} is not a valid file or directory")
        exit()
    LIG = args.ligand_ID
        
    # Set logging level
    
    # Function calls
    DataFrame = LigandSites_DataFrame(SITES_PATH)
    ligand_Pose_energy_stats(DataFrame)
    top_ten_over_all(DataFrame)
    freq_sites_bound(DataFrame)
    best_affinity_at_each_site(DataFrame)
    top_ten_for_ligand(LIG, DataFrame)
    freq_sites_bound_for_ligand(LIG, DataFrame)
    best_affinity_at_each_site_for_ligand(LIG, DataFrame)
    logging.info("\n Finished Analysis Part 2")
    print(" \n Finished Analysis Part 2 \n")
    
    