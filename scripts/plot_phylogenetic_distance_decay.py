import  matplotlib.pyplot as plt

import parse_midas_data
#import pylab
#from pylab import *
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
import os
#import pandas
import parse_patric
import pickle
import sample_utils
import parse_HMP_data

import itertools

import scipy.stats as stats

import ete3
from itertools import combinations
#kegg_pi_Eubacterium_rectale_56927.dat

tree_path = "%stree/RAxML_bestTree.HMP" % parse_midas_data.data_directory
tree = ete3.Tree(tree_path, quoted_node_names=True, format=1)

leaf_name_dict = {}

for leaf in tree.iter_leaf_names():
    leaf_name_dict[leaf.split('|')[0]] = leaf


good_species_list = parse_midas_data.parse_good_species_list()
#good_species_list = ["Eubacterium_rectale_56927"]

kegg_dict_all_species = {}

for species_name in good_species_list:

    kegg_path = '%skegg_pi/kegg_pi_%s.dat' % (parse_midas_data.data_directory, species_name)

    if os.path.exists(kegg_path) == False:
        continue

    kegg_dict = pickle.load( open(kegg_path, "rb" ) )

    pathways = kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys()

    kegg_dict_all_species[species_name] = {}

    for pathway in kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys():

        pathway_array = kegg_dict['fraction_nonsynonymous_per_pathway_core'][pathway]

        n_hosts = pathway_array.shape[0]

        values = pathway_array[numpy.triu_indices(n_hosts, 1)]

        #kegg_dict_all_species[species_name][pathway] = numpy.mean(values[values>0])

        #if pathway not in kegg_dict_all_species:
        #    kegg_dict_all_species[pathway] = []

        if len(values[(values>0) & (values<1)]) <= 10:
            continue

        values = values[(values>0) & (values<1)]

        values = values / (1-values)

        kegg_dict_all_species[species_name][pathway] = numpy.mean(values)

        #kegg_dict_all_species[pathway].append(numpy.mean(values[(values>0) & (values<1)]))



phylogenetic_distances = []
correlation_coefficients = []

species_pairs = list(combinations(kegg_dict_all_species.keys(), 2))

for species_pair in species_pairs:

    if (species_pair[0] not in leaf_name_dict.keys()) or (species_pair[1] not in leaf_name_dict.keys()):
        continue

    dict_1 = kegg_dict_all_species[species_pair[0]]
    dict_2 = kegg_dict_all_species[species_pair[1]]

    intersecting_genes = list(set(dict_1.keys()) & set(dict_2.keys()))

    if len(intersecting_genes) < 5:
        continue

    genes_1 = [ dict_1[g] for g in intersecting_genes]
    genes_2 = [ dict_2[g] for g in intersecting_genes]

    genes_1 = numpy.asarray(genes_1)
    genes_2 = numpy.asarray(genes_2)

    correlation_coefficients.append(numpy.corrcoef(genes_1, genes_2)[1,0])

    phylogenetic_distances.append(tree.get_distance(leaf_name_dict[species_pair[0]] ,leaf_name_dict[species_pair[1]]) )


correlation_coefficients = numpy.asarray(correlation_coefficients)
phylogenetic_distances = numpy.asarray(phylogenetic_distances)


slope, intercept, r_value, p_value, std_err = stats.linregress(phylogenetic_distances, correlation_coefficients)
slope_permuted_list = []
regression_permutations=10000
for permute in range(regression_permutations):
    distances_permuted = numpy.random.permutation(phylogenetic_distances)
    correlations_permuted = numpy.random.permutation(correlation_coefficients)
    slope_permuted, intercept_permuted, r_value_permuted, p_valu_permutede, std_err_permuted = stats.linregress(distances_permuted, correlations_permuted)
    slope_permuted_list.append(slope_permuted)

slope_permuted_list = numpy.asarray(slope_permuted_list)

p_value = len(slope_permuted_list[slope_permuted_list<slope]) / regression_permutations

x_range =  numpy.linspace(0.1, 2, 10000)
y_fit_range = (slope*x_range + intercept)


fig, ax = plt.subplots(figsize = (4, 4))


ax.scatter(phylogenetic_distances, correlation_coefficients, alpha = 0.3 )
ax.plot(x_range, y_fit_range, c='k', lw=2.5, linestyle='--', zorder=2, label='OLS fit')

#ax.title("Distance-de", fontsize=14)
ax.axhline(y=0, color='k', linestyle=':', lw=2, alpha = 0.8, zorder=1, label='Corr.' + r'$=0$')

ax.set_xlabel('Phylogenetic distance bw. species pair' , fontsize = 12)
ax.set_ylabel('Correlation of nonsynonymous ratios\nacross pathways bw. species pair' , fontsize = 12)

#ax.text(0.75, 0.9, r'$y \sim x^{{{}}} \; P < 10^{{-6}}$'.format(str( round(slope, 3) ) ) , fontsize=11, color='k', ha='center', va='center', transform=ax.transAxes  )
ax.text(0.71, 0.9, r'$\beta_{{1}} ={{{}}} \; P < 10^{{-6}}$'.format(str( round(slope, 3) ) ) , fontsize=11, color='k', ha='center', va='center', transform=ax.transAxes  )

ax.legend(loc='lower left',frameon=True)

fig.subplots_adjust(hspace=0.4, wspace=0.35) #hspace=0.3, wspace=0.5
fig_name =  "%scorrelation_distance_decay.png" % parse_midas_data.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
