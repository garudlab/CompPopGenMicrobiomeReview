from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
#from mpl_toolkits.axes_grid.inset_locator import inset_axes
import pylab

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from numpy.random import randint, binomial, choice, poisson
from scipy.stats import poisson as poisson_distribution
import scipy.stats as stats

from collections import defaultdict
import sys
import os
import config
import numpy
from random import sample, shuffle
#import numpy as np

import pickle

import plot_utils
import snps_utils
import sample_utils
import parse_midas_data

from matplotlib.patches import Rectangle


num_bootstraps = 100
#num_bootstraps = 2
n_percentile_permutations = 10000
#n_percentile_permutations = 10

sub_plot_labels = ['a','b','c']


data_dir = config.data_directory
picles_dir = "%spickles/" % data_dir
picles_dir = "%spickles/" % data_dir

#syn_differences.pkl

syn_differences = pickle.load(open("%sbetween_syn_differences.pkl" % (picles_dir), 'rb'))
syn_opportunities = pickle.load(open("%sbetween_syn_opportunities.pkl" % (picles_dir), 'rb'))
syn_pseudocounts = pickle.load(open("%sbetween_syn_pseudocounts.pkl" % (picles_dir), 'rb'))
non_differences = pickle.load(open("%sbetween_non_differences.pkl" % (picles_dir), 'rb'))
non_pseudocounts = pickle.load(open("%sbetween_non_pseudocounts.pkl" % (picles_dir), 'rb'))
non_opportunities = pickle.load(open("%sbetween_non_opportunities.pkl" % (picles_dir), 'rb'))


species_color_map, ordered_species_list = plot_utils.get_species_color_map()


def tp_to_category(tp_pair):
		tpa, tpb = tp_pair
		string = 'MI' if (tpa[0], tpb[0]) == ('I', 'M') else tpa[0]+tpb[0]
		return string

# Simplified version

pSs_by_tp_cat = defaultdict(list)
pNpSs_by_tp_cat = defaultdict(list)
pSs_by_species_tp_cat = defaultdict(dict)
pNpSs_by_species_tp_cat = defaultdict(dict)

all_syn_differences_by_tp_cat = defaultdict(list)
all_syn_opportunities_by_tp_cat = defaultdict(list)
all_non_differences_by_tp_cat = defaultdict(list)
all_non_opportunities_by_tp_cat = defaultdict(list)
#all_core_differences_by_tp_cat = defaultdict(list)
all_core_opportunities_by_tp_cat = defaultdict(list)

for tp_pair in syn_differences:
	try:
		tp_cat = tp_to_category(tp_pair)
	except:
		continue

	for species_name in syn_differences[tp_pair]:

		syn_opps = syn_opportunities[tp_pair][species_name]
		syn_diffs = syn_differences[tp_pair][species_name]
		non_opps = non_opportunities[tp_pair][species_name]
		non_diffs = non_differences[tp_pair][species_name]

		pSs = syn_diffs*1.0/syn_opps
		pNs = non_diffs*1.0/non_opps
		pseudo_pSs = 1.0/(syn_opps/2.0+non_opps)
		pseudo_pNs = 1.0/(syn_opps/2.0+non_opps)

		pNpSs = ((pseudo_pNs+pNs)/(pseudo_pSs+pSs))

		good_idxs = ((syn_diffs+non_diffs)>=10)
		pSs_by_tp_cat[tp_cat] += list(pSs[good_idxs])
		pNpSs_by_tp_cat[tp_cat] += list(pNpSs[good_idxs])

		all_syn_differences_by_tp_cat[tp_cat].extend(syn_diffs[good_idxs])
		all_syn_opportunities_by_tp_cat[tp_cat].extend(syn_opps[good_idxs])
		all_non_differences_by_tp_cat[tp_cat].extend(non_diffs[good_idxs])
		all_non_opportunities_by_tp_cat[tp_cat].extend(non_opps[good_idxs])

		try:
			pSs_by_species_tp_cat[species_name][tp_cat] += list(pSs[good_idxs])
			pNpSs_by_species_tp_cat[species_name][tp_cat] += list(pNpSs[good_idxs])
		except:
			pSs_by_species_tp_cat[species_name][tp_cat] = list(pSs[good_idxs])
			pNpSs_by_species_tp_cat[species_name][tp_cat] = list(pNpSs[good_idxs])

median_pSs_by_tp_cat = defaultdict(list)
median_pNpSs_by_tp_cat = defaultdict(list)

for species_name in pNpSs_by_species_tp_cat:
	for tp_cat in pNpSs_by_species_tp_cat[species_name]:

		median_pS = numpy.median(pSs_by_species_tp_cat[species_name][tp_cat])
		median_pSs_by_tp_cat[tp_cat].append(median_pS)

		median_pNpS = numpy.median(pNpSs_by_species_tp_cat[species_name][tp_cat])
		if median_pNpS > 0.2:
			print(tp_cat + ": " + species_name + " has median pNpS " + str(median_pNpS))
		median_pNpSs_by_tp_cat[tp_cat].append(median_pNpS)




# Bootstrapping dN/dS

avg_cf_ratios = defaultdict(list)
std_cf_ratios = defaultdict(list)
median_cf_ratios = defaultdict(list)
lower_cf_ratios = defaultdict(list)
upper_cf_ratios = defaultdict(list)
avg_sf_ratios = defaultdict(list)
std_sf_ratios = defaultdict(list)

ds = numpy.logspace(-5,-2,50)

#for tp_cat in ['II', 'AA', 'MI', 'MM']:
for tp_cat in ['AA']:

	cf_ratios = [] # cumulative estimates <= total d
	sf_ratios = [] # cumulative estimates >= total d

	sys.stderr.write("Bootstrapping dN/dS...\n")

	for bootstrap_idx in range(num_bootstraps):

		lower_pNpSs, upper_pNpSs = [], []

		all_non_diffs = numpy.array(all_non_differences_by_tp_cat[tp_cat])
		all_non_opps = numpy.array(all_non_opportunities_by_tp_cat[tp_cat])
		all_syn_diffs = numpy.array(all_syn_differences_by_tp_cat[tp_cat])
		all_syn_opps = numpy.array(all_syn_opportunities_by_tp_cat[tp_cat])
		all_NS_differences = all_syn_diffs + all_non_diffs

		# Bootstrap dataset using poisson resampling
		# Each observed difference is considered lambda for Poisson
		# distribution. Resample according to pmf in which output N
		# is weighted by probability that that number of differences
		# occurs within interval lambda, the average number of diffs.

		# When lambda is small, Poisson distribution is highly
		# right skewed. As lambda approaches infinity, become
		# more and more like the binomial distribution

		# Pseudocounts so things w/ 0 counts are not "stuck" in resampling
		# Pseudocounts are chosen w/ dN/dS=1, so should be conservative?
		# (alternatively, we could choose dN/dS=0.1 -- unfair?)

		pseudocount = 0 # 1.0
		bs_non_differences = poisson(all_non_diffs + pseudocount)
		bs_syn_differences = poisson(all_syn_diffs + (all_syn_opps*pseudocount/all_non_opps))

		bs_NS_differences = bs_non_differences + bs_syn_differences

		# Cut down numbers by half on average
		bs_thinned_syn_differences_1 = binomial(bs_syn_differences, 0.5)
		bs_thinned_syn_differences_2 = bs_syn_differences - bs_thinned_syn_differences_1

		# Bootstrapped dS
		bs_divergence = bs_thinned_syn_differences_1 / (all_syn_opps/2.0)

		for d in ds:

			lower_idxs = (bs_divergence <= d)*(all_NS_differences>0.5)*(bs_NS_differences>0.5)
			upper_idxs = (bs_divergence > d)*(all_NS_differences>0.5)*(bs_NS_differences>0.5)

			if lower_idxs.sum()<1.5:
					lower_pNpSs.append(-1)
			else:
					lower_cumulative_non_differences = (bs_non_differences)[lower_idxs].sum()
					lower_cumulative_expected_non_differences = (bs_thinned_syn_differences_2[lower_idxs]*2.0/all_syn_opps[lower_idxs]*all_non_opps[lower_idxs]).sum()
					lower_pNpSs.append( (lower_cumulative_non_differences)/(lower_cumulative_expected_non_differences) )

			if upper_idxs.sum()<1.5:
					upper_pNpSs.append(-1)
			else:
					upper_cumulative_non_differences = (bs_non_differences[upper_idxs]).sum()
					upper_cumulative_expected_non_differences = (bs_thinned_syn_differences_2[upper_idxs]*2.0/all_syn_opps[upper_idxs]*all_non_opps[upper_idxs]).sum()
					upper_pNpSs.append( (upper_cumulative_non_differences)/(upper_cumulative_expected_non_differences) )

		cf_ratios.append(lower_pNpSs)
		sf_ratios.append(upper_pNpSs)

		if bootstrap_idx % 10 == 0:
			print("On bootstrap %i..." % bootstrap_idx)

	cf_ratios = numpy.array(cf_ratios)
	sf_ratios = numpy.array(sf_ratios)

	for i in range(len(ds)):

		ratios = numpy.sort(cf_ratios[:,i])
		good_idxs = (ratios>-0.5)
		if good_idxs.sum()<1.5:
				avg_cf_ratios[tp_cat].append(-1)
				std_cf_ratios[tp_cat].append(0)

		else:
				median_cf_ratios[tp_cat].append(numpy.median(ratios[good_idxs]))
				idx = long(0.025*good_idxs.sum())
				lower_cf_ratios[tp_cat].append( ratios[good_idxs][idx] )
				upper_cf_ratios[tp_cat].append(ratios[good_idxs][-idx-1])

				avg_cf_ratios[tp_cat].append( ratios[good_idxs].mean() )
				std_cf_ratios[tp_cat].append( ratios[good_idxs].std() )

		ratios = sf_ratios[:,i]
		good_idxs = (ratios>-0.5)
		if good_idxs.sum()<1.5:
				avg_sf_ratios[tp_cat].append(-1)
				std_sf_ratios[tp_cat].append(0)
		else:
				avg_sf_ratios[tp_cat].append( ratios[good_idxs].mean() )
				std_sf_ratios[tp_cat].append( ratios[good_idxs].std() )

	avg_cf_ratios[tp_cat] = numpy.array(avg_cf_ratios[tp_cat])
	std_cf_ratios[tp_cat] = numpy.array(std_cf_ratios[tp_cat])
	median_cf_ratios[tp_cat] = numpy.array(median_cf_ratios[tp_cat])
	upper_cf_ratios[tp_cat] = numpy.array(upper_cf_ratios[tp_cat])
	lower_cf_ratios[tp_cat] = numpy.array(lower_cf_ratios[tp_cat])
	avg_sf_ratios[tp_cat] = numpy.array(avg_sf_ratios[tp_cat])
	std_sf_ratios[tp_cat] = numpy.array(std_sf_ratios[tp_cat])

#tp_cat_descrip = {'II': 'infant-infant', 'AA': 'adult-adult', 'MI': 'mother-infant', 'MM': 'mother-mother'}
tp_cat_descrip = {'AA': 'adult-adult'}


# Create a species color legend
all_species = set()
for species in pSs_by_species_tp_cat:
	all_species.add(species)

all_species_ordered = []
colors_ordered = []
for species in ordered_species_list:
	if species in all_species:
		all_species_ordered.append(species)
		colors_ordered.append(species_color_map[species])



gs = gridspec.GridSpec(nrows=2, ncols=2)

fig = plt.figure(figsize = (20, 20))

divergence_axis = fig.add_subplot(gs[0, :])
kegg_axis = fig.add_subplot(gs[1, 0])
percentile_axis = fig.add_subplot(gs[1, 1])


divergence_axis.text(-0.1, 1.07, sub_plot_labels[0], fontsize=16, fontweight='bold', ha='center', va='center', transform=divergence_axis.transAxes)
kegg_axis.text(-0.1, 1.07, sub_plot_labels[1], fontsize=16, fontweight='bold', ha='center', va='center', transform=kegg_axis.transAxes)
percentile_axis.text(-0.1, 1.07, sub_plot_labels[2], fontsize=16, fontweight='bold', ha='center', va='center', transform=percentile_axis.transAxes)


#ax.set_position([box.x0, box.y0, box.width * 0.2, box.height])
#ax.legend(loc='center left', handles=plot_utils.colors_to_legend_elements(colors_ordered, all_species_ordered), ncol=3, fontsize='medium', bbox_to_anchor=(1, 0.5))
#fig.savefig('%s/species_color_legend.pdf' % (config.analysis_directory), bbox_inches='tight')

# Color species
#for tp_cat in ['II', 'AA', 'MI', 'MM']:
#for tp_cat in ['AA']:

tp_cat = 'AA'

#fig, divergence_axis = plt.subplots(figsize=(14, 8))

divergence_axis.set_ylabel('Nonsynonymous ratio, $d_N/d_S$', fontsize=28)
divergence_axis.set_xlabel('Synonymous divergence, $d_S$', fontsize=28)
#divergence_axis.spines['top'].set_visible(False)
#divergence_axis.spines['right'].set_visible(False)
#divergence_axis.get_xaxis().tick_bottom()
#divergence_axis.get_yaxis().tick_left()


# Inset for cumulative dN/dS
#cumulative_axis = inset_axes(divergence_axis, width="25%", height="25%", borderpad=0, bbox_to_anchor=(-0.01,0,1, 1), bbox_transform=divergence_axis.transAxes)

#cumulative_axis.spines['top'].set_visible(False)
#cumulative_axis.spines['right'].set_visible(False)
#cumulative_axis.get_xaxis().tick_bottom()
#cumulative_axis.get_yaxis().tick_left()

#cumulative_axis.set_ylabel('Cumulative $d_N/d_S$')
#cumulative_axis.set_xlabel('Synonymous divergence, $d_S$')

#line, = cumulative_axis.loglog([1e-05,1e-02],[1,1],'k:',linewidth=0.25,zorder=1)
#line.set_dashes((1,1))

#good_idxs = (avg_cf_ratios[tp_cat]>-0.5)
#cumulative_axis.fill_between(ds[good_idxs], lower_cf_ratios[tp_cat][good_idxs], upper_cf_ratios[tp_cat][good_idxs],color='0.7',linewidth=0,zorder=0)
#cumulative_axis.loglog(ds[good_idxs], avg_cf_ratios[tp_cat][good_idxs],'k-',zorder=2)

#cumulative_axis.set_xlim([1e-05,1e-02])
#cumulative_axis.set_ylim([5e-02,2])

# Moving on...
all_pSs = pSs_by_tp_cat[tp_cat]
all_pNpSs = pNpSs_by_tp_cat[tp_cat]

for species in pSs_by_species_tp_cat:
	if tp_cat in pSs_by_species_tp_cat[species]:
		pSs = pSs_by_species_tp_cat[species][tp_cat]
		pNpSs = pNpSs_by_species_tp_cat[species][tp_cat]
		color = species_color_map[species]
		divergence_axis.loglog(pSs, pNpSs, '.', color=color, markersize=6,alpha=0.6,markeredgewidth=0,zorder=0,rasterized=True)

divergence_axis.plot([1e-09],[100], 'o', color='darkgrey', markersize=6,markeredgewidth=0,zorder=0,label='(Species x host x host)')

median_pSs = median_pSs_by_tp_cat[tp_cat]
median_pNpSs = median_pNpSs_by_tp_cat[tp_cat]

divergence_axis.loglog(median_pSs, median_pNpSs, 'kx',markersize=12, markeredgewidth=4, label='Median of each species',alpha=0.9)

#divergence_axis.set_ylim([1e-02,10])
#divergence_axis.set_xlim([1e-06,1e-01])

divergence_axis.set_ylim([0.5e-01,1.7])
divergence_axis.set_xlim([1e-05,1e-01])


# Purifying selection model
# Fitted manually
asymptotic_dNdS = 0.12
dStar = 3e-04
sbymu = 1/dStar/asymptotic_dNdS
print "s/u =", sbymu
print "s =", sbymu*1e-09

def theory_dN(dS):
	return (asymptotic_dNdS+(1-asymptotic_dNdS)*(1-numpy.exp(-sbymu*dS))/(theory_ds*sbymu))*dS

theory_ds = numpy.logspace(-6,-1,100)
theory_dNdSs = theory_dN(theory_ds)/theory_ds

line, = divergence_axis.loglog([1e-06,1e-01],[1,1],'k:',linewidth=3,label='Neutral model')
line.set_dashes((1,1))
divergence_axis.loglog(theory_ds, theory_dNdSs,'k--',linewidth=3,label='Purifying selection model')

divergence_axis.legend(loc='lower left',frameon=True,numpoints=1)

#divergence_axis.set_title("dN/dS vs. dS, all %s sample pairs (n=%i) (%i species)" % (tp_cat_descrip[tp_cat], len(all_pSs), len(median_pSs)))





##################################
# plot dn/ds of functional groups#
##################################

good_species_list = parse_midas_data.parse_good_species_list()
#good_species_list = ["Eubacterium_rectale_56927"]

kegg_dict_all_species = {}

for species in good_species_list:


	kegg_path = '%skegg_pi/kegg_pi_%s.dat' % (parse_midas_data.data_directory, species)

	if os.path.exists(kegg_path) == False:
		continue

	kegg_dict = pickle.load( open(kegg_path, "rb" ) )

	#pathways = kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys()

	for pathway in kegg_dict['fraction_nonsynonymous_per_pathway_core'].keys():

		if (pathway == '') or (pathway == 'Annotated pathways'):
			continue

		pathway_array = kegg_dict['fraction_nonsynonymous_per_pathway_core'][pathway]

		n_hosts = pathway_array.shape[0]

		values = pathway_array[numpy.triu_indices(n_hosts, 1)]

		values = values[(values>0) & (values<1)]

		values = values / (1-values)

		#if pathway not in kegg_dict_all_species:
		#	kegg_dict_all_species[pathway] = []

		if len(values[(values>0) & (values<1)]) <= 20:
			continue
		#kegg_dict_all_species[pathway].append(numpy.mean(values[(values>0) & (values<1)]))

		if pathway not in kegg_dict_all_species:
			kegg_dict_all_species[pathway] = {}
			kegg_dict_all_species[pathway]['mean_dnds_list'] = []
			kegg_dict_all_species[pathway]['species_list'] = []
			kegg_dict_all_species[pathway]['colors_list'] = []

		kegg_dict_all_species[pathway]['mean_dnds_list'].append(numpy.mean(values))
		kegg_dict_all_species[pathway]['species_list'].append(species)
		kegg_dict_all_species[pathway]['colors_list'].append(species_color_map[species])

		#kegg_dict_all_species[pathway].append()


fraction_nonsynonymous_per_pathway_core_means = []

pathway_names = kegg_dict_all_species.keys()

paired_means = [ [pathway_name, numpy.mean(kegg_dict_all_species[pathway_name]['mean_dnds_list']), numpy.var(kegg_dict_all_species[pathway_name]['mean_dnds_list']) ] for pathway_name in kegg_dict_all_species.keys() if (len(kegg_dict_all_species[pathway_name]['mean_dnds_list']) > 10) ]

paired_means_sorted = sorted(paired_means, key = lambda x: float(x[1]))


# axis
species_dict = {}
# make dict for each species
for pathway, pathway_dict in kegg_dict_all_species.items():

	mean_dnds_list_pathway = kegg_dict_all_species[pathway]['mean_dnds_list']

	for species_j, mean_dnds_j in zip(kegg_dict_all_species[pathway]['species_list'], mean_dnds_list_pathway):

		if species_j not in species_dict:
			species_dict[species_j] = {}
			species_dict[species_j]['pathways_list'] = []
			species_dict[species_j]['mean_dnds_list'] = []

		species_dict[species_j]['mean_dnds_list'].append(mean_dnds_j)
		species_dict[species_j]['pathways_list'].append(pathway)



mean_dnds_across_species_iter_dict = {}
for i in range(n_percentile_permutations): # n_percentile_permutations
	mean_dnds_across_species_dict = {}
	for species_i, dnds_dict_i in species_dict.items():
		dnds_list_permuted = dnds_dict_i['mean_dnds_list']

		shuffle(dnds_list_permuted)

		for dnds_i_permuted, pathway_i_permuted in zip(dnds_list_permuted, dnds_dict_i['pathways_list']):

			if pathway_i_permuted not in mean_dnds_across_species_dict:
				mean_dnds_across_species_dict[pathway_i_permuted] = []

			mean_dnds_across_species_dict[pathway_i_permuted].append(dnds_i_permuted)

	for pathway_i, dnds_i_list in mean_dnds_across_species_dict.items():

		mean_dnds_i = numpy.mean(dnds_i_list)

		if pathway_i not in mean_dnds_across_species_iter_dict:
			mean_dnds_across_species_iter_dict[pathway_i] = []

		mean_dnds_across_species_iter_dict[pathway_i].append(mean_dnds_i)



dnds_null_025_list = []
dnds_null_975_list = []
for pathway_i, dnds_i in mean_dnds_across_species_iter_dict.items():

	dnds_i.sort()

	dnds_null_025_list.append(dnds_i[int(n_percentile_permutations*0.025)])
	dnds_null_975_list.append(dnds_i[int(n_percentile_permutations*0.975)])

dnds_null_025 = numpy.mean(dnds_null_025_list)
dnds_null_975 = numpy.mean(dnds_null_975_list)


#all_dnds_null_025 = [ mean_dnds_across_species_dict[k][int(n_percentile_permutations*0.025)] for k in  mean_dnds_across_species_dict.keys()]
#all_dnds_null_975 = [ mean_dnds_across_species_dict[k][int(n_percentile_permutations*0.975)] for k in  mean_dnds_across_species_dict.keys()]


#all_dnds_null = [item for sublist in all_dnds_null for item in sublist]
#all_dnds_null.sort()
#dnds_null_025 = all_dnds_null[int(len(all_dnds_null) * 0.025)]
#dnds_null_975 = all_dnds_null[int(len(all_dnds_null) * 0.975)]

mean_dnds_all_species = []
var_dnds_all_species = []

for paired_means_i_idx, paired_means_i in enumerate(paired_means_sorted):

	mean_dnds_list_kegg = kegg_dict_all_species[paired_means_i[0]]['mean_dnds_list']
	colors_list_kegg = kegg_dict_all_species[paired_means_i[0]]['colors_list']

	mean_kegg_pathway = numpy.mean(mean_dnds_list_kegg)
	var_kegg_pathway = numpy.var(mean_dnds_list_kegg)

	mean_dnds_all_species.append(mean_kegg_pathway)
	var_dnds_all_species.append(var_kegg_pathway)

	#for zip_i in zip(mean_dnds_list_kegg, colors_list_kegg):

	location = numpy.asarray([paired_means_i_idx] * len(mean_dnds_list_kegg))


	kegg_axis.scatter(mean_dnds_list_kegg, location,
	        s=18, \
			linewidths=0, \
			rasterized=True, \
	        linewidth=1, facecolors=colors_list_kegg, \
	        marker='o', \
	        alpha=0.7, zorder=3)

	#kegg_axis.scatter(zip_i[0], paired_means_idx,
	#        s=6, \
	#		linewidths=0, \
	#		rasterized=True, \
	#        linewidth=1, facecolors=zip_i[1], \
	#        marker='.', \
	#        alpha=0.7, zorder=1)

	kegg_axis.scatter(mean_kegg_pathway, paired_means_i_idx,
	        s = 20, \
	        linewidth=1, facecolors='k', \
	        edgecolors='k', marker='o', \
	        alpha=1, zorder=4)


kegg_axis.set_xlim(0.02, 1.7)
kegg_axis.set_ylim(-1, len(paired_means_sorted)+1)

kegg_axis.axvline(x=1, color='k', linestyle=':', alpha = 0.8, linewidth=3, label='Neutral model', zorder=1)
kegg_axis.axvline(x=dnds_null_025, color='k', linestyle='--', alpha = 0.8, linewidth=3, zorder=1)
kegg_axis.axvline(x=dnds_null_975, color='k', linestyle='--', alpha = 0.8, linewidth=3, label='95% CI', zorder=1)

kegg_axis.set_xscale('log', basex=10)
#kegg_axis.set_xlabel('Mean $d_N/d_S$ across hosts', fontsize=26)
kegg_axis.set_xlabel('Nonsynonymous ratio, $d_N/d_S$', fontsize=26)


labels = [x[0] for x in paired_means_sorted]

kegg_axis.tick_params(axis='y', which='y',length=0)


kegg_axis.set_yticks(numpy.arange(len(labels)))
kegg_axis.set_yticklabels(labels, rotation=30)
kegg_axis.yaxis.set_tick_params(labelsize=7)



kegg_axis.plot([1e-09],[100], 'o', color='darkgrey', markersize=6,markeredgewidth=0,zorder=0,label='Mean of all hosts for a species')
kegg_axis.plot([1e-09],[100], 'o', color='k', markersize=6,markeredgewidth=0,zorder=0,label='Mean of all species')

kegg_axis.legend(loc='upper left',frameon=True, numpoints=1)






means = [x[1] for x in paired_means_sorted]
variances = [x[2] for x in paired_means_sorted]


ax_taylors = inset_axes(kegg_axis, width="100%", height="100%", loc='lower right', bbox_to_anchor=(0.68,0.07,0.3,0.3), bbox_transform=kegg_axis.transAxes)

ax_taylors.scatter(means, variances, s = 50, \
        linewidth=1.5, facecolors='k', \
        edgecolors='k', marker='o', \
        alpha=0.9)


slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(means), numpy.log10(variances))


x_log10_fit_range =  numpy.linspace(numpy.log10(min(means) * 0.5), numpy.log10(max(means) * 1.5), 10000)

y_fit_range = 10 ** (slope*x_log10_fit_range + intercept)
ax_taylors.plot(10**x_log10_fit_range, y_fit_range, c='k', lw=3, linestyle='--', zorder=2)

ax_taylors.text(0.33,0.9, r'$\sigma^{{2}}_{{ d_{{N}}/d_{{S}} }} \propto \left \langle d_{{N}}/d_{{S}} \right \rangle^{{{}}}$'.format(str(round(slope, 2)) ), fontsize=11, color='k', ha='center', va='center', transform=ax_taylors.transAxes  )

ax_taylors.set_xscale('log', basex=10)
ax_taylors.set_yscale('log', basey=10)

#ax_taylors.set_xlim(0.02, 1.7)


ax_taylors.set_xlabel('Mean $d_N/d_S$, ' +  r'$\left \langle d_{N}/d_{S} \right \rangle$', fontsize=10)
ax_taylors.set_ylabel('Variance of $d_N/d_S$, ' +  r'$\sigma^{2}_{d_{N}/d_{S}}$', fontsize=10)


# patch for y axis label
kegg_axis.add_patch(Rectangle((0.57 , 0.1), 0.072, 0.25, alpha=1, facecolor='white', transform=kegg_axis.transAxes, zorder=5))


# patch for x axis label
kegg_axis.add_patch(Rectangle((0.73 , 0.018), 0.18, 0.032, alpha=1, facecolor='white', transform=kegg_axis.transAxes, zorder=5))


##################################
# plot percentile of each species
##################################


#stats.percentileofscore([1, 2, 3, 4], 3)
percentile_dict = {}

for pathway, pathway_dict in kegg_dict_all_species.items():

	mean_dnds_list_pathway = kegg_dict_all_species[pathway]['mean_dnds_list']

	if len(mean_dnds_list_pathway) <10:
		continue

	for species_j, mean_dnds_j in zip(kegg_dict_all_species[pathway]['species_list'], mean_dnds_list_pathway):

		if species_j not in percentile_dict:
			percentile_dict[species_j] = {}
			percentile_dict[species_j]['percentiles_list'] = []
			percentile_dict[species_j]['pathways'] = []

		# we want fractile
		#percentile = stats.percentileofscore(mean_dnds_list_pathway, mean_dnds_j) / 100

		#percentile_dict[species_j]['percentiles_list'].append(percentile)
		percentile_dict[species_j]['percentiles_list'].append(mean_dnds_j)

		percentile_dict[species_j]['pathways'].append(pathway)



# generate null distribution
mean_percentile_null_dict = {}
for i in range(n_percentile_permutations):

	percentile_null_dict = {}

	for pathway, pathway_dict in kegg_dict_all_species.items():

		mean_dnds_list_pathway = kegg_dict_all_species[pathway]['mean_dnds_list']
		if len(mean_dnds_list_pathway) <10:
			continue

		shuffle(mean_dnds_list_pathway)

		for species_j, mean_dnds_j in zip(kegg_dict_all_species[pathway]['species_list'], mean_dnds_list_pathway):

			if species_j not in percentile_null_dict:
				percentile_null_dict[species_j] = {}
				percentile_null_dict[species_j]['percentiles_list'] = []
				percentile_null_dict[species_j]['pathways'] = []

			#percentile = stats.percentileofscore(mean_dnds_list_pathway, mean_dnds_j) / 100

			#percentile_null_dict[species_j]['percentiles_list'].append(percentile)
			percentile_null_dict[species_j]['percentiles_list'].append(mean_dnds_j)

			percentile_null_dict[species_j]['pathways'].append(pathway)

	percentile_means_null = [ [species_j, numpy.mean(percentile_null_dict[species_j]['percentiles_list']) ] for species_j in percentile_null_dict.keys() ]

	for species_mean_percentile_list in percentile_means_null:

		if species_mean_percentile_list[0] not in mean_percentile_null_dict:
			mean_percentile_null_dict[species_mean_percentile_list[0]] = []

		mean_percentile_null_dict[species_mean_percentile_list[0]].append(species_mean_percentile_list[1])


# generate confidence interval dict
conf_interval_percentile_null_dict = {}
CI_025 = []
CI_975 = []
for species_i, species_i_null_list in mean_percentile_null_dict.iteritems():

	species_i_null_list.sort()

	conf_interval_percentile_null_dict[species_i] = {}
	CI_025_i = species_i_null_list[int(n_percentile_permutations*0.025)]
	CI_025.append(CI_025_i)
	CI_975_i = species_i_null_list[int(n_percentile_permutations*0.975)]
	CI_975.append(CI_975_i)

	conf_interval_percentile_null_dict[species_i]['CI_025'] = CI_025_i
	conf_interval_percentile_null_dict[species_i]['CI_975'] = CI_975_i

CI_025 = numpy.asarray(CI_025)
CI_975 = numpy.asarray(CI_975)
CI_025_mean = numpy.mean(CI_025)
CI_975_mean = numpy.mean(CI_975)




for species_i, percentile_dict_i in percentile_dict.items():

	if len(percentile_dict_i['percentiles_list']) < 10:
		del percentile_dict[species_i]



percentile_means = [ [species_j, numpy.mean(percentile_dict[species_j]['percentiles_list']) ] for species_j in percentile_dict.keys() ]


percentile_means_sorted = sorted(percentile_means, key = lambda x: float(x[1]))

for percentile_means_i_idx, percentile_means_i in enumerate(percentile_means_sorted):


	percentile_means_i_list = percentile_dict[percentile_means_i[0]]['percentiles_list']

	location = numpy.asarray([percentile_means_i_idx] * len(percentile_means_i_list))

	#percentile_axis.scatter(percentile_means_i_list, location,
	#        s=23, \
	#		linewidths=0.5, \
	#		rasterized=True, \
	#        linewidth=1, facecolors=species_color_map[percentile_means_i[0]], \
	#        marker='o', \
	#        alpha=0.5, zorder=2)


	#percentile_axis.scatter(percentile_means_i[1], percentile_means_i_idx,
	#        s = 25, \
	#        linewidth=1, facecolors='k', \
	#        edgecolors='k', marker='o', \
	#        alpha=1, zorder=3

	percentile_axis.scatter(percentile_means_i[1], percentile_means_i_idx,
	        s = 52, \
	        linewidth=0, facecolors=species_color_map[percentile_means_i[0]], \
	        edgecolors='k', marker='o', \
	        alpha=1, zorder=3)









percentile_axis.axvline(x=CI_025_mean, color='k', linestyle='--', alpha = 0.8, linewidth=3, zorder=1)
percentile_axis.axvline(x=CI_975_mean, color='k', linestyle='--', alpha = 0.8, linewidth=3, label='95% CI', zorder=1)


#percentile_axis.axvline(x=0.5, color='k', linestyle=':', alpha = 0.8, linewidth=3, label='Null expectation', zorder=1)

#percentile_axis.set_xlabel('Mean $d_N/d_S$ fractile across pathways', fontsize=26)
#percentile_axis.set_xlabel('$d_N/d_S$ fractile', fontsize=26)


#ax_taylors.set_ylabel('Variance of $d_N/d_S$ across hostss, ' +  r'$\sigma^{2}_{d_{N}/d_{S}}$', fontsize=10)


#percentile_axis.plot([1e-09],[100], 'o', color='darkgrey', markersize=6,markeredgewidth=0,zorder=0,label='Mean of all hosts for a pathway')
#percentile_axis.plot([1e-09],[100], 'o', color='k', markersize=6,markeredgewidth=0,zorder=0,label='Mean of all pathways')
percentile_axis.plot([1e-09],[100], 'o', color='darkgrey', markersize=6,markeredgewidth=0,zorder=0,label='Mean of all pathways')


percentile_axis.axvline(x=1, color='k', linestyle=':', alpha = 0.8, linewidth=3, label='Neutral model', zorder=1)



labels = [plot_utils.latex_species_name_dict[x[0]] for x in percentile_means]

percentile_axis.tick_params(axis='y', which='y',length=0)


#percentile_axis.set_xlim(-0.02, 1.02)
percentile_axis.set_xlim(0.02, 1.7)
percentile_axis.set_ylim(-1, len(percentile_means_sorted))

percentile_axis.set_xscale('log', basex=10)

percentile_axis.set_xlabel('Nonsynonymous ratio, $d_N/d_S$', fontsize=26)


percentile_axis.set_yticks(numpy.arange(len(labels)))
percentile_axis.set_yticklabels(labels, rotation=30)
percentile_axis.yaxis.set_tick_params(labelsize=7.5)


percentile_axis.legend(loc='lower right',frameon=True,numpoints=1)





fig.subplots_adjust(hspace=0.18, wspace=0.24) #hspace=0.3, wspace=0.5
# bbox_inches = "tight",, pad_inches = 0.5, dpi = 600
# , dpi = 600
#fig.savefig("%s/dnds_by_species_%s_high_res.pdf" % (config.analysis_directory, tp_cat), bbox_inches = "tight", pad_inches = 0.5, dpi = 600)

fig.savefig("%s/dnds_by_species_%s.pdf" % (config.analysis_directory, tp_cat), bbox_inches = "tight", pad_inches = 0.5)
