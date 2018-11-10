#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# Input: 3 .vcf files (Stop, Syn, Non-Syn) + metadata file 
# Output: a .tsv file containing the number individuals, the number of of LoF SNPs taken into account in the analysis, Sigma^2/Va, and the result of the statistical test for each population 


RED = "#EB6231"
BLUE = "#5D80B4"
GREEN = "#8FB03E"


def build_array_from_vcf(path, vcf_name, _dict_pop, _dict_sample, _dict_sample_index, _cutoff):
    array_pop_dict = dict()
    for population, sample_list in _dict_pop.items():
        array_pop_dict[population] = []

    header_list = []
    vcf_file = open("{0}/{1}".format(path, vcf_name), 'r')
    for vcf_line in vcf_file:
        if vcf_line[0] == "#":
            if vcf_line[1] != "#":
                header_list = vcf_line.replace("\n", "").split("\t")
        else:
            line_list = vcf_line.replace("\n", "").split("\t")

            line_pop_dict = dict()
            for population, sample_list in _dict_pop.items():
                line_pop_dict[population] = np.zeros(len(sample_list), dtype=np.int8)

            for header_value, value_str in zip(header_list, line_list):
                if header_value in _dict_sample:
                    population = _dict_sample[header_value]

                    value = value_str.count("1")
                    if value > 0:
                        sample_index = _dict_sample_index[header_value]
                        line_pop_dict[population][sample_index] = value

            for population, line_pop in line_pop_dict.items():
                if _cutoff >= np.sum(line_pop) > 0:
                    array_pop_dict[population].append(line_pop)

    vcf_file.close()

    for population in array_pop_dict:
        array_pop_dict[population] = np.array(array_pop_dict[population], dtype=np.int8)

    return array_pop_dict


def epistasis(array):
    va = np.sum([np.var(col) for col in array])
    s2 = np.var([np.sum(row) for row in array.T])
    return s2 / va


def bootstrap(nbr_snps, array, bootstrap_resample=1000):
    """
    input:
        nbr_snps: size of the subsamples
        array: vcf
        bootstrap_resample: number of bootstrap samples
    ouput: list_ratio
        
    """
    row_nb, col_nb = np.shape(array)
    assert(nbr_snps<=row_nb) #checks if the instance is correct
    
    list_ratio = [] #initialization of the output
    for i in range(bootstrap_resample): 
        sample_idx = np.random.choice(np.array(range(row_nb)), nbr_snps, replace = False) # randomly select nbr_snps rows in array 
        
        subarray = np.zeros((nbr_snps, col_nb)) # initialization of subarray
        for m in range(col_nb):
            for n in range(nbr_snps):
                subarray[n][m]=array[sample_idx[n]][m]
        list_ratio.append(epistasis(subarray))
        
    return np.array(list_ratio)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--stop', required=True, type=str,
                        dest="stop", metavar="<vcf>",
                        help="The relative name of the stop .vcf file")
    parser.add_argument('-s', '--syn', required=True, type=str,
                        dest="syn", metavar="<vcf>",
                        help="The relative name of the synonymous .vcf file")
    parser.add_argument('-n', '--nonsyn', required=True, type=str,
                        dest="non_syn", metavar="<vcf>",
                        help="The relative name of the non-synonymous .vcf file")
    parser.add_argument('-p', '--panel', required=True, type=str,
                        dest="p", metavar="<panel>",
                        help="The relative name of the .panel file")
    parser.add_argument('-c', '--cutoff', required=False, default=1, type=int,
                        dest="c", metavar="<panel>",
                        help="The cut-off for minor allele count")

    args = parser.parse_args()

    dict_pop = {}
    dict_sample_index = {}
    dict_sample = {}

    call_samples = open("{0}/{1}".format(os.getcwd(), args.p), 'r')
    call_samples.readline()
    for line in call_samples:
        split_line = line.split("\t")
        pop = split_line[1]
        sample = split_line[0]
        if pop not in dict_pop:
            dict_pop[pop] = list()
        dict_sample_index[sample] = len(dict_pop[pop])
        dict_pop[pop].append(sample)
        dict_sample[sample] = pop
    call_samples.close()

    stop_array_dict = build_array_from_vcf(os.getcwd(), args.stop, dict_pop, dict_sample, dict_sample_index, args.c)
    print("Computed stop variants array for the whole population.")
    syn_array_dict = build_array_from_vcf(os.getcwd(), args.syn, dict_pop, dict_sample, dict_sample_index, args.c)
    print("Computed synonymous variants array for the whole population.")
    non_syn_array_dict = build_array_from_vcf(os.getcwd(), args.non_syn, dict_pop, dict_sample, dict_sample_index,
                                              args.c)
    print("Computed non-synonymous variants array for the whole population.")

    tsv_file = open("{0}/meta_analysis_{1}.tsv".format(os.getcwd(), args.c), 'w')
    tsv_file.write('\t'.join(['Population', 'NbrIndividuals', 'NbrStops', 's^2/Va', 'Pvalue', 'Significant']) + '\n')

    for pop, stop_array in stop_array_dict.items():
        print("{0} population with {1} individuals".format(pop, len(dict_pop[pop])))

        stop_epistasis = epistasis(stop_array)

        nbr_stops = stop_array.shape[0]
        if nbr_stops > 0:
            syn_bootstrap = bootstrap(nbr_stops, syn_array_dict[pop])
            non_syn_bootstrap = bootstrap(nbr_stops, non_syn_array_dict[pop])

            my_dpi = 128
            fig = plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            p_val = len([1 for epi in syn_bootstrap if epi < stop_epistasis]) / syn_bootstrap.shape[0]

            result_line = [pop, str(len(dict_pop[pop])), str(nbr_stops), str(stop_epistasis), str(p_val)]
            if p_val < 0.05:
                result_line += ["True"]
                print("\t p-value={0:3g}, the statistical test is significant ".format(p_val))
            else:
                result_line += ["False"]
                print("\t p-value={0:3g}, the statistical test is not significant ".format(p_val))
            tsv_file.write('\t'.join(result_line) + '\n')

            plt.title('{0} ($p={1:3g}$)'.format(pop, p_val))
            bins = 25
            syn_hist, _, _ = plt.hist(syn_bootstrap, bins, density=1, facecolor=BLUE, alpha=0.4, label='Syn')
            non_syn_hist, _, _ = plt.hist(non_syn_bootstrap, bins, density=1, facecolor=GREEN, alpha=0.4,
                                          label='NonSyn')
            y_max = 1.2 * max((max(syn_hist), max(non_syn_hist)))
            plt.ylim((0, y_max))
            plt.plot((stop_epistasis, stop_epistasis), (0, y_max), label="LoF", linewidth=3, color=RED)
            plt.xlabel('$\sigma^2/V_{A}$')
            plt.ylabel('density')
            plt.legend()
            plt.tight_layout()
            plt.savefig("{0}/analysis_{1}_{2}.svg".format(os.getcwd(), pop, args.c), format="svg")
            plt.savefig("{0}/analysis_{1}_{2}.png".format(os.getcwd(), pop, args.c), format="png")
            plt.close()
        else:
            print('\tNo stop variants for this population')

    tsv_file.close()
    print("Analysis completed")
