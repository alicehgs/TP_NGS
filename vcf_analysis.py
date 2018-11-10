#!/usr/bin/env python3

import numpy as np
import argparse
import os

# Input: a .vcf file  
# Output: a .tsv file containing the number of SNPs taken into account in the analysis, the number of individuals, Va, Sigma^2 and the ratio Sigma^2/Va


def var(array_1d, n):
    """
    Returns variance of an array of length n
    """
    x_2_bar = sum([x * x for x in array_1d]) / n
    x_bar = sum(array_1d) / n
    return x_2_bar - x_bar * x_bar


def variance_additive(array):
    """
    Returns the variance of the sum
    """
    # Calcul de Va qui est la variance de la somme
    ligne, colonne = array.shape
    liste_var = []
    for i in range(ligne):
        liste_var.append(var(array[i,:], colonne))
    return sum(liste_var)


def sigma_squared(array):
    """
    Returns the sum of the variances
    """
    ligne, colonne = array.shape
    liste_score = []
    for j in range(colonne):
        liste_score.append(sum(array[:,j]))
    return var(liste_score, colonne)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--vcf', required=True, type=str,
                        dest="v", metavar="<vcf>",
                        help="The relative name of the .vcf file")
    parser.add_argument('-c', '--cutoff', required=False, default=1, type=int,
                        dest="c", metavar="<cutoff>",
                        help="The cut-off for rare alleles")
    # Utilise le cut-off par défaut: nombre de 1 voulus
    parser.add_argument('-f', '--freq', required=False, default=0, type=float,
                        dest="f", metavar="<freq>",
                        help="The frequency cut-off for rare alleles")
    # Utilise la fréquence si on le demande: on veut tant de % de variants dans la population
    args = parser.parse_args()

    variant_list = []

    path = os.getcwd()

    vcf_file = open("{0}/{1}".format(path, args.v), 'r')
    nbr_indiv = 0
    for vcf_line in vcf_file:
        # Pour chaque SNP, faire la boucle: regarde SNP par SNP
        if vcf_line[0] == "#":
            if vcf_line[1] != "#":
                nbr_indiv = len(vcf_line.split("\t")) - 9
                print("{0} individuals found in the vcf file.".format(nbr_indiv))
        else:
            line_list = vcf_line.replace('\n', '').split("\t")[9:9 + nbr_indiv]

            sparse_line = np.zeros(nbr_indiv, dtype=np.int8)
            for i, genotype in enumerate(line_list):
                # Pour chaque individu, regarder son génotype au SNP étudié
                value = genotype.count("1")
                # Compte le nb de 1: 0 (0/0 homozygote référence),1 (0/1 ou 1/0 hétérozygote), 2 (1/1 homozygote varaint)
                if value > 0:
                    # Ajoute le génotype à value et si value sup à 0, conserver la valeur de l'individu et l'ajouter à la somme de tous les individus
                    sparse_line[i] = value
            if args.f > 0:
                if args.f * nbr_indiv > np.sum(sparse_line) > 0:
                    variant_list.append(sparse_line)
                    # Si la somme pour le SNP est comprise entre notre fréquence et 0, garder le SNP
            elif args.c >= np.sum(sparse_line) > 0:
                # Si la somme pour le SNP est comprise entre le cut-off et 0, garder le SNP
                variant_list.append(sparse_line)
    vcf_file.close()

    sparse_array = np.array(variant_list, dtype=np.int8)
    nbr_snps = sparse_array.shape[0]
    print("{0} SNPs found in the .vcf file".format(nbr_snps))

    tsv_file = open("{0}/{1}_{2}.analysis.tsv".format(path, args.v[:-4], args.c), 'w')
    tsv_file.write('\t'.join(['SNP', 'sample_size', 's^2', 'Va', 's^2/Va']) + '\n')

    va = variance_additive(sparse_array)

    print("Va={0}".format(va))

    s2 = sigma_squared(sparse_array)

    print("\sigma^2={0}".format(s2))
    print("\sigma^2/Va={0}".format(s2 / va))
    tsv_file.write('\t'.join([str(nbr_snps), str(nbr_indiv), str(s2), str(va), str(s2 / va)]) + '\n')
    tsv_file.close()
    print("Analysis performed")


