#!/usr/bin/env python3

import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--panel', required=True, type=str,
                        dest="p", metavar="<panel>",
                        help="The relative name of the .panel file")
    parser.add_argument('-k', '--key', required=True, type=str,
                        dest="k", metavar="<keyword>",
                        help="The keyword that must match each line of the .panel file.")
    args = parser.parse_args()

    path = os.getcwd()
    call_samples = open("{0}/{1}".format(path, args.p), 'r')
    call_samples.readline()
    output_file = open("{0}/{1}.txt".format(path, args.k), 'w')
    nbr_indiv = 0
    for line in call_samples:
        if args.k in line.replace("\n", "").split("\t"):
            nbr_indiv += 1
            output_file.write(line.split("\t")[0]+"\n")
    output_file.close()
    call_samples.close()
    print("{0} individuals extracted from {1} population".format(nbr_indiv, args.k))
