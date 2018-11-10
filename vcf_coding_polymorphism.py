#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import argparse
import os


class Cds(object):
    def __init__(self, chromosome, strand, name):
        self.chromosome = chromosome
        self.strand = strand
        self.name = name
        self.exons = []

    def add_exon(self, start_exon, end_exon):
        if int(start_exon) <= int(end_exon):
            if self.strand == "+":
                self.exons.append((int(start_exon), int(end_exon)))
            else:
                self.exons.insert(0, (int(start_exon), int(end_exon)))
            for i in range(len(self.exons) - 1):
                if not self.exons[i][1] < self.exons[i + 1][0]:
                    print("At least one exon is overlapping with an other")

    def nt_position(self, position):
        nt = -1
        running_sum = 0
        if self.strand == "+":
            for start, end in self.exons:
                if start <= position <= end:
                    nt = position - start + running_sum
                    break
                running_sum += end - start + 1
        else:
            for start, end in reversed(self.exons):
                if start <= position <= end:
                    nt = end - position + running_sum
                    break
                running_sum += end - start + 1
        return nt

    
    def amino_acid(self, seq, position, nt_ref, nt_alt):
        """
        Returns the codon to which the SNP belongs and the corresponding amino 
        acid both for the reference and the individual of interest. 
        """
        assert (self.seq_length() % 3 == 0)
        assert (self.seq_length() == len(seq))
        # assert: stop if the conditions are not fulfilled
        
        position_relative = self.nt_position(position)
        
        if position_relative == -1:
            return("", "", "", "")
        
        # To read the sequence in the good direction
        if self.strand == "-":
            nt_ref = complement[nt_ref]
            nt_alt = complement[nt_alt]
            
        reste = position_relative % 3
        ref_codon = str(seq[position_relative - reste : position_relative - reste + 3])
        alt_codon = list(ref_codon)
        alt_codon[reste] = nt_alt
        
        if nt_ref != ref_codon[reste]: # si le VCF est en désaccord avec le FASTA
            return('Error FASTA VS VCF', "", "", "") # renvoie une sortie unique que l'on peut identifier comme une erreur
            
        ref_aa = codontable[ref_codon]
        alt_aa = codontable["".join(alt_codon)]
        
        return ref_aa, alt_aa, ref_codon, alt_codon


    def empty_exon(self):
        return sum([1 for exon_len in self.exons_length() if exon_len == 1]) > 0

    def exons_length(self):
        return [j - i + 1 for i, j in self.exons]

    def seq_length(self):
        return sum(self.exons_length())


def build_dict_cds(data_path, file_name):
    gtf_file = open("{0}/{1}".format(data_path, file_name), 'r')
    dico_cds = {}
    nf_tr_id = set()
    for line in gtf_file:
        line_split = line.replace('\n', '').split('\t')
        if len(line_split) > 7:
            info = line_split[8]
            if line_split[2] == 'CDS':
                transcript_find = info.find('transcript_id')
                if transcript_find != -1 and info.find('CCDS') != -1 and info.find('ccds_id') != -1:
                    tr_id = info[transcript_find + 15:].split("\"")[0]
                    if info.find('cds_start_NF') != -1 or info.find('cds_end_NF') != -1:
                        if tr_id not in nf_tr_id:
                            nf_tr_id.add(tr_id)
                    if tr_id not in dico_cds:
                        dico_cds[tr_id] = Cds(line_split[0], line_split[6], tr_id)
                    dico_cds[tr_id].add_exon(line_split[3], line_split[4])

    gtf_file.close()
    return dico_cds, nf_tr_id


def most_common(lst):
    """
    lst = liste indiquant la nature de la mutation pour chaque transcrit
    """
    return max(set(lst), key=lst.count)


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
codontable = defaultdict(lambda: "-")
codontable.update({
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'})

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta', required=True, type=str,
                        dest="f", metavar="<fasta>",
                        help="The relative name of the .fasta file")
    parser.add_argument('-g', '--gtf', required=True, type=str,
                        dest="g", metavar="<gtf>",
                        help="The relative name of the .gtf file")
    parser.add_argument('-v', '--vcf', required=True, type=str,
                        dest="v", metavar="<vcf>",
                        help="The relative name of the .vcf file")
    args = parser.parse_args()

    path = os.getcwd()
    dict_cds, not_confirmed_tr = build_dict_cds(path, args.g)

    dict_fasta = {}
    for fasta in SeqIO.parse(open("{0}/{1}".format(path, args.f), 'r'), 'fasta'):
        dict_fasta[fasta.id.split(".")[0]] = fasta.seq[:-3]

    stop_filename = '{0}/{1}.Stop.vcf'.format(path, args.v[:-4])
    syn_filename = '{0}/{1}.Syn.vcf'.format(path, args.v[:-4])
    nonsyn_filename = '{0}/{1}.NonSyn.vcf'.format(path, args.v[:-4])
    error_filename = '{0}/{1}.errors.tsv'.format(path, args.v[:-4])

    stop_file = open(stop_filename, 'w')
    syn_file = open(syn_filename, 'w')
    nonsyn_file = open(nonsyn_filename, 'w')
    error_file = open(error_filename, 'w')

    dict_cat_info = {"Syn": "{0} SNPs are synonymous variations",
                     "NonSyn": "{0} SNPs are non-synonymous variations",
                     "Stop": "{0} SNPs are stop variations",
                     "NotInFasta": "{0} Tr_id are not in fasta",
                     "RefStop": "{0} SNPs have stop codon as reference amino-acid",
                     "RefDiff": "{0} SNPs retrieved from the fasta are not equal to the reference",
                     "NotIdentified": "{0} SNPs have non-identified reference or alternate amino-acid",
                     "NotInCds": "{0} SNPs are not inside the CDS"}

    dict_cat_nbr = {}
    for snp_cat in dict_cat_info:
        dict_cat_nbr[snp_cat] = 0

    vcf_file = open("{0}/{1}".format(path, args.v), 'r')
    for line in vcf_file:
        if line[0] != '#':
            split_line = line.split("\t")
            if len(split_line[3]) == 1 and len(split_line[4]) == 1 and split_line[4] != ".":
                snp_id, chromosome, pos, ref_nt, alt_nt = split_line[2], split_line[0], int(split_line[1]), \
                                                          split_line[3], split_line[4]
                snp_types = []
                for tr_id in split_line[-1].replace('\n', '').replace('; ', ', ').split(","):
                    cds = dict_cds[tr_id]
                    if tr_id in dict_fasta:
                        ref_aa, alt_aa, ref_codon, alt_codon = cds.amino_acid(dict_fasta[tr_id], pos, ref_nt, alt_nt)
                        if ref_aa == '':
                            snp_types.append("NotInCds")
                        elif ref_aa == "Error FASTA VS VCF":
                            snp_types.append("RefDiff")
                        elif ref_aa == '-' or alt_aa == '-':
                            snp_types.append("NotIdentified")
                        elif ref_aa == 'X':
                            snp_types.append("RefStop")
                        elif alt_aa == 'X':
                            snp_types.append("Stop")
                        elif alt_aa == ref_aa:
                            snp_types.append("Syn")
                        else:
                            snp_types.append("NonSyn")
                    else:
                        print(tr_id + " not found in Fasta")
                        snp_types.append("NotInFasta")
                        continue

                assert len(snp_types) > 0
                if "NotInCds" in snp_types:
                    dict_cat_nbr["NotInCds"] += 1
                elif "NotInFasta" in snp_types:
                    dict_cat_nbr["NotInFasta"] += 1
                elif "RefDiff" in snp_types:
                    dict_cat_nbr["RefDiff"] += 1
                elif "NotIdentified" in snp_types:
                    dict_cat_nbr["NotIdentified"] += 1
                elif "RefStop" in snp_types:
                    dict_cat_nbr["RefStop"] += 1
                else:
                    max_type = most_common(snp_types)
                    if max_type == "Stop":
                        stop_file.write(line)
                        dict_cat_nbr["Stop"] += 1
                    elif max_type == "Syn":
                        syn_file.write(line)
                        dict_cat_nbr["Syn"] += 1
                    elif max_type == "NonSyn":
                        nonsyn_file.write(line)
                        dict_cat_nbr["NonSyn"] += 1
        else:
            stop_file.write(line)
            syn_file.write(line)
            nonsyn_file.write(line)
    vcf_file.close()

    nbr_snp_total = sum(dict_cat_nbr.values())
    error_file.write("{0} SNPs in total".format(nbr_snp_total))
    for cat, nbr in dict_cat_nbr.items():
        error_file.write("\n\n" + dict_cat_info[cat].format(nbr) + " ({0:.3f}%)".format(nbr * 100. / nbr_snp_total))

    stop_file.close()
    syn_file.close()
    nonsyn_file.close()
    error_file.close()

    print("{0} variants analyzed in total".format(nbr_snp_total))
    print("File containing {0} stop variants in {1}".format(dict_cat_nbr["Stop"], stop_filename))
    print("File containing {0} synonymous variants in {1}".format(dict_cat_nbr["Syn"], syn_filename))
    print("File containing {0} non-synonymous variants in {1}".format(dict_cat_nbr["NonSyn"], nonsyn_filename))
    nbr_errors = nbr_snp_total - dict_cat_nbr["Stop"] - dict_cat_nbr["Syn"] - dict_cat_nbr["NonSyn"]
    print("File containing {0} errors variants in {1}".format(nbr_errors, error_filename))
