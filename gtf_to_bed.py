#!/usr/bin/env python3

import argparse
import os


class Cds(object):
    # Définit de l'objet CDS
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

    def empty_exon(self):
        return sum([1 for exon_len in self.exons_length() if exon_len == 1]) > 0

    def exons_length(self):
        return [j - i + 1 for i, j in self.exons]

    def seq_length(self):
        return sum(self.exons_length())

    def test(self):
        """
        Test if the length of coding sequences is 3-divisible
        
        input: Cds object
        output: a boolean
        """
        return self.seq_length()%3==0


    def befile_lines(self):
        lines = []
        for start, end in self.exons:
            lines.append("{0}\t{1}\t{2}\t{3}\t0\t{4}\n".format(self.chromosome, start, end, self.name, self.strand))
        return lines


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--gtf', required=True, type=str,
                        dest="g", metavar="<gtf>",
                        help="The relative name of the .gtf file")
    args = parser.parse_args()
    path = os.getcwd()

    dict_errors_info = {"NotConfirmed": "GTF file: {0} cds have a start or end that could not be confirmed",
                        "EmptyExon": "GTF file: {0} cds have at least one empty track",
                        "LengthError": "GTF file: {0} cds length is not a multiple of 3"}

    dict_errors_cds = {}
    for error in dict_errors_info:
        dict_errors_cds[error] = []

    dict_cds, not_confirmed_tr = build_dict_cds(path, args.g)

    bedfile = open("{0}/{1}.bed".format(path, args.g[:-4]), 'w')
    for tr_id, cds in dict_cds.items():

        if tr_id in not_confirmed_tr:
            dict_errors_cds["NotConfirmed"].append(tr_id)
            continue

        if cds.empty_exon():
            dict_errors_cds["EmptyExon"].append(tr_id)
            continue

        if not cds.test():
            # Exclude CDS that do not fulfill the condition on coding sequence length
            dict_errors_cds["LengthError"].append(tr_id)
            continue

        for exon in cds.befile_lines():
            bedfile.write(exon)

    bedfile.truncate()
    bedfile.close()

    nbr_cds_errors = sum([len(errors) for errors in dict_errors_cds.values()])
    nbr_cds_total = len(dict_cds)

    error_filename = '{0}.gtf_to_bed_errors.txt'.format(args.g[:-4])
    error_header = "{0} errors out of {1} coding sequences ({2:.3f}%)".format(
        nbr_cds_errors, nbr_cds_total, nbr_cds_errors * 100. / nbr_cds_total)
    print(error_header)
    print("Errors written {0}".format(error_filename))

    error_file = open('{0}/{1}'.format(path, error_filename), 'w')
    error_file.write(error_header)
    if nbr_cds_errors > 0:
        for error, list_cds in dict_errors_cds.items():
            nbr_cds = len(list_cds)
            if nbr_cds > 0:
                error_file.write("\n\n" + dict_errors_info[error].format(nbr_cds))
                error_file.write("\n" + "\t".join(list_cds))
    error_file.close()
