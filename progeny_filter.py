#Daniel Ence
#May 1, 2018
#This script filters a set of diploid variant calls
#into two groups


import argparse
from cyvcf2 import VCF


def main(args):

    tmp_susc_handle = open(args.susc_samples)
    susc_list = []
    for line in tmp_susc_handle:
        line = line.strip()
        susc_list.append(line)

    tmp_rest_handle = open(args.rest_samples)
    rest_list = []
    for line in tmp_rest_handle:
        line = line.strip()
        rest_list.append(line)

    tmp_progeny_filter = progeny_filter(args.vcf, susc_list, args.max_susc_maf, rest_list, args.min_rest_maf)
    tmp_progeny_filter.filter_susc_samples()
    tmp_progeny_filter.filter_other_samples()
    #tmp_progeny_filter.print_filtered_rest_results()
    variant_list = tmp_progeny_filter.intersect_filtered_sets()
    for variant in variant_list:
        print str(variant).strip()


class progeny_filter(object):

    def __init__(self):
        self.__vcf_file = ""

    def __init__(self,filename, susc_list, max_susc_maf, rest_list, min_rest_maf):
        self.__vcf_file = filename

        self.__susc_list = susc_list
        self.__filtered_susc_variants = []
        self.max_susc_maf = max_susc_maf
        self.__rest_list = rest_list
        self.min_rest_maf = min_rest_maf
        self.__filtered_rest_variants = []

    #maybe add setters for the lists and the MAF's
    def intersect_with_variant_list(self, ext_variant_list):
        tmp_list = []

        for curr_variant in intersect_filtered_sets():
            for ext_variant in ext_variant_list:
                if(self.variants_overlap(curr_variant, ext_variant) ==1):
                    tmp_list.append(curr_variant)
        return tmp_list

    def variants_overlap(self, variant1, variant2):

        if(variant1.CHROM == variant2.CHROM):
            if(variant1.start >= variant2.start and variant1.end <= variant2.end):
                return 1
        return 0

    def get_filtered_rest_results(self):
        return self.__filtered_rest_variants

    def print_filtered_rest_results(self):
        for variant in self.__filtered_rest_variants:
            print str(variant).strip()

    def intersect_filtered_sets(self):

        passed_variants = []
        for curr_variant in self.__filtered_rest_variants:
            for other_variant in self.__filtered_susc_variants:
                if(curr_variant.CHROM == other_variant.CHROM):
                    if(curr_variant.start >= other_variant.start and curr_variant.end <= other_variant.end):
                        passed_variants.append(curr_variant)

        return passed_variants

    def filter_other_samples(self):
        group_members_list = self.setup_group_members_list(self.__vcf_file, self.__rest_list)

        tmp_VCF = VCF(self.__vcf_file)

        for variant in tmp_VCF:
            if(len(variant.ALT) == 1 and self.passed_MAF_filter(variant, group_members_list, self.min_rest_maf, 0) == 1):
                self.__filtered_rest_variants.append(variant)

    def filter_susc_samples(self):
        group_members_list = self.setup_group_members_list(self.__vcf_file, self.__susc_list)

        tmp_VCF = VCF(self.__vcf_file)

        for variant in tmp_VCF:
            #filter down to just the biallelic sites
            if(len(variant.ALT) == 1 and self.passed_MAF_filter(variant, group_members_list, self.max_susc_maf, 1) == 1):
                #print variant
                self.__filtered_susc_variants.append(variant)


    #below_or_above is a boolean variable indicating whethe the comparison should
    #be below (1) or above the maf variable
    def passed_MAF_filter(self, variant, members_list, MAF, below_or_above):

        group_count = 0
        alt_allele_count = 0

        for i in range(len(variant.gt_types)):
            if(members_list[i] == 1):
                if(variant.gt_types[i] ==  0):
                    group_count = group_count + 1
                elif(variant.gt_types[i] == 1):
                    group_count = group_count + 1
                    alt_allele_count = alt_allele_count + 1
                elif(variant.gt_types[i] == 3):
                    group_count = group_count + 1
                    alt_allele_count = alt_allele_count + 2
                #if gt_types[i] == 2, then it is missing/unknown and we don't count it

        if(group_count == 0):
            tmp_maf = 0
        else:
            tmp_maf = float(alt_allele_count) / float(group_count * 2)

        if(below_or_above == 1):
            return tmp_maf <= MAF
        else:
            return tmp_maf >= MAF

    def setup_group_members_list(self, filename, select_list):

        tmp_VCF = VCF(filename)

        group_members_hash = {}
        for tmp_member in select_list:
            group_members_hash.setdefault(tmp_member,1)

        group_members_list = []

        for sample in tmp_VCF.samples:
            if(group_members_hash.get(sample) == 1):
                group_members_list.append(1)
            else:
                group_members_list.append(0)

        return group_members_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf",type=str,help="vcf file of diploid")
    parser.add_argument("-susc_samples",type=str,help="file containing a list of susceptible samples")
    parser.add_argument("-max_susc_maf",type=float,help="maximum allele frequency for susceptible samples")
    parser.add_argument("-rest_samples", type=str, help="file containing a list of resistant samples")
    parser.add_argument("-min_rest_maf", type=float, help="minimum allele frequency for resistant samples")

    args = parser.parse_args()
    main(args)