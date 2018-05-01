#Daniel Ence
#May 1, 2018
#This script identifies heterzygyous variants in a set of 
#megagametophyte (haploid) 


import os
import sys
from numpy.random import binomial
import argparse
from cyvcf2 import VCF
from scipy.stats import binom

def main(args):
    tmp_get_hets = get_hets(args.vcf)
    print tmp_get_hets.get_proper_hets()


class get_hets(object):

    def __init__(self):
        #bare constructor.
        #need to do a lot of
        self.__vcf_file = ""
        self.__min_depth = 0
        self.__max_depth = 0
        self.__missing_var = 0
        self.__missing_indv = 0
        self.__num_samples = 0


    def __init__(self,fileName):
        self.__vcf_file = VCF(fileName)
        self.__set_num_samples(fileName)
        self.__min_depth = 1
        self.__max_depth = 100
        self.__missing_var = 0.2
        self.__missing_indv = 0.2

    #Loaders, setters
    def load_vcf_file(self, fileName):
        self.__vcf_file = VCF(fileName)
        self.__set_num_samples(fileName)

    def __set_num_samples(self, fileName):
        tmp_vcf_file = VCF(fileName)

        for variant in tmp_vcf_file:
            self.__num_samples = len(variant.gt_types)
            break

    def set_min_depth(self, depth):
        self.min_depth = depth
    def set_max_depth(self, depth):
        self.max_depth = depth
    def set_missing_var(self, missing):
        self.missing_var = missing
    def set_missing_indv(self, missing):
        self.missing_indv = missing

    def get_proper_hets(self):

        #just get the binomial added in, add filters later
        tmp_interval = binom.interval(0.95, self.__num_samples, 0.5)
        interval = [tmp_interval[0] * self.__num_samples, tmp_interval[1] * self.__num_samples]


        proper_variants = []
        for variant in self.__vcf_file:
            if(len(variant.ALT) == 1 and
                    variant.INFO.get('AF') >= interval[0] and
                    variant.INFO.get('AF') <= interval[1]):
                proper_variants.append(variant)

        return proper_variants


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf",type=str,help="vcf file of haploid variant calls")
    parser.add_argument("-min_depth", type=int, default=1, help="minimum mean depth per variant to be included")
    parser.add_argument("-max_depth", type=int, default=100, help="maximuum mean depth per variant to be included")
    parser.add_argument("-missing_var",type=float, default=0.2, help="maximum percent missing data per variant")
    parser.add_argument("-missing_indv",type=float, default=0.2, help="maximum percent missing data per per sample")
    args = parser.parse_args()
    main(args)