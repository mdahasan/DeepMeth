"""
First finds the common gene that exists in both ordered gene list and the
methylation file for each cancer. The based on the common cancer gene list
gene by sample matrix is build for each cancer that contains methylation
values per gene per sample.
"""

import numpy as np

class SampleMaking(object):
    
    def __init__(self, data):
        super(SampleMaking, self).__init__()

        self.data = data
        self.all_cancer_gene_methylation_dict = self.data.getNameWiseCancerMethylationData()
        self.ordered_gene_list = self.data.getOrderedGeneByLocation()
        self.by_cancer_common_gene_list_dict = findCommonGenes(self.all_cancer_gene_methylation_dict, self.ordered_gene_list)
        self.per_caner_methylation_matrix_dict = \
            buildMethylationMatrixPerCancer(self.all_cancer_gene_methylation_dict, self.by_cancer_common_gene_list_dict)

    def getCommonGeneList(self):
        return self.by_cancer_common_gene_list_dict

    def getPerCanerMethylationMatrix(self):
        return self.per_caner_methylation_matrix_dict


def findCommonGenes(cancer_gene_methylation_dict, gene_list):
    cancer_common_gene_dict = dict()

    for cancer_name in cancer_gene_methylation_dict:
        common_gene_list = list()
        gene_methylation_dict = cancer_gene_methylation_dict[cancer_name]
        for gene in gene_list:
            if gene in gene_methylation_dict:
                common_gene_list.append(gene)

        cancer_common_gene_dict[cancer_name] = common_gene_list

    return cancer_common_gene_dict

def buildMethylationMatrixPerCancer(all_cancer_gene_methylation_dict, by_cancer_common_gene_dict):

    by_cancer_methylation_sample_matrix_dict = dict()

    for cancer_name in all_cancer_gene_methylation_dict:
        gene_methylation_dict = all_cancer_gene_methylation_dict[cancer_name]
        cancer_common_gene = by_cancer_common_gene_dict[cancer_name]

        gene_by_sample_list = list()

        for gene in cancer_common_gene:

            # potential feature selection: to get rid of some gene that
            # shows less variation. Or some other kind of feature selection method
            # to reduce the number of genes

            gene_by_sample_list.append(gene_methylation_dict[gene])

        by_cancer_methylation_sample_matrix_dict[cancer_name] = np.array(gene_by_sample_list)

    return by_cancer_methylation_sample_matrix_dict
