"""
Reads the methylation file. Spearate the cancer and normal samples
based on the sample name. 01-10 (cancer), > 10 (normal). Then stores
mean methylation values per gene per cancer in a dictionary for both
cancer sample and normal sample.

Reads ordered gene list from Biomart file
"""
import pandas as pd

class ProcessData(object):

    def __init__(self, param):
        super(ProcessData, self).__init__()
        self.param = param

        self.dataset_stat, self.cancer_all_gene_methylation_dict, self.normal_all_gene_methylation_dict = \
            processData(self.param.getEachCancerMethylationFile())

    def getDatasetStat(self):
        return self.dataset_stat

    def getNameWiseCancerMethylationData(self):
        return self.cancer_all_gene_methylation_dict

    def getNameWiseNormalMethylationData(self):
        return self.normal_all_gene_methylation_dict

    def getOrderedGeneByLocation(self):
        return processGeneLocationData(self.param.getGeneLocationFile())


def processData(methylationFilePathList):
    cancer_sample_stat_dict = dict()
    cancer_all_gene_methylation_dict = dict()
    normal_all_gene_methylation_dict = dict()

    for filePath in methylationFilePathList:
        cancer = (filePath.split('/')[-1]).split('.')[0]
        D = pd.read_table(filePath, sep='\t', header='infer', skiprows=[1])

        cancer_sample_indices = list()
        normal_sample_indices = list()

        sample_names = D.columns.values[1:]  # careful about the index when reading methylation values

        for index, sample in enumerate(sample_names):
            class_label_id = int(sample.split('-')[-1])
            if 0 < class_label_id < 10:
                cancer_sample_indices.append(index)
            else:
                normal_sample_indices.append(index)

        cancer_sample_stat_dict[cancer] = [cancer_sample_indices, normal_sample_indices]

        # read gene wise methylation value
        cancer_methylation_value_dict = dict()
        normal_methylation_value_dict = dict()

        for index, row in D.iterrows():
            gene_name = row['Hybridization REF']
            methylation_values = row[1:]
            cancer_methylation_values = methylation_values[cancer_sample_indices]
            normal_methylation_values = methylation_values[normal_sample_indices]

            cancer_methylation_value_dict[gene_name] = cancer_methylation_values
            normal_methylation_value_dict[gene_name] = normal_methylation_values

        cancer_all_gene_methylation_dict[cancer] = cancer_methylation_value_dict
        normal_all_gene_methylation_dict[cancer] = normal_methylation_value_dict

    return cancer_sample_stat_dict, cancer_all_gene_methylation_dict, normal_all_gene_methylation_dict


def processGeneLocationData(geneLocationFile):
    ordered_list_of_gene = list()

    D = pd.read_table(geneLocationFile, sep='\t', header='infer')
    for index, row in D.iterrows():
        ordered_list_of_gene.append(row['Gene name'])

    return ordered_list_of_gene
