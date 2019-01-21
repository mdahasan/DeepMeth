"""
Process methlyation matrix into image for each sample. Each column in the methylation
matrix represent each sample methylation values across all genes.
Some feature (gene) reduction may be applied to reduce the number of gene and retain
only important genes.
List of gene methylation values are converted into square matrix (image) of size n X n
0 padded at the end of the square matrix (image) if necessary

returns a dictionary of per sample image per cancer
"""

import math
import numpy as np

class BuildIMage(object):

    def __init__(self, sample):
        super(BuildIMage, self).__init__()

        self.sample = sample
        self.by_cancer_reduced_gene_sample_matrix_dict = gene_reduction(self.sample.getPerCanerMethylationMatrix())

    def getSampleImagePerCancer(self):
        return processMethylationMatrixToImage(self.by_cancer_reduced_gene_sample_matrix_dict)

# reduce the size of the genes with some feature reduction
def gene_reduction(methylation_matrix_dict):
    # apply row wise (each gene) feature reduction

    reduced_methylation_matrix_dict = methylation_matrix_dict

    return reduced_methylation_matrix_dict


# returns a list of images of samples per cancer in a dictionary
def processMethylationMatrixToImage(methylation_matrix_dict):
    by_cancer_sample_image_dict = dict()

    for cancer_name in methylation_matrix_dict:
        methylation_matrix = methylation_matrix_dict[cancer_name]
        genes, samples = methylation_matrix.shape
        sample_image_list = list()
        for sample in range(samples):
            methylation_value_per_sample = methylation_matrix[:, sample]
            image_matrix = sampleToImageMatrix(methylation_value_per_sample.tolist())
            sample_image_list.append(image_matrix)

        by_cancer_sample_image_dict[cancer_name] = sample_image_list

    return by_cancer_sample_image_dict


# return a nXn matrix of methylation value list
def sampleToImageMatrix(methylation_value_list):
    # first pad zeros at the end to make a square matrix if necessary
    dim = int(math.ceil(math.sqrt(len(methylation_value_list))))
    methylation_value_list.extend([0] * ((dim * dim) - len(methylation_value_list)))

    return np.reshape(np.array(methylation_value_list), (dim, dim))
