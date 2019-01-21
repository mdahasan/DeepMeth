"""
This tool classify the cancer type based on the methylation data
Methylation data is transformed into images per sample and the images
is classifier using convolutional neural network
author: Md Abid Hasan
University of California Riverside
Contact: mhasa006@ucr.edu
"""

import sys
import time

from AnalyzeParameters import *
from DataProcessing import *
from TransformSample import *
from SampleImage import *


def main():
    methylationFilesDir = sys.argv[1]  # directory containing all methylation files
    geneLocationFile = sys.argv[2]  # file containing gene name and location in genome
    experimentName = sys.argv[3]  # name of the experiment

    param = ParameterProcessing(methylationFilesDir, geneLocationFile)
    param.getEachCancerMethylationFile()
    data = ProcessData(param)
    sample = SampleMaking(data)
    image = BuildIMage(sample)

    image.getSampleImagePerCancer()


if __name__ == '__main__':
    main()
