"""
Analyzes the parameters. Generates list of the name of the
cancers used in the experiment. Also process the file path
for the methylation file directory and the ordered gene
file directory
"""

import os

class ParameterProcessing(object):

    def __init__(self, methylationFileDir, geneLocationFile):
        super(ParameterProcessing, self).__init__()
        self.methylationDir = methylationFileDir
        self.geneLocationFile = geneLocationFile

    def getListofCancerNames(self):
        CancerNameList = list()
        fileNameList = os.listdir(self.methylationDir)

        for file in fileNameList:
            cancer = file.split('.')[0]
            if len(cancer) > 1:
                CancerNameList.append(cancer)

        return CancerNameList

    def getEachCancerMethylationFile(self):
        filePathList = list()
        methyDirPath = os.path.abspath(self.methylationDir)
        filePaths = os.listdir(self.methylationDir)

        for file in filePaths:
            if 'meth' in file:
                filePathList.append(os.path.join(methyDirPath, file))

        return filePathList

    def getGeneLocationFile(self):
        return self.geneLocationFile
