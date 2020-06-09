from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def datacleaning(path='',dtype='fasta'):
    sequences=SeqIO.parse(path,dtype)
    for record in sequences:
        data=str(record.seq.upper())
    return data

def split(str1,lim):
    v=[]
    for i in range(0,len(str1),lim):
        v.append(str1[i:i+lim])
    return np.array(v)

class dotMatrix:

    def __init__(self):
        self.seq_1 = None
        self.seq_2 = None
        self.windowS = 1
        self.threshold = 100
        self.percentU = 100

    def loadData(self,data1='',data2=''):
        self.seq_1 = datacleaning(data1)
        self.seq_2 = datacleaning(data2)

    def setup(self,_windowS=10,_threshold=25):
        self.windowS = _windowS
        self.threshold = _threshold
        
        mod1 = len(self.seq_1)%_windowS
        mod2 = len(self.seq_2)%_windowS
        
        for i in range(_windowS-mod1):
            self.seq_1 += '0'
        for j in range(_windowS-mod2):
            self.seq_2 += '1'

        self.seq_1 = split(self.seq_1,_windowS)
        self.seq_2 = split(self.seq_2,_windowS)

        self.percentU=100/_windowS

    def thresHold(self,sub1,sub2):
        v=0.0
        for i in range(self.windowS):
            if sub1[i] == sub2[i]:
                v+=self.percentU
            if v >= self.threshold:
                break

        if v >= self.threshold:
            return 1
        else:
            return 0


    def generateMatrix(self):

        self.newMatrix=[ [] for i in range(len(self.seq_1))]

        i=0
        for sub1 in self.seq_1:
            for sub2 in self.seq_2:
                self.newMatrix[i].append(self.thresHold(sub1,sub2))
            i+=1
        return np.array(self.newMatrix)

    def plotMatrix(self):
        arr1 = []
        arr2 = []
        rows,columns = len(self.seq_1),len(self.seq_2)
        for i in range(rows):
            for j in range(columns):
                if self.newMatrix[i][j] == 1:
                    arr1.append(i)
                    arr2.append(j)

        plt.plot(arr1,arr2,'b.',markersize=0.5)
        plt.show()

#gmatrix=dotMatrix()
#gmatrix.loadData('P21333.fasta','Q8BTM8.fasta')
#gmatrix.loadData('P0C6X7.fasta','P0DTC1.fasta')
#gmatrix.setup(10,25)

#newmatrix=gmatrix.generateMatrix()
#gmatrix.plotMatrix()

