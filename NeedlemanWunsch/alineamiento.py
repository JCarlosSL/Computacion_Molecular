from Bio import SeqIO
import numpy as np

def datacleaning(path,dtype='fasta'):
    sequences = SeqIO.parse(path,dtype)
    for record in sequences:
        data = str(record.seq.upper())
    return data


class Alligment:
    def __init__(self,_cadenaA,_cadenaB):
        
        self.matrixS = {'A' : {'A' : 2, 'C' : -7, 'G' : -5 , 'T' : -7},
                'C' : {'A' : -7, 'C' : 2, 'G' : -7, 'T' : -5},
                'G' : {'A' : -5, 'C' : -7, 'G' : 2, 'T' : -7},
                'T' : {'A' : -7, 'C' : -5, 'G' : -7, 'T' : 2}}
        
        self.cadenaA ="ATAGA" #datacleaning(_cadenaA)
        self.cadenaB ="TTCAG" #datacleaning(_cadenaB)
        self.matrixF = np.zeros((len(self.cadenaB)+1, len(self.cadenaA)+1))
        self.match = 1
        self.missmatch = -1
        self.gap = -5

    def _local(self, _match=1, _missmatch=-1, _gap=-1):
        self.match = _match
        self.missmatch = _missmatch
        self.gap = _gap
        self.makeFLocalMatrix()
        return self.localBackTracking()

    def _global(self, _match=1, _missmatch=-1, _gap=-1):
        self.match = _match
        self.missmatch = _missmatch
        self.gap = _gap
        self.makeFGlobalMatrix()
        return self.globaltraceback()

    def makeFLocalMatrix(self):
        filas, columnas = len(self.cadenaB), len(self.cadenaA)
        for i in range(filas):
            for j in range(columnas):
                self.matrixF[i+1][j+1] = max(0, self.matrixF[i+1][j]+self.gap,
                    self.matrixF[i][j+1] + self.gap,
                    self.matrixF[i][j] + self.matrixS[self.cadenaB[i]][self.cadenaA[j]])
        print(self.matrixF)


    def makeFGlobalMatrix(self):
        filas, columnas = len(self.cadenaB), len(self.cadenaA)
        self.matrixF[0][0] = 0

        for i in range(1,filas+1):
            self.matrixF[i][0] = self.matrixF[i-1][0] + self.gap
        for j in range(1,columnas+1):
            self.matrixF[0][j] = self.matrixF[0][j-1] + self.gap

        for i in range(filas):
            for j in range(columnas):
                self.matrixF[i+1][j+1] = max(self.matrixF[i+1][j] + self.gap,
                    self.matrixF[i][j+1] + self.gap,
                    self.matrixF[i][j] + self.matrixS[self.cadenaB[i]][self.cadenaA[j]])

        print(self.matrixF)
    def localBackTracking(self):
        newcadA = ""
        newcadB = ""

    def globaltraceback(self):
        vecA = []
        vecB = []
        newcadA = ""
        newcadB = ""
        n = len(self.cadenaA) - 1
        m = len(self.cadenaB) - 1
        self.globalBackTracking(m,n,newcadA,newcadB,vecA,vecB)
        return vecA,vecB

    def  globalBackTracking(self,m,n,newcadA,newcadB,vecA,vecB):
        count = 0
        temporalA = ""
        temporalB = ""
        while( m >= 0 and n >= 0):
            score = self.matrixF[m+1,n+1]            
            scorediag = self.matrixF[m,n]
            scoreup = self.matrixF[m+1,n]
            scoreleft = self.matrixF[m,n+1]
            
            if score == scorediag + self.matrixS[self.cadenaB[m]][self.cadenaA[n]]:
                count = 1
                newcadA = self.cadenaA[n] + newcadA
                newcadB = self.cadenaB[m] + newcadB
                m-=1
                n-=1
            if score == scoreleft + self.gap:
                if count==1:
                    self.globalBackTracking(m,n+1,"-"+temporalA,
                            self.cadenaB[m+1]+temporalB,vecA,vecB)
                else:
                    newcadB = self.cadenaB[m] + newcadB
                    newcadA = "-"+newcadA
                    m-=1
            if score == scoreup + self.gap:
                if count == 1:
                    self.globalBackTracking(m+1,n,
                            self.cadenaA[n+1]+temporalA,
                            "-"+temporalB,vecA,vecB)
                else:
                    newcadB = "-" +newcadB
                    newcadA = self.cadenaA[n] + newcadA
                    n-=1

            count = 0
            temporalA = newcadA
            temporalB = newcadB
        
        while(m >= 0):
            newcadB = self.cadenaB[m] + newcadB
            newcadA = "-" + newcadA
            m-=1
        
        while(n >= 0):
            newcadB = "-" + newcadB
            newcadA = self.cadenaA[n] + newcadA
            n-=1
        
        vecA.append(newcadA)
        vecB.append(newcadB)



