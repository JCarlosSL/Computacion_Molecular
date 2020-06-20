from Bio import SeqIO
import numpy as np

def datacleaning(path,dtype='fasta'):
    sequences = SeqIO.parse(path,dtype)
    for record in sequences:
        data = str(record.seq.upper())
    return data


class Alineacion:
    def __init__(self,_cadenaA,_cadenaB, _d):
        
        self.matrixS = {'A' : {'A' : 2, 'C' : -7, 'G' : -5 , 'T' : -7},
                        'C' : {'A' : -7, 'C' : 2, 'G' : -7, 'T' : -5},
                        'G' : {'A' : -5, 'C' : -7, 'G' : 2, 'T' : -7},
                        'T' : {'A' : -7, 'C' : -5, 'G' : -7, 'T' : 2}}
        
        self.cadenaA ="AAG" #datacleaning(_cadenaA)
        self.cadenaB ="AGC" #datacleaning(_cadenaB)
        self.d = _d
        self.matrixF = np.zeros((len(self.cadenaB)+1, len(self.cadenaA)+1))

    def makeFmatrix(self):
        filas, columnas = len(self.cadenaB), len(self.cadenaA)
        self.matrixF[0][0] = 0

        for i in range(1,filas+1):
                self.matrixF[i][0] = self.matrixF[i-1][0] + self.d
        for j in range(1,columnas+1):
          self.matrixF[0][j] = self.matrixF[0][j-1] + self.d

        for i in range(filas):
          for j in range(columnas):
                self.matrixF[i+1][j+1] = max(self.matrixF[i+1][j]+self.d,
                        self.matrixF[i][j+1]+self.d,
                        self.matrixF[i][j]+
                        self.matrixS[self.cadenaB[i]][self.cadenaA[j]])
        print(self.matrixF)
    def backTracking(self):
        newcadB = ""
        newcadA = ""
        m = len(self.cadenaB)-1
        n = len(self.cadenaA)-1

        while( m >= 0 and n >= 0):
            score = self.matrixF[m+1,n+1]
            
            scorediag = self.matrixF[m,n]
            scoreup = self.matrixF[m+1,n]
            scoreleft = self.matrixF[m,n+1]
            
            if score == scorediag + self.matrixS[self.cadenaB[m]][self.cadenaA[n]]:
                newcadA = self.cadenaA[n] + newcadA
                newcadB = self.cadenaB[m] + newcadB
                m-=1
                n-=1
            elif score == scoreleft + self.d:
                newcadB = self.cadenaB[m] + newcadB
                newcadA = "-"+newcadA
                m-=1
            else:
                #score == scoreup + d
                newcadB = "-" +newcadB
                newcadA = self.cadenaA[n] + newcadA
                n-=1
        
        print(m,n)
        
        while(m >= 0):
            newcadB = self.cadenaB[m] + newcadB
            newcadA = "-" + newcadA
            m-=1
        
        while(n >= 0):
            newcadB = "-" + newcadB
            newcadA = self.cadenaA[n] + newcadA
            n-=1
        
        return newcadA,newcadB