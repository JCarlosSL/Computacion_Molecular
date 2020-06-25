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
        """
        self.matrixS = {'A':{'A':2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'C':{'A':-2,'C':2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'D':{'A':-2,'C':-2,'D':2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'E':{'A':-2,'C':-2,'D':-2,'E':2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'F':{'A':-2,'C':-2,'D':-2,'E':-2,'F':2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'G':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'H':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'I':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'K':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'L':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'M':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'N':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'P':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'Q':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'R':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'S':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'T':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':2,'V':-2,'W':-2,'Y':-2},
                        'V':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':2,'W':-2,'Y':-2},
                        'W':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':2,'Y':-2},
                        'Y':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':2}}

        """
        self.cadenaA = _cadenaA
        self.cadenaB = _cadenaB
        self.matrixF = np.zeros((len(self.cadenaB)+1, len(self.cadenaA)+1))
        self.gap = 0
        
    def updateMatrixS(self,match,missmatch):
        for key in self.matrixS:
            for key2 in self.matrixS[key]:
                if(key == key2):
                    self.matrixS[key][key2] = match
                else:
                    self.matrixS[key][key2] = missmatch

    def _local(self, _match=0, _missmatch=0, _gap=0):
        self.updateMatrixS(_match,_missmatch)
        self.gap = _gap
        self.makeFLocalMatrix()
        return self.localBackTracking(_match, _missmatch, _gap)

    def makeFLocalMatrix(self):
        max_valor = 0
        list_posiciones = []
        filas, columnas = len(self.cadenaB), len(self.cadenaA)
        
        for i in range(filas):
            for j in range(columnas):
                self.matrixF[i+1][j+1] = max(
                    0, 
                    self.matrixF[i+1][j]+self.gap,
                    self.matrixF[i][j+1] + self.gap, 
                    self.matrixF[i][j] + self.matrixS[self.cadenaB[i]][self.cadenaA[j]])
                if self.matrixF[i+1][j+1] > max_valor:
                    max_valor = self.matrixF[i+1][j+1]
                    list_posiciones = []
                    list_posiciones.append([i,j])
                elif self.matrixF[i+1][j+1] == max_valor:
                    list_posiciones.append([i,j])
        return max_valor, list_posiciones
    
    def score(self, match, missmatch, gap, cadA, cadB):
       score_=0
       for i in range(len(cadA)):
           if cadA[i] == cadB[i]:
               score_=score_ + match
           elif cadA[i]=='-' or cadB[i]=='-':
               score_=score_ + gap
           elif cadA[i] != cadB[i]:
               score_=score_ + missmatch
       return score_
    
    def localBackTracking(self, match, missmatch, gap):
        _max_valor, _list_posiciones = self.makeFLocalMatrix()
        print(self.matrixF)
        num = len(_list_posiciones)/2 + 1
        num = int(num)
        newcadA = ""
        newcadB = ""
        for i in range (num):
            newcadA = self.cadenaA[_list_posiciones[i][1]]
            newcadB = self.cadenaB[_list_posiciones[i][0]] 
            m=_list_posiciones[i][0]
            n=_list_posiciones[i][1]
            while( m >= 0 and n >= 0 and self.matrixF[m][n] != 0 ):
                newcadA = self.cadenaA[n-1] + newcadA
                newcadB = self.cadenaB[m-1] + newcadB
                m-=1
                n-=1
            print(newcadA)
            print(newcadB)
            print(self.score(match, missmatch, gap, newcadA,newcadB))
        
alinear = Alligment('ATACTGGG','TGACTGAG')
alinear._local(1,-1,-2)

#alinear = Alligment('AAG','AGC')
#alinear._local(2,-7,-5)
