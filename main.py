#from DotMatrix.dotMatrix import dotMatrix
from NeedlemanWunsch.alineamiento import Alineacion as DP

"""
Compilacion del Algoritmo Dot Matrix
gmatrix=dotMatrix()
#gmatrix.loadData('P21333.fasta','Q8BTM8.fasta')
gmatrix.loadData('P0C6X7.fasta','P0DTC1.fasta')
gmatrix.setup(20,75)

newmatrix=gmatrix.generateMatrix()
gmatrix.plotMatrix()

"""


#Compilacion del Algoritmo Needleman-Wunsvh
gmatrix = DP('P21333.fasta','Q8BTM8.fasta',-1)

gmatrix.makeFmatrix()

vecA,vecB = gmatrix.traceback()
print(vecA,vecB)

