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
gmatrix = DP('P21333.fasta','Q8BTM8.fasta',-5)

gmatrix.makeFmatrix()

cad1, cad2 = gmatrix.backTracking()

print(cad1)
print(cad2)
