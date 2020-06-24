#from DotMatrix.dotMatrix import dotMatrix
from alineamiento import Alligment

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
alinear = Alligment('P21333.fasta','Q8BTM8.fasta')

seq1, seq2 = alinear.alLocal(1,-1,-2)
seq1, seq2 = alinear.alGlobal(1,-1,-2)

#print(seq1)
#print(seq2)
