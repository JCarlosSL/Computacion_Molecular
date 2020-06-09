from dotMatrix import dotMatrix

gmatrix=dotMatrix()
gmatrix.loadData('P21333.fasta','Q8BTM8.fasta')
#gmatrix.loadData('P0C6X7.fasta','P0DTC1.fasta')
gmatrix.setup(10,25)

newmatrix=gmatrix.generateMatrix()
gmatrix.plotMatrix()
