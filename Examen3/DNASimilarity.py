from Bio import SeqIO
import numpy as np
import cv2 as cv
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
import plotly.figure_factory as ff

def datacleaning(path,dtype='fasta'):
    sequences = SeqIO.parse(path,dtype)
    key = None
    for record in sequences:
        data = str(record.seq.upper())
        key = record.id
        #key = record.description
    return data,key

def convertimage(data):
    sequences = []
    j=0
    for i in range(70,len(data)+70,70): #add 70
        sequences.append(data[j:i])
        j+=70
    hight = len(sequences)
    width = 69
    newimg = np.zeros((hight,width))
    NT = {'A':{'A':1,'G':17,'C':34,'T':51},
          'G':{'A':68,'G':85,'C':102,'T':119},
          'C':{'A':136,'G':153,'C':170,'T':187},
          'T':{'A':204,'G':221,'C':238,'T':255}}

    for i in range(hight):
        for j in range(1,len(sequences[i])):
            newimg[i,j-1]=NT[sequences[i][j-1]][sequences[i][j]]
    return newimg

def histogram(matrix):
    height,width = matrix.shape
    values = [1,17,34,51,68,85,102,119,136,153,170,187,204,221,238,255]
    hist = {0:0,1:0,17:0,34:0,51:0,68:0,85:0,102:0,119:0,136:0,
            153:0,170:0,187:0,204:0,221:0,238:0,255:0}
    for x in range(height):
        for y in range(width):
            f = matrix[x,y]
            hist[f] += 1
    p={}
    mn=height*width
    for key in hist:
        p[key]=hist[key]/mn
    return hist,p

def _mu(p):
    #valor medio dado
    """nivel de intensidad promedio de la imagen"""
    mu = 0
    for key in p:
        mu += key*p[key]
    return mu

def _sigma(mu,p):
    #varianza
    """variación de la densidad alrededor de la media"""
    sigma2 = 0
    for key in p:
        sigma2 +=  np.power((key-mu),2)*p[key]
    return sigma2

def mu_3(sigma,mu,p):
    #skewness
    """cero si el histograma es simetrico alrededor de la media"""
    mu3 = 0
    for key in p:
        mu3 += np.power((key - mu),3)*p[key]
    return mu3/(sigma*np.sqrt(sigma))

def mu_4(sigma,mu,p):
    #Kurtosis
    """medida de planeidad del histograma"""
    mu4 = 0
    for key in p:
        mu4 += np.power((key - mu),4)*p[key]
    return (mu4/np.power(sigma,2))-3

def Energy(p):
    #energia
    """medida de la uniformidad del histograma"""
    s = 0
    for key in p:
        s += np.power(p[key],2)
    return s

def Entropy(p):
    #Entropia
    """estimacion de aletoriedad, mide su textura"""
    s = 0
    for key in p:
        s += p[key]*np.log2(p[key])
    return -s

def procesamiento(data):
    image = convertimage(data)
    hist,prob = histogram(image)
    #plt.hist(image.ravel(),255,[0,255])
    #plt.show()
    mu = _mu(prob)
    sigma = _sigma(mu,prob)
    mu3 = mu_3(sigma,mu,prob)
    mu4 = mu_4(sigma,mu,prob)
    E = Energy(prob)
    H = Entropy(prob)
    return [mu3,mu4,E,H]

if __name__=='__main__':
    v = []
    names = []

    for i in range(13):
        data,key = datacleaning('sequence_'+str(i)+'.fasta')
        v.append(procesamiento(data))
        names.append(key)

    y = pdist(v,'euclidean')
    x = linkage(y,'average','euclidean')

    fig = ff.create_dendrogram(x,orientation='left',labels=names)
    fig.update_layout(width=800,height=800)
    fig.show()
"""
data,key = datacleaning('sequence_10.fasta')
print(key,procesamiento(data))
image = convertimage(data)
cv.imwrite('image.png',image)
plt.imshow(image,'gray')
plt.show()
"""
