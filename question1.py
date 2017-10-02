# -*- coding: utf-8 -*-
#######################################################
### Maximisation du chevauchement & calcul de score ###
### entre deux séquences à partir d'un fichier.     ###
###                                                 ###
### Autheurs: Andre Lalonde, Maude Sabourin         ###
#######################################################

### numpy est utilisé pour les numerical python array #
#   qui sont nettement plus efficaces que les array   #
#   standard. os est utiliser pour le path.         ###
########################
import os
import numpy as np
np.set_printoptions(threshold=np.nan)
########################
### valeur des scores
########################
global MATCH
global MISMATCH
global INDEL
global WANTEDSEQUENCES
MATCH = 4
MISMATCH = -4
INDEL = -8
WANTEDSEQUENCES = 'last' # Les options sont 'first' 'last' et 'all'
########################
### functions
########################
### Lecture des séquences d'un fichier et mise en liste
def fetchSequences(filename):
#  if not os.path.filename(filename):
#    raise Exception("Le fichier n'a pas été trouvé! Veuillez vérifier le pathname.")
  if ".fq" in filename:
    return fastaSequences(filename)
  with open(filename, 'r') as f:
    sequences = f.readlines()
    sequences[0] = sequences[0].strip('\n')
    return sequences

def fastaSequences(filename):
  sequences = []
  with open(filename, 'r') as f:
    for line in f:
      if line[0] in '@+':
        continue
      elif line[0] in 'ACGT':
        sequences.append(line.strip("\n"))
      else:
        continue
  return sequences

### initiation de la matrice
def sequenceMatrix(len1, len2):
  if (len1+len2)*4 > 2**32:
    matrix = np.zeros(shape=(len1+1,len2+1), dtype=np.uint64)
  elif (len1+len2)*4 > 2**16:
    matrix = np.zeros(shape=(len1+1,len2+1), dtype=np.uint32)
  elif (len1+len2)*4 > 2**8:
    matrix = np.zeros(shape=(len1+1,len2+1), dtype=np.uint16)
  else:
    matrix = np.zeros(shape=(len1+1,len2+1), dtype=np.uint8)
  return matrix

### remplissage du tableau
def fillMatrix(matrix, seq1, seq2):
  ### Ici on laisse la première ligne et la première  #
  #   colonne remplis avec des zéros et on applique   #
  #   l'algorithme de remplissage de score.           #
  #   Match = +4, Mismatch = -4, Indel = -8         ###
  for i in range(1,len(matrix)):
    for j in range(1,matrix[i].size):
      match = matching(seq1[i-1],seq2[j-1])
      matrix[i][j] = max(0, matrix[i-1][j-1]+match, matrix[i][j-1]+INDEL, matrix[i-1][j]+INDEL)
  return matrix

### match mismatch or indel? return score
def matching(data1, data2):
  if data1 == data2:
    return MATCH
  else:
    return MISMATCH

### sequence pathing
def sequencePath(matrice, pos, seq1, seq2):
  x = pos[0]
  y = pos[1]
  path = []
  current = matrice[x][y]
  while ((x > 0) and (y > 0)):
    #if (current!=0):
      if ((current - MISMATCH == matrice[x-1][y-1]) and (seq1[x-1] != seq2[y-1])) or ((current - MATCH == matrice[x-1][y-1]) and (seq1[x-1] == seq2[y-1])):
        x -= 1
        y -= 1
        path.append([1,1])
      elif (current - INDEL == matrice[x-1][y]):
        x -= 1
        current = matrice[x][y]
        path.append([1,0])
      elif (current - INDEL == matrice[x][y-1]):
        y -= 1
        current = matrice[x][y]
        path.append([0,1])
      else:
        break
      current = matrice[x][y]
  return path, np.array([x,y])

### alignement des séquences selon un chemin connue
def alignSequences(start, path, seq1, seq2, end, size):

  path.reverse()
  temp = seq2
  seq2 = seq1
  seq1 = temp

  #Création de fonctions pour séparer les étapes
  """Indels de départ"""
  seq1align,seq2align,x,y,z = genIndelStart(start,seq1,seq2)

  """Séquences de chemin"""
  seq1align,seq2align,y,z,i = genSeqPath(seq1align,seq2align,path,y,z,seq1,seq2)

  """Complétion de la séquence"""
  seq1align,seq2align= genIndelEnd(seq1align,seq2align,end,size,seq1,seq2)

  return [seq1align, seq2align], i

def genIndelStart(start, seq1, seq2, ):
  y = 0
  z = 0
  seq1align = ""
  seq2align = ""
  ### initialisation des indels pour le préfixe

  #IMPORTANT : L'erreur est ici! Les indels doivent être dans AGTCA.. et non dans GGGGTT
  if start[0] > 0:
    x = start[0]
    while x > 0:
      seq1align += "-"
      x -= 1
      seq2align += seq2[y]
      y += 1
  if start[1] > 0:
    x = start[1]
    while x > 0:
      seq2align += "-"
      x -= 1
      seq1align += seq1[z]
      z += 1
  else:
    pass
  return seq1align,seq2align,x,y,z

#s1=seqalign1,s2=seqalign2,sq1=seq1,sq2=seq2
def genSeqPath(s1,s2,p,y,z,sq1,sq2):
  i = 0
  while i < len(p):
    if p[i][0] == 1:
      s2 += sq2[y]
      y += 1
    else:
      s2 += "-"

    if p[i][1] == 1:
      s1 += sq1[z]
      z += 1
    else:
      s1 += "-"
    i += 1

  return s1,s2,y,z,i

#s1=seqalign1,s2=seqalign2,sq1=seq1,sq2=seq2
def genIndelEnd(s1,s2,end,size,sq1,sq2):
  size_ligne = size[0] - 1
  size_col = size[1] - 1

  if end[0] < size_ligne:
    x = end[0]
    while x < size_ligne:
      s1 += "-"
      s2 += sq2[x]
      x += 1
  elif end[1] < size_col:
    x = end[1]
    while x < size_col:
      s2 += "-"
      s1 += sq1[x]
      x += 1

  return s1,s2





"""Génère la matrice 20 par 20 des séquences reads"""
def genMatrix2020(sequences):
    # SCORE : M[RX,RY] SI RX SUFFIXE, RY PREFIXE
    #TODO:Remove hardcoded
  matrix = np.zeros(shape=(20, 20))
  for i in range (len(sequences)):
      for j in range (len(sequences)):
        #On ne calcule pas la diagonale
        if not(i==j):
            if (i<j):
                #On calcule le meilleur alignement de la colonne et de la ligne (Rx,Ry) et (Ry,Rx)
                bestLigne,bestCol = genMatrixAlignement(sequences[i], sequences[j], False)
                #TODO: iNVERSE? Pas sure si bestLigne en premier
                if bestCol>bestLigne and bestCol>=80:
                  matrix[i][j]=bestCol
                  print("seq", i , " seq", j, " ", bestCol)
                elif bestLigne>=80:
                  matrix[j][i]=bestLigne
                  print("seq", j, " seq", i, " ", bestLigne)
                #matrix[i][j] = bestCol #Rx suffixe, Ry prefixe
                #matrix[j][i] = bestLigne #Rx prefixe, Ry suffixe
                #print ("Seq ",i, " avec Seq ", j, "| Score:", bestCol, "|  Inverse: ", bestLigne)

  print (matrix)

def reverseSeq(seq):
  seqrev=""
  for s in seq:
    if s=="A":
      seqrev="T"
    elif s=="C":
      seqrev="G"
    elif s=="T":
      seqrev="A"
    elif s=="G":
      seqrev="C"

  seqrev=seqrev[::-1]
  return seqrev

def genMatrixAlignement(seq1, seq2, show):
  alignValue = '1'
  matrice = sequenceMatrix(len(seq1), len(seq2))
  matrice = fillMatrix(matrice, seq1, seq2)

  # Trouve le score et la position totale, en plus de la valeur maximale de la ligne et colonne
  maxLigne, maxCol, score, posFinal = findMax(matrice)
  if show:
      #On affiche que le chevauchement optimal
      print (matrice)
      path, end = sequencePath(matrice, posFinal, seq1, seq2)
      aligned, alignValue = alignSequences(end, path, seq1, seq2, posFinal, matrice.shape)
      print("Sequence 1: " + aligned[0])
      print("Sequence 2: " + aligned[1])
      print("Chevauchement: " + str(alignValue))
      print("Score: " + str(score))

  #Retourne le max ligne colonne pour la matrice 20x20
  return maxLigne, maxCol

def findMax(matrice):
  shape = matrice.shape
  maxLigne,maxCol = 0,0
  posLigne, posCol, posFinal = (0, 0),(0,0),(0,0)
  #Shape a une valeur en trop
  indexLigne = shape[0] - 1
  indexCol = shape[1] - 1

  # Parcourt les colonnes de la dernière ligne et cherche la plus haute valeur
  for i in range(indexCol + 1):
    if matrice[indexLigne][i] > maxLigne:
      maxLigne = matrice[indexLigne][i]
      posLigne = (indexLigne, i)

  # Parcourt les lignes de la dernière col et cherche la plus haute valeur
  for i in range(indexLigne + 1):
    if matrice[i][indexCol] > maxCol:
      maxCol = matrice[i][indexCol]
      posCol = (i, indexCol)

  # Compare les scores pour obtenir le plus élevé
  if (maxLigne > maxCol):
    score = maxLigne
    posFinal = posLigne
  else:
    score = maxCol
    posFinal = posCol
  return maxLigne, maxCol, score, posFinal

#Supprime les espacements des séquences
def stripSeq(seqList):
  for i in range (len(seqList)):
    seqList[i] = seqList[i].strip("\n")
  return seqList

### main
def main():
  while True:
    #!IMPORTANT: Pour tester le bug, faire 3
      res = input("1 pour comparer deux séquences, 2 pour voir la matrice d'adjacence 20x20, 3 POUR VOIR LE BUG \n")
      if res == "1":
          file = input("Veuillez entrer le nom du fichier avec son extension (ex: test.txt) \n")
          sequences1 = fetchSequences(file)
          sequences1 = stripSeq(sequences1)
          genMatrixAlignement(sequences1[0], sequences1[1], True)
      elif res == "2":
          sequences2 = fetchSequences("reads.fq")
          sequences2 = stripSeq(sequences2)
          genMatrix2020(sequences2)
      else:
          break

  sequences3 = fetchSequences("test.txt")
  sequences4 = fetchSequences("test2.txt")
  sequences5 = fetchSequences("test3.txt")
  sequences6 = fetchSequences("test5.txt")
  sequences7 = ["GTAGACC", "AGCGTAGA"]
  sequences6 = stripSeq(sequences7)
  genMatrixAlignement(sequences6[0],sequences6[1],True)


  return None


if __name__ == "__main__":
  main()