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
def fillMatrix(matrix, sequences):
  ### Ici on laisse la première ligne et la première  #
  #   colonne remplis avec des zéros et on applique   #
  #   l'algorithme de remplissage de score.           #
  #   Match = +4, Mismatch = -4, Indel = -8         ###
  for i in range(1,len(matrix)):
    for j in range(1,matrix[i].size):
      match = matching(sequences[0][i-1],sequences[1][j-1])
      matrix[i][j] = max(0, matrix[i-1][j-1]+match, matrix[i][j-1]+INDEL, matrix[i-1][j]+INDEL)
  return matrix

### match mismatch or indel? return score
def matching(data1, data2):
  if data1 == data2:
    return MATCH
  else:
    return MISMATCH

### Case de la meilleure séquence si plusieurs résultats sont disponibles
def starter(list):
  n = 0
  startPos = [0,0]
  while n <= len(list):
    if sum(startPos) < sum(list[n]):
      startPos = list[n]
    n += 1
  return startPos

### Point de départ de la séquence chevauchée
def startingPos(matrice):
  ### on trouve le(s) point de départ selon le choix.
  start = np.argwhere(matrice == np.amax(matrice)) ##
  #print ("starting position ",start)
  if len(start) > 1:
    if WANTEDSEQUENCES == 'last':
      start = starter(start)
    elif WANTEDSEQUENCES == 'first':
      start = x.argmax(axis=0)
    else:
      raise Exception("Not implemented")
  return start

### sequence pathing
def sequencePath(matrice, pos):
  x = pos[0,0]
  y = pos[0,1]
  path = []
  current = matrice[x][y]
  while ((x > 0) or (y > 0)):
    if (current!=0):
      if (current - MISMATCH == matrice[x-1][y-1]) or (current - MATCH == matrice[x-1][y-1]):
        x -= 1
        y -= 1
        path.append([1,1])
        #print (matrice[x-1][y-1])
      elif (current - INDEL == matrice[x-1][y]):
        x -= 1
        current = matrice[x][y]
        path.append([1,0])
      elif (current - INDEL == matrice[x][y-1]):
        y -= 1
        current = matrice[x][y]
        path.append([0,1])
      else:
        raise Exception("Score match pas avec algo")
      current = matrice[x][y]
    else:
      if x>0:
        x -= 1
        path.append([1,0])
      else:
        y -= 1
        path.append([0,1])
  return path, np.array([x,y])

### alignement des séquences selon un chemin connue
def alignSequences(start, path, seqs):
  seqLength = sum(start) + len(path)
  seq1 = seqs[1].strip("\n")
  seq2 = seqs[0].strip("\n")
  seq1align = ""
  seq2align = ""
  y = 0
  z = 0
  ### initialisation des indels pour le préfixe
  if start[0] > 0:
    length = len(seq2)
    x = start[0]
    while x > 0:
      seq1align = seq1align + "-"
      x -= 1
      seq2align = seq2align + seq2[y]
      y += 1

  elif start[1] > 0:
    length = len(seq1)
    x = start[1]
    while x > 0:
      seq2align = seq2align + "-"
      x -= 1
      seq1align = seq1align + seq1[z]
      z += 1
  else:
    pass

  ### ajout des séquences de chemin
  i = 0
  while i < len(path):
    if path[i][0] == 1:
      seq1align = seq1align + seq1[z]
      z += 1
    else:
      seq1align = seq1align + "-"

    if path[i][1] == 1:
      seq2align = seq2align + seq2[y]
      y += 1
    else: seq2align = seq2align + "-"
    i +=1

  ### complétion de la séquence

  return [seq1align, seq2align]

### main
def main():
  #sequences1 = fetchSequences("test.txt")
  sequences1 = ["GTAGACC","AGCGTAGA"]
  matrice = sequenceMatrix(len(sequences1[0]), len(sequences1[1]))
  matrice = fillMatrix(matrice, sequences1)
  print (matrice)
  start = startingPos(matrice)
  sequences2 = fetchSequences("reads.fq")
  path,end = sequencePath(matrice,start)
  print (path)
  aligned = alignSequences(end, path, sequences1)
  print(aligned[0])
  print(aligned[1])

  return None


if __name__ == "__main__":
  main()