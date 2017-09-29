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
  while n < len(list):
    if sum(startPos) < sum(list[n]):
      startPos = list[n]
    n += 1
  return startPos

### Point de départ de la séquence chevauchée
def startingPos(matrice):
  ### on trouve le(s) point de départ selon le choix.
  start = np.argwhere(matrice == np.amax(matrice))
  size = matrice.shape
  final_strt=0
  if len(start) > 1:
    if WANTEDSEQUENCES == 'last':
      final_strt = starter(start)
    elif WANTEDSEQUENCES == 'first':
      final_strt = x.argmax(axis=0)
    else:
      raise Exception("Not implemented")
  if start[0][0] == size[0]-1:
    suffixeprefix = False
  elif start[0][1] == size[1]-1:
    suffixeprefix = True
  else:
    raise Exception ("pas un suffixe prefixe ni un prefixe suffixe")
  return (suffixeprefix,start)

### sequence pathing
def sequencePath(matrice, pos, seq1, seq2):
  x = pos[0,0]
  y = pos[0,1]
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
  print (path)
  return path, np.array([x,y])

### alignement des séquences selon un chemin connue
def alignSequences(start, path, seqs, end, size):
  path.reverse()
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
      seq1align += "-"
      x -= 1
      seq2align += seq2[y]
      y += 1
  elif start[1] > 0:
    length = len(seq1)
    x = start[1]
    while x > 0:
      seq2align += "-"
      x -= 1
      seq1align += seq1[z]
      z += 1
  else:
    pass

  ### ajout des séquences de chemin
  i = 0
  while i < len(path):
    if path[i][0] == 1:
      seq2align += seq2[y]
      y += 1
    else:
      seq2align += "-"

    if path[i][1] == 1:
      seq1align += seq1[z]
      z += 1
    else: seq1align = seq1align + "-"
    i +=1
  ### complétion de la séquence
  #TODO: On gère le cas PRÉFIXE/SUFFIXE mais doit modifier la size si SUFFIZE/PREFIXE car matrice non symétrique
  size_int = size[0] - 1
  #TODO: Trouver une condition moins dégueux
  if end[0][0]<size_int:
    x = end[0][0]
    while x<size_int:
      seq1align += "-"
      seq2align += seq2[x]
      x += 1
  elif end[0][1]<size_int :
    x = end[0][1]
    while x<size_int:
      seq2align += "-"
      seq1align += seq1[x]
      x+=1

  return [seq1align, seq2align], i
#derniere colonne : suffixe prefixe

### main
def main():
  cheval = '1'
  sequences1 = fetchSequences("test.txt")
  #sequences1 = ["GTAGACC", "AGCGTAGA"]
  #enlève le \n
  seq2 = sequences1[1].split("\n")
  matrice = sequenceMatrix(len(sequences1[0]), len(seq2[0]))
  sequences1[1] = seq2[0]
  matrice = fillMatrix(matrice, sequences1)
  start = startingPos(matrice)[1]
  suffixPrefix = startingPos(matrice)[0]
  score = np.amax(matrice)
  print (matrice)
  sequences2 = fetchSequences("reads.fq")
  path,end = sequencePath(matrice,start, sequences1[0], sequences1[1])
  aligned, cheval = alignSequences(end, path, sequences1, start, matrice.shape)
  print("Sequence 1: " + aligned[0])
  print("Sequence 2: " + aligned[1])
  print("Chevauchement: " + str(cheval))
  print("Score: " + str(score))
  return None


if __name__ == "__main__":
  main()