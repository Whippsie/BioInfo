
class Case:
    def __init__(self):
        self.valeur = 0
        self.listeFleches = []


seq1 = "AGTCAGGCA"
seq2 = "GCTACCA"

poidsMatch = 4
poidsMismatch = -4
poidsIndel = -8

"""seq1 = input ("Entrez la premiere sequence")"""
"""seq2 = input ("Entrez la deuxieme sequence")"""

lenSeq1 = len(seq1)
lenSeq2 = len(seq2)

"""matriceTotale = list([range(0, lenSeq1+1), range(0, lenSeq2+1)])"""
matriceTotale = [[Case() for j in range(lenSeq1 + 1)] for i in range(lenSeq2 + 1)]
from pprint import pprint
pprint(matriceTotale)
#matriceTotale = list([Case() for j in range(lenSeq1 + 1), for i in range(lenSeq2 + 1)])



def main():
    initialiserCol()
    initialiserLig()
    calculerScore()
    prettyPrint()


def calculerScore():
    """write me plz"""
    for i in range (1, lenSeq1+1):
        for j in range (1, lenSeq2+1):
            print (matriceTotale[i-1][j].valeur)
            res1 = matriceTotale[i-1][j].valeur+poidsIndel
            res2 = matriceTotale[i][j-1].valeur+poidsIndel
            res3 = matriceTotale[i-1][j-1].valeur

            if seq1[i] == seq1[j]:
                res3 += poidsMatch
            else:
                res3 += poidsMismatch

            matriceTotale[i][j] = max(res1,res2,res3)
            current = matriceTotale[i][j]

            #fleches
            if res1 == current.valeur:
               current.listeFleches.append((i-1, j))
            if res2 == current.valeur:
                current.listeFleches.append((i, j-1))
            if res3 == current.valeur:
                current.listeFleches.append((i-1, j-1))


def initialiserCol():
    for i in range(0,lenSeq1+1):
        matriceTotale[0][i] = 0


def initialiserLig():
    for i in range (0,lenSeq2+1):
        matriceTotale[i][0] = i*poidsIndel

def prettyPrint():
    for i in range(0, lenSeq1 + 1):
        for j in range(0, lenSeq2 + 1):
            print(matriceTotale[i][j].valeur," ")
        print ("\n")


if __name__ == "__main__":
    main()

