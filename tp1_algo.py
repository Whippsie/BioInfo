
class Case:
    def __init__(self, val=0):
        self.valeur = val
        self.listeFleches = []


#seq1 = "AGTCAGGCA"
#seq2 = "GCTACCA"

seq1 = "AGT"
seq2 = "GCTA"

poidsMatch = 4
poidsMismatch = -4
poidsIndel = -8

"""seq1 = input ("Entrez la premiere sequence")"""
"""seq2 = input ("Entrez la deuxieme sequence")"""

lenSeq1 = len(seq1)
lenSeq2 = len(seq2)

matriceTotale = [[Case() for i in range(lenSeq1+1)] for j in range(lenSeq2 + 1)]


def main():
    initialiserCol()
    initialiserLig()
    calculerScore()
    print("\n")
    prettyPrint()


def calculerScore():
    prettyPrint()
    for i in range (1, lenSeq1+2):
        for j in range (1, lenSeq2):
            res1 = matriceTotale[i-1][j].valeur+poidsIndel
            res2 = matriceTotale[i][j-1].valeur+poidsIndel
            res3 = matriceTotale[i-1][j-1].valeur
            if seq1[j-1] == seq2[i-1]:
                res3 += poidsMatch
            else:
                res3 += poidsMismatch

            matriceTotale[i][j].valeur = max(res1,res2,res3)
            #TOFO: vérifier quand plus qu'un égal
            current = matriceTotale[i][j]

            #fleches
            if res1 == current.valeur:
               current.listeFleches.append((i-1, j))
            if res2 == current.valeur:
                current.listeFleches.append((i, j-1))
            if res3 == current.valeur:
                current.listeFleches.append((i-1, j-1))


def initialiserCol():
    for i in range(lenSeq1+1):
        matriceTotale[0][i] = Case()


def initialiserLig():
    for i in range (lenSeq2+1):
        matriceTotale[i][0] = Case(i*poidsIndel)


def prettyPrint():
    for i in range(0, lenSeq1+2):
        for j in range(0, lenSeq2):
            #print("Rangée :", i, " Colonne :", j, " Valeur :", matriceTotale[i][j].valeur," \n ")
            print(matriceTotale[i][j].valeur," ")

    #for val in matriceTotale:
    #matrix = "\n".join([" ".join([matriceTotale[i][j].valeur for i in range(lenSeq1)]) for y in range(lenSeq2)])
    #print (matrix)


if __name__ == "__main__":
    main()

