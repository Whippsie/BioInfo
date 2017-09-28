
def fetchFullSeq(filename):
  fullsequence=""
  with open(filename, 'r') as f:
        line = f.readlines()
        fullsequence = line[1]
  return fullsequence

"""Return 0 if no stop is found, position of the stop if found"""
def splitSeq(fullsequence,start):
    result = [fullsequence[i:i+3] for i in range(start, len(fullsequence), 3)]
    stop = ["TAA","TAG","TGA","*"]
    count=0
    print (result)
    for ele in result:
        #is - identity testing, == - equality testing
        for tag in stop:
            if ele == tag :
                print ("in")
                return count
        count+=1
    return 0

### main
def main():
  splitSequence = []
  """The value equal to 1 is the right division"""
  splitSequence.append(splitSeq(fetchFullSeq("sequence.fasta"),0))
  splitSequence.append(splitSeq(fetchFullSeq("sequence.fasta"),1))
  splitSequence.append(splitSeq(fetchFullSeq("sequence.fasta"),2))

  for i in range (len(splitSequence)):
      print (splitSequence[i])

  return None

if __name__ == "__main__":
  main()