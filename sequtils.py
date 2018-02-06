#Function to reverse a seq
#input: ACG+TGA
#output: AG+TGCA
def reverse(seq):
        #Make a string from the seq list
        seq = ''.join(seq)
        #Reverse the seq
        seqr = seq[::-1]
        #Find the positions of + characters but shifted one to the left
        positions = [i-1 for i, plus in enumerate(seqr) if plus == '+']
        #Replace the + chars with X
        seqr = seqr.replace("+",'X')
        seqrN = ''
        for i,c in enumerate(seqr):
                #Remove the X's
                if c == 'X':
                        c=''
                #Put in a plus at the positions
                if(i in positions):
                        c='+'+c
                seqrN = seqrN+c
        return toList(seqrN)


#Function to make the complementary sequence of the input sequence
#The complementary base of a LNA base is a normal base.
#Input: ACG+TATC
#Output:TGCATAG
def makeComp(seq):
        seqN = []
        for c in seq:
                if(c == 'A' or c == "+A"):
                        seqN.append('T')
                if(c == 'T' or c == "+T"):
                        seqN.append('A')
                if(c == 'G' or c == "+G"):
                        seqN.append('C')
                if(c == 'C' or c == "+C"):
                        seqN.append('G')
        return seqN

#Function to check if two sequences are complementary to each other
#Input: ACGT,TGCT
#Output: True
def checkComp(seq1,seq2):
        if(makeComp(seq1) == seq2 ):
                 return True
        else:
                return False

#Function to convert a sequence in string format to a list for ease of use
#Input: A+CGT
#Output: A,+C,G,T
def toList(seq):
        LNA = ''
        seqL = []
        for c in seq:
                if c == '+':
                        LNA = '+'
                else:
                        c=LNA+c
                        LNA = ''
                        seqL.append(c)
        return seqL

