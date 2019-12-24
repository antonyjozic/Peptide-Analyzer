import re

AA = {'A':0,'R':1,'N':2,'D':3,'C':4,'E':5,'Q':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
Print_AA = {0:'A', 1:'R', 2:'N', 3:'D', 4:'C', 5:'E', 6:'Q', 7:'G', 8:'H', 9:'I', 10:'L', 11:'K', 12:'M', 13:'F', 14:'P', 15:'S', 16:'T', 17:'W', 18:'Y', 19:'V'}

#Read in sequences, return as a list
def get_seq(filename):
    #need to take in a tab-deliminated file in the format TAG *tab* SEQUENCE, single tag&sequence on one line, so .txt or .xlsx should work
    seq = []
    counter=0
    f = open(filename, 'r')
    counter = 0
    for line in f:
        counter+=1
        seq += line.strip().split('\t')
        if seq[len(seq)-1] == 'a':
            counter-=1
    f.close()
   
    final_seq = seq[1::2]
    
    return final_seq,counter

#initialize the final list, a matrix to present positional occurances of AA residues in a given sequence
def init_final_list(seq_length):
    final_list = []
    for x in range(20):
        l1 = []
        for y in range(seq_length):
            l2 = int()
            l1.append(l2)
        final_list.append(l1)
    return final_list

#count occurances of AA residues in each sequence position, add counts to the final list
def count_residues(seq, final_list):
    counter = 0
    for x in seq:
        c2 = 0
        for y in seq[counter]:
            if y in AA:
                final_list[AA[y]][c2] += 1
            c2 += 1
        counter += 1
    return final_list

#print the final list to the command line window
#
#   SAME THING AS THE NEXT FUNCTION
#
#
def print_list(final_list):
    print()
    print("POS:", end = ' ')
    for i in range(7):
        print('{:<3}   ' .format(i + 1), end = ' ')
    print()
    print("-------------------------------------------------")

    for x in range(len(final_list)):
        if x in Print_AA:
            print('{:s}:  ' .format(Print_AA[x]), end = ' ')
        for y in range(len(final_list[x])):
            print('{:<3}   ' .format(final_list[x][y]), end = ' ')
        print()
    return None

#print the divided list (probabilities of occurance) to the command line window
#
#   NEED TO CHANGE THIS TO OUTPUT TO A TXT FILE, TAB-DELIMINATED TO MAKE GRAPHING EASIER
#   WE CAN ALSO TRY TO IMPLEMENT AN R SCRIPT TO CREATE A HEATMAP FROM IT
#
def print_divided_list(divided_list):
    print()
    print("POS:", end = ' ')
    for i in range(7):
        print('{:<8}   ' .format(i + 1), end = ' ')
    print()
    print("-------------------------------------------------------------------------------------")

    for x in range(len(divided_list)):
        if x in Print_AA:
            print('{:s}:  ' .format(Print_AA[x]), end = ' ')
        for y in range(len(divided_list[x])):
            print('{:<8.3f}   ' .format(divided_list[x][y]), end = ' ')
        print()
    return None

#add normalization function that divides each count in final list by 20
def normalize_list(final_list, seqlength, numseq):
    divided_list = []
    for x in range(20):
            l1 = []
            for y in range(seqlength):
                l2 = float()
                l1.append(l2)
            divided_list.append(l1)

    for x in range(20):
        for y in range(7):
            divided_list[x][y] = final_list[x][y] / numseq
    return divided_list

#add function that creates a sequence based on most occuring AA in each position
#
#   PRINT TO OUTPUT FILE
#

def new_seq(final_list):
    newSeq = ""
    for x in range(7):
        temp = 0
        for y in range(19):
            if final_list[y + 1][x] > final_list[temp][x]:
                temp = y + 1
        newSeq += Print_AA[temp]
    print()
    print("New Sequence: ", newSeq)
    print()
    return None
