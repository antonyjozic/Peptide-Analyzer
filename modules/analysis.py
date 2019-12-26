import re

AA = {'A':0,'R':1,'N':2,'D':3,'C':4,'E':5,'Q':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
Print_AA = {0:'A', 1:'R', 2:'N', 3:'D', 4:'C', 5:'E', 6:'Q', 7:'G', 8:'H', 9:'I', 10:'L', 11:'K', 12:'M', 13:'F', 14:'P', 15:'S', 16:'T', 17:'W', 18:'Y', 19:'V'}

# Read in sequences, return list of sequences and total # of sequences
def get_seq(filename):
    with open(filename, 'r') as f:
        seq = 0
        counter = 0
        for line in f:
            counter += 1
            seq += line.strip()
        f.close()
        return final_seq, counter

# initialize the final list, a matrix to present positional occurances of AA residues in a given sequence
def init_final_list(seq_length):
    final_list = []
    for x in range(20):
        l1 = []
        for y in range(seq_length):
            l2 = int()
            l1.append(l2)
        final_list.append(l1)
    return final_list

# count occurances of AA residues in each sequence position, add counts to the final list
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

# returns sequence comprised of highest occuring amino acids at each position
def new_seq(final_list, seqlength):
    newSeq = ""
    for x in range(seqlength):
        temp = 0
        for y in range(19):
            if final_list[y + 1][x] > final_list[temp][x]:
                temp = y + 1
        newSeq += Print_AA[temp]
    return newSeq

# prints the matrices to an output.txt file, if no file is present, it will create one
def print_list(a_list, seqlength, mode):
    with open("output.txt", mode) as f:
        if mode == 'a':
            newseq = new_seq(a_list, seqlength)
            f.write('SEQUENCE COMPRISED OF HIGHEST OCCURING AMINO ACIDS IN EACH POSITION: ' + newseq + '\n'
            f.write('\n\n')
            f.write('COUNTS NORMALIZED TO TOTAL # OF SEQUENCES: \n')
        elif mode == 'w':
            f.write('\n\n')
            f.write('COUNTS OF POSITIONAL OCCURENCES OF AMINO ACIDS: \n')
        f.write("POS: \t")
        for i in range(seqlength):
            f.write(str(i + 1) + '\t')
        f.write('\n')
        f.write("-------------------------------------------------------------------------------------\n")

        for x in range(len(a_list)):
            if x in Print_AA:
                f.write(str(Print_AA[x]) + '\t')
            for y in range(len(a_list[x])):
                f.write(str(a_list[x][y]) + '\t')
            f.write('\n')
        return None


# normalization function that divides each count in final list by the total number of sequences
def normalize_list(final_list, seqlength, numseq):
    divided_list = []
    for x in range(20):
            l1 = []
            for y in range(seqlength):
                l2 = float()
                l1.append(l2)
            divided_list.append(l1)

    for x in range(20):
        for y in range(seqlength):
            divided_list[x][y] = final_list[x][y] / numseq
    return divided_list
