import matplotlib.pyplot as plt
import matplotlib.font_manager
import numpy as np
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment as MSA


AA = {'A':0,'R':1,'N':2,'D':3,'C':4,'E':5,'Q':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
PRINT_AA = {0:'A', 1:'R', 2:'N', 3:'D', 4:'C', 5:'E', 6:'Q', 7:'G', 8:'H', 9:'I', 10:'L', 11:'K', 12:'M', 13:'F', 14:'P', 15:'S', 16:'T', 17:'W', 18:'Y', 19:'V'}

def get_seq(filename: str) -> list:
    '''
    Takes in the path to a filename, and extracts each sequence from the file.
    File format is a .txt file with a single sequence per line.

    Sequences are returned in a list.
    '''
    with open(filename, 'r') as f: #open the file
        seq = []
        for line in f:
            seq.append(line.strip())
        f.close()
    return seq

def init_final_list(seq_length: int) -> list:
    '''
    Initializes a 2-dimensional matrix of null int values that will eventually hold the counts of each amino acid at each position.
    The matrix contains 20 rows (one for each amino acid) and X columns, where X is the length of the sequences.
    '''
    final_list = []
    for x in range(20):
        l1 = []
        for y in range(seq_length):
            l2 = int()
            l1.append(l2)
        final_list.append(l1)
    return final_list

def count_residues(seq: str, final_list: list) -> list:
    '''
    Counts the residues at each position in each sequence, and iterates the respective i,j value in the count matrix.
    '''
    counter = 0
    for x in seq:
        c2 = 0
        for y in seq[counter]:
            if y in AA:
                final_list[AA[y]][c2] += 1
            c2 += 1
        counter += 1
    return final_list

def consensus_seq(final_list: list, seqlength: int) -> str:
    '''
    Determines the consensus sequence of all the peptide sequences.
    This is defined as the sequence comprised of the most occuring amino acids at each position.
    '''
    newSeq = ""
    for x in range(seqlength):
        temp = 0
        for y in range(19):
            if final_list[y + 1][x] > final_list[temp][x]:
                temp = y + 1
        newSeq += PRINT_AA[temp]
    return newSeq

def print_list(a_list: list, seqlength: int, mode: str, out_path: str) -> None:
    '''
    Prints the final matrix to a log file, if no file is present, one will be created.
    '''
    path = out_path + "/peptide-analysis-output.txt"
    with open(path, mode) as f:
        if mode == 'a':
            newseq = consensus_seq(a_list, seqlength)
            f.write('SEQUENCE COMPRISED OF HIGHEST OCCURING AMINO ACIDS IN EACH POSITION: ' + newseq + '\n')
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
            if x in PRINT_AA:
                f.write(str(PRINT_AA[x]) + '\t')
            for y in range(len(a_list[x])):
                f.write(str(a_list[x][y]) + '\t')
            f.write('\n')
        return None

def normalize_list(final_list: list, seqlength: int, numseq: int) -> list:
    '''
    Normalizes the count matrix to the total number of sequences, this provides a proportion of occurence for each amino acid in each position.
    '''
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

def plot_heatmap(final_list: list, seqlength: int, outpath: str, fname: str) -> None:
    '''
    Plots the proportional occurence of each amino acid in each position.
    '''
    plot_list = []
    for i in range(len(final_list)):
        plot_list.append(np.array(i))
    
    plot_list = np.array(plot_list)
    plt.imshow(final_list, cmap=plt.cm.OrRd, interpolation='nearest')
    
    xticks = []
    for x in range(seqlength):
        xticks.append(str(x+1))

    plt.yticks([*PRINT_AA],[*AA])
    plt.xticks(range(seqlength),xticks)
    
    plt.title('Proportional Occurence of Amino Acids')
    plt.xlabel('Position in ' + str(seqlength) +'mer')
    plt.ylabel('Amino Acid')

    plt.colorbar()

    plt.show()
    return None

def check_input_file(fname: str) -> None:
    '''
    Checks if the input file is properly formatted. If it is not, an exception is raised.
    '''
    with open(fname) as f:
        for line in f:
            x = line.strip()
            for char in x:
                if char not in AA:
                    raise Exception("Input file is not properly formatted")
                    break
        return None

def comp(seq: str) -> float:
    '''
    Computes sequence complexity.

    TO-DO:
        Use this as a measure in conjunction with other measures to perform clustering?
    '''
    for x in seq:
        AA[x] +=1
    
    count = 0
    temp = 0.0
    running = 0.0
    for k in AA:
        if count == 0:
            temp = math.factorial(AA[k])
            count+=1
        else:
            running = math.factorial(AA[k])*temp
            temp = running

    comp = (1/7)*math.log(((math.factorial(7))/running),20)
    return comp

def average(li: list) -> float:
    '''
    Compute an average from a list of floats.
    '''
    run = 0.0
    for x in li:
        run += x

    return run/len(li)


#what more information can be gained?
#   - alignment can give information about sequence homology -- MSA
#   - MD can give information about how the peptide sequence might interact with a receptor
#       - can also use the approach data to fine-tune a specific sequence for further testing
#   - clustering the AA based on local alignments can help differentiate the into peptide characteristics (hydrophobicity, etc)
#       - find a way to plot this quantitatively, then use k-means clustering.

def main(fname,out_path):
    seqs = get_seq(fname)
    numseq = len(seqs)
    seqlength = len(seqs[0])
    final_list = init_final_list(seqlength)
    counts = count_residues(seqs,final_list)
    normalized = normalize_list(final_list,seqlength,numseq)
    print_list(counts,seqlength, 'w',out_path)
    print_list(normalized,seqlength, 'a',out_path)
    plot_heatmap(normalized, seqlength, out_path, fname)
    return None
