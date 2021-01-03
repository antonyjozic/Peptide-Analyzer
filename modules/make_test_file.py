import random

random.seed()

Print_AA = {0:'A', 1:'R', 2:'N', 3:'D', 4:'C', 5:'E', 6:'Q', 7:'G', 8:'H', 9:'I', 10:'L', 11:'K', 12:'M', 13:'F', 14:'P', 15:'S', 16:'T', 17:'W', 18:'Y', 19:'V'}

seqs = []
for x in range(30):
    seq = ''
    for y in range(12):
        rint = random.randint(0,19)
        seq += (Print_AA[rint])
    seqs.append(seq)

with open('test.txt','w') as of:
    for seq in seqs:
        of.write(f"{seq}\n")
