import editdistance
import random
import itertools
import multiprocessing as mp
import matplotlib.pyplot as plt

k = 2 #size of the alphabet
length = 10
alphabet = ['0','1','2','3']

outputs = []
for i in list(itertools.product(alphabet[:k], repeat=length)):
    outputs.append("".join(m for m in i))
random.shuffle(outputs)

ed = 6
pos = [0]*length

for s in outputs:
    neighbors = []
    for o in outputs:
        if editdistance.eval(o, s) == ed:
            neighbors.append(o)
    for n in neighbors:
        for letter in range(len(s)):
            if n[letter] != s[letter] :
                pos[letter] += 1
                
print (pos)
filename = 'ed ' + str(ed) 
plt.plot(pos)
plt.ylabel('error')
plt.title('ed ' + str(ed))
plt.savefig(filename)                