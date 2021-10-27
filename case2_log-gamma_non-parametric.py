#!/usr/bin/env python



import pandas as pd
import numpy as np
import scipy 
import random
from sklearn.preprocessing import StandardScaler
import scipy.stats
import matplotlib.pyplot as plt
import shutil, os



# Load data and select first column

from sklearn import datasets
import argparse


parser = argparse.ArgumentParser(description='coverage size input')

	#-> required arguments
parser.add_argument('-Y', action='store', dest='COVERAGE_STORAGE', required=True, help='the input genome file in fasta format')

arg = parser.parse_args()



x=[12,3,6,13,5,1,12,11,10,6,3,5,25,10,5,12,8,3,17,16,24,3,16,10,9,18,13,17,14,13,17,11,22,4,21,19,7,24,8,9,16,21,10,5,12,5,23,15,7,5,15,3,7,12,25,12,14,10,3,12,13,7,18,4,1,12,15,2,11,6,18,14,12,21,25,8,30,8,18,11,10,19,20,12,27,6,14,6,3,3,22,1,14,18,13,8,25,9,2,31,11,29,18,16,3,4,24,24,9,6,19,19,5,2,12,18,23,15,5,11,10,20,19,16,8,18,19,20,17,17,22,1,33,10,23,15,4,15,17,17,37,17,15,1,16,19,21,2,6,10,13,14,18,17,10,9,26,2,12,10,19,15,17,26,18,20,3,20,2,14,12,4,14,13,16,22,8,22,4,3,17,21,18,2,10,16,21,9,2,8,20,12,6,15,7,27,10,15,3,3,22,17,9,16,21,13,5,19,40,26,10,13,20,4,1,20,13,16,4,3,10,13,2,9,12,15,17,14,5,16,6,3,20,9,10,9,17,19,17,2,29,20,17,19,16,16,12,5,7,2,5,31,13,18,14,10,14,25,17,29,15,2,18,13,3,5,23,9,15,34,22,13,20,1,15,14,1,16,15,14,17,9,3,23,15,15,40,17,18,6,10,6,16,4,22,13,12,9,15,2,12,17,12,13,6,26,14,18,20,11,20,9,8,16,19,13,16,12,16,27,6,4,18,12,19,16,16,13,10,15,12,12,19,1,12,17,13,15,21,6,12,6,20,31,24,17,7,19,19,17,14,4,17,1,16,20,32,12,8,8,33,20,4,21,11,15,29,34,14,9,10,2,8,18,4,4,17,4,12,5,18,16,14,12,3,20,12,19,30,9,3,19,11,7,4,17,4,17,34,13,14,11,16,2,20,16,22,9,14,8,3,13,27,14,18,19,20,12,15,21,10,8,29,17,14,22,18,2,2,12,18,2,18,11,15,31,24,19,22,23,22,9,15,28,5,9,20,31,14,16,21,21,3,22,23,24,3,13,17,25,11,11,5,5,36,13,16,3,13,13,11,33,3,24,10,9,13,21,16,20,6,4,17,21,9,13,9,20,15,6,1,27,16,3,6,15,6,8,15,17,17,20,11,10,17,19,22,6,5,9,11,24,26,19,14,5,14,12,14,23,17,20,13,21,20,14,12,9,7,11,16,16,16,7,1,18,13,2,19,25,10,5,1,12,14,15,14,22,25,17,20,4,8,14,16,20,15,19,12,16,15,19,21,17,8,13,25,8,20,4,16,18,9,22,9,20,3,7,2,19,1,25,8,15,7,3,6,13,3,24,14,7,17,13,18,19,19,11,11,23,15,20,10,22,6,2,18,15,22,1,12,15,3,20,12,11,8,19,11,24,15,12,18,26,28,17,2,27,14,27,18,12,18,9,21,17,10,5,5,11,11,8,15,1,12,14,17,6,21,3,10,8,8,9,10,36,8,21,14,34,14,11,6,13,18,13,6,19,20,12,20,19,23,22,6,25,6,15,3,10,4,14,12,21,14,15,10,16,11,16,6,18,7,14,16,18,21,14,20,15,1,5,17,13,21,9,34,12,22,18,19,21,19,12,4,21,6,16,23,8,11,11,7,13,0,0,0,0,0,0,0,0,0]














# Finding out from data
y_df = pd.DataFrame(x, columns=['Data'])
# y_df.describe()

cov= y_df.mean()/y_df.std()

# coverage determination
input_mean1=arg.COVERAGE_STORAGE

#print('input_mean',int(input_mean))



#print('g',g)

out_std1=int(input_mean1)/cov

#construct the simulated sampling distribution

# Statistical test
import random

sample_props1=[]
sample_data1=[]
#sample_data_props=[]
for _ in range(100000):
    sample=np.random.choice(x,size=677)
    if (sample.mean()>=13.354356 and sample.mean()<=13.930597):
        #if(sample.std()>=45.39 and sample.std()<=45.40):
        if(sample.mean()/sample.std()>=1.813 and sample.mean()/sample.std()<=1.815):
            sample_data1.append(sample)
            #sample_data_props.append(sample.mean())
    sample_props1.append(sample.mean())

k1=random.choice(sample_data1)
out1=k1.astype(int)

print('input_mean1',int(input_mean1))

sample_props=[]
sample_data=[]
#sample_data_props=[]
pb_runA = out1 / np.sum(out1, dtype=float)  # portion of copies of each sequence
print('pb_runA',len(pb_runA))

for _ in range(1000):
    _q = np.random.choice(677, size=int(input_mean1)*677, replace=True, p=pb_runA.astype(np.float))
    sim = np.zeros(len(pb_runA))
    for idx in _q:
        sim[idx] += 1
    if (int(input_mean1)>14):
        if (sim.mean()/sim.std()>=1.750 and sim.mean()/sim.std()<=1.850):
                sample_data.append(sim)
        sample_props.append(sim.mean())
    elif (int(input_mean1)<14):
        if (sim.mean()/sim.std()>=1.10 and sim.mean()/sim.std()<=1.80):
                sample_data.append(sim)
        sample_props.append(sim.mean())
    
if (int(input_mean1)!=14): 
   k=random.choice(sample_data)
   out=k.astype(int)

if (int(input_mean1)==14):
   k=random.choice(sample_data1)
   out=k.astype(int)

#p=out.tolist()
#y  = [val for sublist in p for val in sublist] 


y=out

file=open("log-gamma.txt","w+")
content=str(y)
file.write(content)
file.close()


print('y',y)


def simulate_indelsubs(read, sub_prob = 0.0, del_prob = 0.0, ins_prob = 0.0):
    '''
    add iid indels and substitions to read
    '''
    char_list = [c for c in read]
    pos_in_char_list = 0
    new_char_list = []
    alphabet = {}
    alphabet['all'] = ['A','C','G','T']
    alphabet['A'] = ['C','G','T']
    alphabet['C'] = ['A','G','T']
    alphabet['G'] = ['C','A','T']
    alphabet['T'] = ['C','G','A']
    while True:
        ins = (np.random.random_sample()<ins_prob)
        if ins:
            new_char_list.append(np.random.choice(alphabet['all']))
        else:
            if pos_in_char_list == len(char_list):# end of original read and not inserting
                break
            _del = (np.random.random_sample()<del_prob) 
            if _del:
                pos_in_char_list += 1
            else:
                sub = (np.random.random_sample()<sub_prob)
                if sub:
                    #q=np.random.choice(alphabet[char_list[pos_in_char_list]])
                    #new_char_list=new_char_list.get(q)
                    try:
                      new_char_list.append(np.random.choice(alphabet[char_list[pos_in_char_list]]))
                    except KeyError:
                      print "hi"
                else:
                    new_char_list.append(char_list[pos_in_char_list])
                pos_in_char_list += 1
    return ''.join(new_char_list)


i=0
b=677
for i in  range(677):
  r=y[i]
  #print('r:',r)

  j=1
  if (r>1):
     for j in range(1,r):
        shutil.copy("/home/subhasankar3593/DeepSimulator/sampled_readtry_DeepSimu/sampled_read_processed_genome_"+str(i)+".fasta", "/home/subhasankar3593/DeepSimulator/sampled_readtry_DeepSimu/sampled_read_processed_genome_"+str(b)+".fasta")
        f=open("/home/subhasankar3593/DeepSimulator/sampled_readtry_DeepSimu/sampled_read_processed_genome_"+str(b)+".fasta","rw+")
        
        lines =f.readlines()
        
        syn_seq = simulate_indelsubs(lines[1], sub_prob = 0.004, del_prob = 0.0085, ins_prob = 0.0005)
        #print syn_seq
        
        content=str(syn_seq)
       
        f.seek(0)
        f.truncate()
        p=str(b+1)
        L=[">Seq",p,"\n",content,"\n"]
       
        f.writelines(L)
        
        f.close()
        m=1
        b=b+m
  elif (r==0):
        os.remove("/home/subhasankar3593/DeepSimulator/sampled_readtry_DeepSimu/sampled_read_processed_genome_"+str(i)+".fasta")
   


