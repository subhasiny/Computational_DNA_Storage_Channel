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








x=[7,12,22,12,8,9,16,10,6,22,16,8,23,8,10,14,5,23,17,18,2,17,6,6,13,5,1,13,8,11,13,17,3,10,19,4,4,14,6,26,11,11,17,2,10,37,16,16,4,20,21,9,15,6,22,4,8,2,8,11,2,21,6,8,11,13,8,20,7,2,3,4,8,10,3,3,6,6,17,14,23,9,20,2,3,7,4,12,17,11,16,13,9,16,3,21,16,6,7,15,20,1,10,17,7,6,24,29,14,13,15,3,16,2,1,2,22,21,17,9,2,26,4,12,1,5,3,7,12,40,29,4,17,9,9,2,17,14,20,17,2,18,2,18,14,6,8,24,2,6,18,23,10,14,13,10,7,21,4,14,8,7,1,8,20,11,4,1,9,5,13,3,8,9,9,10,14,4,5,20,14,14,13,16,18,27,4,9,7,18,7,21,2,15,11,13,3,17,5,38,7,1,4,10,2,13,3,20,3,7,19,3,11,21,14,29,13,22,19,13,16,5,2,7,4,3,2,4,7,14,2,11,12,19,6,19,10,10,20,2,16,14,2,16,17,4,16,5,32,11,13,2,24,5,19,12,19,28,5,12,5,8,4,4,21,12,19,14,26,2,18,24,8,3,11,2,27,7,2,1,2,6,16,1,6,8,2,10,4,1,6,9,4,3,4,15,15,15,10,11,17,11,11,13,23,13,25,18,7,18,9,1,5,1,4,7,1,13,25,4,19,7,6,19,17,15,1,20,5,6,9,9,3,6,6,25,3,1,3,25,1,14,21,19,13,18,16,6,4,15,4,16,5,16,6,23,6,12,17,9,8,18,2,12,8,7,1,18,11,2,2,28,7,13,1,18,21,3,5,12,4,8,9,20,3,23,7,6,26,18,2,20,8,8,19,21,10,15,9,19,25,15,16,6,3,4,15,5,15,20,11,14,3,4,11,7,34,15,2,7,8,9,11,7,11,20,12,14,14,14,5,2,13,1,19,6,21,21,18,15,26,10,1,19,1,4,10,13,12,2,3,13,22,11,15,27,23,10,24,5,21,23,15,35,8,11,15,14,23,10,17,5,2,9,12,10,6,12,9,13,5,16,7,8,11,11,20,12,7,21,4,13,4,19,4,22,16,7,27,23,7,2,14,5,20,14,11,7,4,21,11,18,10,4,3,3,37,7,11,12,10,19,11,13,40,10,21,1,5,15,22,2,10,16,21,11,16,4,2,3,15,10,6,10,4,14,21,7,2,3,3,16,3,12,13,8,23,10,4,13,2,7,5,14,2,14,7,16,16,1,2,13,9,2,12,1,22,2,2,6,3,11,13,7,17,13,5,10,8,3,19,31,13,31,1,17,14,2,31,9,11,10,17,18,21,5,22,5,20,9,7,5,18,3,10,19,28,42,10,10,34,10,12,17,5,6,15,3,9,10,4,10,9,10,3,7,19,7,7,1,12,10,10,16,10,25,9,26,15,9,14,2,9,3,14,9,1,8,12,18,24,14,9,17,19,23,11,19,1,1,12,18,9,12,8,10,24,2,22,7,12,18,15,8,17,10,15,16,9,30,3,12,14,3,15,7,10,4,18,1,1,15,6,8,14,4,7,4,30,5,5,26,14,8,6,15,21,8,16,17,13,26,9,18,11,14,18,11,7,17,1,16,4,20,12,15,10,13,12,7,1,11,31,8,16,10,11,25,3,17,13,19,9,15,12,13,25,18,1,15,3,14,33,19,22,6,28,1,9,11,12,21,10,12,21,5,17,24,14,16,20,29,15,15,11,13,24,1,3,23,9,15,3,8,16,20,20,20,15,2,16,7,25,4,11,3,22,14,1,5,14,5,12,10,15,26,20,11,26,11,3,1,9,17,5,6,9,3,12,4,5,3,8,6,12,3,4,14,18,14,18,20,18,8,7,2,19,11,6,16,12,13,3,11,1,8,2,1,22,10,2,6,2,17,0,0,0,0,0,0,0,0,0,0,0,0,0]











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
    if (sample.mean()>=11.10987 and sample.mean()<=11.695418):
        #if(sample.std()>=45.39 and sample.std()<=45.40):
        if(sample.mean()/sample.std()>=1.49 and sample.mean()/sample.std()<=1.50):
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
    if (int(input_mean1)>11):
        if (sim.mean()/sim.std()>=1.49 and sim.mean()/sim.std()<=1.50):
                sample_data.append(sim)
        sample_props.append(sim.mean())
    elif (int(input_mean1)<11):
        if (sim.mean()/sim.std()>=1.10 and sim.mean()/sim.std()<=1.50):
                sample_data.append(sim)
        sample_props.append(sim.mean())
    
if (int(input_mean1)!=11): 
   k=random.choice(sample_data)
   out=k.astype(int)

if (int(input_mean1)==11):
   k=random.choice(sample_data1)
   out=k.astype(int)


#p=out.tolist()
#y  = [val for sublist in p for val in sublist] 


y=out

file=open("beta.txt","w+")
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
   


