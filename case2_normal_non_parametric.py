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

#x=[5,10,12,9,8,10,11,10,9,2,4,13,15,2,11,1,11,9,5,3,3,8,1,1,1,11,2,9,12,10,9,10,2,8,2,6,9,3,4,1,4,2,1,4,10,3,1,5,10,4,7,16,8,4,5,2,1,8,10,2,3,2,10,10,6,17,7,11,8,10,6,15,1,8,6,10,7,7,9,8,4,12,3,5,9,5,6,8,5,13,4,10,14,2,8,3,7,7,5,9,15,4,3,5,8,8,4,3,6,2,8,6,6,9,2,19,5,5,4,14,7,6,3,10,3,5,5,1,4,8,8,7,6,15,5,4,7,5,13,9,11,4,5,7,14,11,9,7,3,12,10,4,7,3,9,3,8,1,7,8,5,11,9,3,1,9,9,9,7,10,2,5,3,16,10,20,4,6,4,8,3,12,4,6,8,5,2,6,5,8,3,7,6,6,3,5,8,5,7,10,8,5,1,7,13,2,3,5,3,6,6,5,13,7,3,11,7,8,9,6,9,12,7,6,10,4,6,6,4,2,14,7,6,10,7,10,10,12,5,8,9,9,12,6,16,3,12,6,14,4,8,4,9,4,8,8,5,10,2,10,1,4,5,3,6,10,6,14,9,6,12,4,7,5,9,6,4,8,4,14,2,8,8,5,5,8,10,5,3,6,6,4,8,7,11,5,2,8,8,7,7,9,13,5,1,2,7,17,13,5,15,4,4,4,1,6,10,4,8,4,3,12,15,12,8,11,4,3,5,15,11,3,1,11,5,8,13,11,7,4,17,9,1,7,1,5,11,2,9,7,7,8,2,10,13,6,6,9,1,8,7,7,2,5,6,6,2,9,7,5,9,7,5,6,3,4,15,1,9,14,7,10,10,10,7,10,12,5,12,4,9,10,7,9,1,12,6,5,4,5,12,9,17,11,13,12,8,4,10,8,10,8,9,3,5,2,2,12,8,9,8,8,4,3,6,9,8,7,4,10,9,11,1,3,3,10,1,7,7,15,7,8,4,5,3,8,3,12,7,8,4,4,9,7,4,7,8,6,3,5,3,5,1,9,9,1,8,4,7,9,5,4,10,11,12,7,3,7,10,4,5,4,11,3,5,12,3,7,1,4,10,8,6,5,10,7,10,8,8,6,9,8,5,1,7,11,8,1,11,2,12,2,6,9,7,18,1,1,4,2,11,8,3,3,4,10,14,7,1,9,8,4,3,9,17,6,7,9,8,1,10,3,8,4,4,11,11,9,11,6,2,9,8,11,6,10,6,6,10,3,8,11,5,7,1,5,10,6,2,10,8,4,9,12,4,10,2,2,8,7,2,9,4,4,3,5,7,10,9,7,8,11,8,7,13,12,13,5,4,7,13,5,10,2,11,14,12,13,8,8,6,11,8,3,2,4,5,6,10,2,6,6,14,7,6,4,5,6,1,1,11,6,9,3,3,4,5,5,6,6,1,6,5,12,11,5,8,4,12,18,8,5,5,13,2,2,6,9,9,10,7,6,6,1,10,13,10,11,9,6,22,7,9,5,2,2,10,11,4,4,3,12,3,3,6,6,2,9,13,2,7,6,8,1,7,8,12,8,7,11,16,2,14,10,9,15,4,6,5,5,11,3,11,12,8,4,5,7,3,2,13,3,5,5,15,2,2,8,13,11,8,8,4,8,6,6,14,4,11,4,5,2,7,3,8,9,9,1,5,1,1,8,7,9,3,3,23,6,9,8,9,6,7,1,8,7,4,5,4,24,7,6,13,7,5,2,5,4,7,18,7,7,5,7,10,2,3,7,2,3,7,10,5,3,14,5,7,7,7,6,2,2,7,7,1,6,3,4,3,3,16,13,3,4,6,10,6,7,4,3,6,8,13,2,4,10,8,7,11,3,7,18,11,8,9,13,4,6,8,11,9,2,5,8,9,8,7,11,9,8,6,8,1,9,10,7,7,8,11,16,7,9,9,6,10,16,5,4,4,18,5,2,10,4,6,6,7,7,4,6,6,6,7,1,11,10,5,6,9,7,4,9,16,6,4,10,5,4,8,6,11,2,9,10,6,5,3,9,13,7,13,7,10,7,11,6,6,4,6,14,9,4,9,8,11,9,8,3,2,8,14,11,8,4,10,12,6,5,5,7,11,9,2,13,16,7,4,4,10,4,3,10,2,7,12,10,10,4,7,2,3,11,10,4,6,6,7,3,3,11,1,5,11,9,5,4,8,5,4,4,5,5,2,2,7,14,5,6,9,6,4,23,8,3,8,2,4,8,3,5,8,6,1,6,4,9,6,10,2,2,9,1,9,5,1,8,8,6,10,13,11,2,1,5,7,7,11,5,3,6,4,5,1,3,14,15,6,4,3,5,1,5,2,10,4,6,8,3,10,2,5,6,18,13,7,5,8,2,5,6,7,5,10,3,10,11,5,17,2,7,10,7,10,11,12,8,8,11,7,16,6,4,5,4,1,5,7,4,15,11,9,8,8,3,5,2,5,8,11,7,4,5,4,10,13,9,11,4,9,5,8,13,6,3,8,2,7,17,1,1,6,10,4,6,6,4,8,7,7,4,13,8,6,5,8,5,12,7,6,15,6,12,7,5,7,5,9,11,7,5,10,4,2,6,7,10,8,11,13,8,6,1,9,1,4,3,5,11,8,5,5,7,10,9,7,1,4,7,7,7,5,10,11,5,5,7,6,5,11,14,4,5,11,21,12,6,11,2,10,5,4,16,17,7,8,15,6,4,13,3,4,5,8,5,15,2,6,7,5,7,7,13,12,3,4,9,9,1,3,4,4,2,10,8,4,4,4,2,14,4,2,6,9,11,6,1,13,4,4,10,10,2,8,9,8,7,5,5,8,3,9,2,5,7,6,12,8,8,7,3,16,6,4,8,4,10,10,18,2,6,9,3,9,8,10,6,13,6,8,6,10,3,7,7,10,8,12,6,9,7,4,7,10,4,11,9,7,9,4,8,8,2,13,10,8,7,8,2,14,5,7,4,6,4,7,6,15,3,11,4,9,9,8,3,4,8,7,6,

x=[7,8,6,3,13,3,17,4,7,4,9,11,11,5,1,3,12,5,6,11,14,8,2,10,2,1,5,2,4,10,10,8,7,11,8,6,7,7,10,1,9,9,9,2,8,3,10,7,1,11,2,7,5,7,4,2,13,2,5,8,11,13,7,2,7,4,14,8,9,8,2,3,3,9,3,2,13,1,2,10,12,8,2,5,4,1,19,5,3,6,1,3,12,8,9,6,3,21,12,3,8,9,1,11,7,6,3,5,4,4,7,11,2,9,8,4,5,8,1,7,7,4,4,9,8,3,4,8,1,4,8,8,8,6,14,6,7,12,6,8,9,6,10,8,1,5,4,11,15,1,3,2,2,5,7,4,17,3,5,1,31,7,8,8,12,5,12,11,7,2,15,4,2,13,6,10,6,6,12,9,7,17,6,5,8,4,10,1,15,5,6,9,7,2,12,1,12,2,9,4,3,6,8,5,6,5,10,10,6,7,5,4,8,9,6,12,5,3,4,6,8,5,12,14,4,14,14,11,1,11,9,13,7,5,6,5,11,11,4,2,11,13,7,2,6,6,18,11,7,7,8,6,8,1,11,2,4,17,8,7,3,4,14,7,10,16,8,4,6,12,5,7,5,9,7,7,6,3,4,6,7,5,12,6,9,8,1,11,8,14,2,8,3,1,11,14,8,9,8,5,4,4,12,10,6,7,9,11,3,8,6,6,2,5,9,7,10,13,7,10,7,10,3,9,7,6,4,12,4,12,2,18,6,9,9,11,1,26,6,1,3,2,5,9,6,4,4,12,10,11,10,8,2,6,4,9,4,5,8,12,8,9,6,3,6,7,3,14,10,12,5,10,2,10,9,5,11,7,10,6,6,2,9,12,7,2,3,8,2,11,6,1,3,11,10,6,5,7,10,4,7,4,15,2,8,12,8,5,14,8,9,11,4,3,10,7,8,7,10,4,2,6,4,5,5,20,7,3,2,15,6,8,4,5,8,9,13,8,7,7,2,6,4,5,1,11,12,4,9,7,4,5,2,3,7,5,10,5,10,14,3,7,3,6,1,10,2,2,9,9,2,10,5,4,10,3,7,15,5,6,6,3,8,8,5,6,2,3,13,8,8,5,1,9,4,15,10,3,2,9,8,6,15,13,8,2,8,14,13,5,17,6,2,11,8,3,5,6,2,8,3,4,7,5,14,10,12,12,4,3,3,2,9,13,7,16,5,23,11,12,4,4,3,11,6,4,6,6,4,3,8,9,7,15,5,4,5,2,5,9,9,9,4,9,9,12,11,11,2,7,12,6,2,18,6,3,10,8,6,13,4,6,12,13,8,3,10,8,7,4,7,3,10,10,12,8,12,2,7,4,1,9,17,7,9,4,8,3,1,4,5,16,3,5,4,14,5,16,8,6,10,9,6,18,12,7,5,10,1,4,7,1,19,6,8,14,5,9,11,9,6,10,12,3,8,8,4,12,10,9,8,2,10,13,7,10,12,4,5,5,2,6,1,3,1,6,11,10,6,11,9,6,7,6,3,3,4,4,9,14,5,7,11,6,13,7,5,10,1,14,8,2,9,7,1,2,5,3,10,4,2,1,15,13,4,7,1,1,8,6,4,6,5,5,6,3,3,2,10,9,12,13,4,10,13,2,7,7,8,8,9,14,10,9,7,6,8,9,4,3,4,10,8,6,4,8,6,1,7,4,1,5,5,5,3,8,13,6,3,4,1,11,5,10,2,26,3,9,1,3,8,13,17,7,9,4,5,4,1,1,4,8,10,4,9,13,9,7,7,2,11,9,3,2,4,6,3,6,4,10,1,9,6,5,3,5,6,7,5,6,5,8,5,6,4,3,6,4,3,4,15,6,3,6,2,5,8,9,3,5,5,2,8,8,9,14,9,6,8,11,6,3,13,2,8,6,8,12,7,4,9,11,6,8,3,8,3,2,6,5,7,3,10,3,8,2,4,2,12,1,9,14,12,10,8,2,4,5,8,8,6,9,7,11,5,14,7,1,14,8,1,6,4,12,5,8,2,8,3,9,2,5,10,5,2,9,20,5,1,7,21,10,8,6,7,2,7,8,8,27,7,9,8,1,10,8,4,1,5,6,3,10,5,9,3,8,11,1,2,7,3,6,17,7,4,9,5,12,14,5,5,4,11,6,9,10,10,12,8,1,12,7,8,8,2,4,1,8,4,16,11,2,3,7,1,6,9,8,8,7,10,2,5,8,11,9,7,8,6,3,4,6,3,1,8,6,5,9,5,6,6,8,4,14,6,6,9,10,10,7,6,3,9,8,5,6,14,15,7,7,8,5,6,4,6,1,10,10,6,1,1,5,2,11,10,4,5,9,11,18,6,3,15,1,8,18,10,10,2,15,6,8,9,4,4,7,14,1,5,5,9,14,8,13,8,20,7,8,11,11,18,4,12,4,8,13,16,13,5,6,7,6,7,6,2,11,3,6,6,4,8,14,8,12,1,4,9,10,3,5,10,12,5,6,2,12,8,9,7,3,4,9,12,3,15,6,6,6,10,4,5,3,1,2,9,2,7,14,7,8,7,7,9,7,4,12,11,9,4,7,3,4,9,4,4,11,5,11,2,10,9,10,7,22,9,4,10,2,4,9,6,3,8,11,4,6,9,9,5,6,6,2,11,9,3,5,9,2,6,7,4,9,11,12,3,2,8,5,6,8,6,7,3,9,5,20,5,12,5,8,6,3,9,3,7,15,3,7,9,8,10,10,3,8,5,3,2,18,16,11,4,2,6,8,8,2,10,4,3,7,5,4,15,7,8,1,12,6,5,6,12,10,10,4,10,4,3,9,6,2,1,6,2,5,5,6,8,7,7,9,3,5,10,13,14,14,4,6,3,7,12,4,4,2,5,10,4,7,7,5,8,6,7,3,13,6,5,9,8,6,7,2,5,1,4,13,5,12,3,1,13,10,6,4,9,10,1,5,9,4,4,5,12,9,2,7,5,9,10,12,11,15,4,5,4,7,4,4,8,5,11,8,1,11,9,6,1,2,5,7,9,6,13,8,6,9,6,4,5,4,5,5,3,7,5,6,3,13,11,4,4,3,14,11,8,3,2,15,8,6,14,4,3,4,2,2,7,11,10,1,15,8,1,11,7,1,9,3,6,14,20,2,3,2,6,8,7,2,2,7,5,1,8,3,23,2,7,7,7,5,15,11,2,6,4,11,23,1,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]




# Finding out from data
y_df = pd.DataFrame(x, columns=['Data'])
# y_df.describe()

cov= y_df.mean()/y_df.std()

# coverage determination
input_mean1=arg.COVERAGE_STORAGE

print('input_mean1',int(input_mean1))



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
    if (sample.mean()>=6.661704 and sample.mean()<=6.982027):
        #if(sample.std()>=45.39 and sample.std()<=45.40):
        if(sample.mean()/sample.std()>=1.61 and sample.mean()/sample.std()<=1.63):
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
    if (int(input_mean1)>7):
        if (sim.mean()/sim.std()>=1.61 and sim.mean()/sim.std()<=1.63):
                sample_data.append(sim)
        sample_props.append(sim.mean())
    elif (int(input_mean1)<7):
        if (sim.mean()/sim.std()>=1.10 and sim.mean()/sim.std()<=1.77):
                sample_data.append(sim)
        sample_props.append(sim.mean())
  
if (int(input_mean1)!=7): 
   k=random.choice(sample_data)
   out=k.astype(int)

if (int(input_mean1)==7):
   k=random.choice(sample_data1)
   out=k.astype(int)


#p=out.tolist()
#y  = [val for sublist in p for val in sublist] 


y=out
#np.savetxt("normal_coverage5.txt",y)

file=open("normal_cov7_f1200_e0.txt","w+")
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
   


