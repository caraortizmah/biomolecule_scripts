import numpy as np
import sys
import os

"""
@CAOM
Program for aligning sequences:
Find a sequence in another longer one and return the longer sequence organized such as the coincidence in the shorter one. Void spaces are replaced by the characters '-':

input.txt:                      --->    output.txt:
RSIFVMIQIRDHIVM   FVMIQIRDH             ---RSIFVMIQIRDHIVM---
ELMKFVYLVQTENRL   YLVQTENRL             ELMKFVYLVQTENRL------
RDTSMLHWFNRRSVA   MLHWFNRRS             --RDTSMLHWFNRRSVA----
QWHFTEMIRHHGGKW   FTEMIRHHG             ---QWHFTEMIRHHGGKW---

Use the script in python 3 executing:
>>> python align_sequence.py input.txt

Script returns an output.txt file that contains the alignment
"""

#1 difference of lenght between two columns of the input
#2 find() function returns the index when shorter sequence were matched
#3 operation in the range() gives the value of the ending no mathing places

arr_r = np.loadtxt(sys.argv[1], dtype=str)
f=open("output.txt","w+")
j=0
tr_r=[]
for i in arr_r:
  diff=len(i[0])-len(i[1]) #1
  pre=pos=''
  for k in range(diff-i[0].find(i[1])): #2
    pre+='-'
  for k in range(diff-(len(i[0])-(len(i[1])+i[0].find(i[1])))): #3
    pos+='-'
  tr_r.append(pre+i[0]+pos)
  f.write(" %s \n" % tr_r[j])
  j+=1

f.close()