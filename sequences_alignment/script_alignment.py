import numpy as np
import sys
import os

"""
@CAOM
Program for aligning sequences:
Find a sequence in another longer one and return the matching sequence followed by the characters '-' in the no matching spaces:

input.txt:                      --->    output.txt:
RSIFVMIQIRDHIVM   FVMIQIRDH             ---FVMIQIRDH---
ELMKFVYLVQTENRL   YLVQTENRL             ------YLVQTENRL
RDTSMLHWFNRRSVA   MLHWFNRRS             ----MLHWFNRRS--
QWHFTEMIRHHGGKW   FTEMIRHHG             ---FTEMIRHHG---

Use the script in python 3 executing:
>>> python script_alignment.py input.txt

Script returns an output.txt file that contains the alignment
"""

#1 find() function returns the index when shorter sequence were matched
#2 operation in the range() gives the value of the ending no mathing places

arr_r = np.loadtxt(sys.argv[1], dtype=str)
f=open("output.txt","w+")
j=0
tr_r=[]
for i in arr_r:
  pre=pos=''
  for k in range(i[0].find(i[1])): #1
    pre+='-'
  for k in range(len(i[0])-(len(i[1])+i[0].find(i[1]))): #2
    pos+='-'
  tr_r.append(pre+i[1]+pos)
  f.write(" %s \n" % tr_r[j])
  j+=1

f.close()