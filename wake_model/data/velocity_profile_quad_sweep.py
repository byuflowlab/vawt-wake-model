import numpy as np
import matplotlib.pyplot as plt
# from velocity_profile_read_skew_def_quad import fit
from velocity_optimize_SMG import fit
import csv
from os import path

s = np.array([])
t = np.array([])

solidity = np.array(['s1'])#,'s2','s3','s4','s5'])
tsr = np.linspace(150,700,23)

tsr_write = tsr/100.
solidity_write = np.array([0.15])#,0.25,0.5,0.75,1.0])

for i in range(np.size(solidity)):
    for j in range(np.size(tsr)):
        s = np.append(s,solidity[i])
        t = np.append(t,str(int(tsr[j])))
        
s1length = np.array([210,210,205,196,185,178,170,165,160,145,140,123,115,112,108,101,101,90,85,80,78,75,70])
# s1length = np.array([197,193,185,176,140,146,126,114,103,96,100,86,77,72,70,68,60,64,54,50,47,45,44])
# s1length = np.array([185,150,100,95,83,76,72,63,60,49,50,41,39,36,34,33,30,31,28,30,29,28,27])
# s1length = np.array([145,100,73,60,53,44,42,37,38,30,33,26,22,24,23,21,21,19,24,23,22,21,20])
# s1length = np.array([78,70,52,43,37,32,29,27,26,23,20,20,23,21,20,19,19,18,18,16,16,15,14])

slength = np.array([])
slength = np.append(slength,s1length)
# slength = np.append(slength,s2length)
# slength = np.append(slength,s3length)
# slength = np.append(slength,s4length)
# slength = np.append(slength,s5length)
slength = slength*1.

men1 = np.array([])
men2 = np.array([])
men3 = np.array([])
spr1 = np.array([])
spr2 = np.array([])
spr3 = np.array([])
spr4 = np.array([])
rat1 = np.array([])
rat2 = np.array([])
tns1 = np.array([])
tns2 = np.array([])
scl1 = np.array([])
scl2 = np.array([])
scl3 = np.array([])

for i in range(np.size(s)):
    mend,sprd,scld,ratd,tnsd = fit(s[i],t[i],slength[i],False,'mac',4,False)
    
    men1 = np.append(men1,mend[0])
    men2 = np.append(men2,mend[1])
    men3 = np.append(men3,mend[2])
    spr1 = np.append(spr1,sprd[0])
    spr2 = np.append(spr2,sprd[1])
    spr3 = np.append(spr3,sprd[2])
    spr4 = np.append(spr4,sprd[3])
    rat1 = np.append(rat1,ratd[0])
    rat2 = np.append(rat2,ratd[1])
    tns1 = np.append(tns1,tnsd[0])
    tns2 = np.append(tns2,tnsd[1])
    scl1 = np.append(scl1,scld[0])
    scl2 = np.append(scl2,scld[1])
    scl3 = np.append(scl3,scld[2])
    print s[i],t[i]

s1_men1 = np.array([])
s1_men2 = np.array([])
s1_men3 = np.array([])
s1_spr1 = np.array([])
s1_spr2 = np.array([])
s1_spr3 = np.array([])
s1_spr4 = np.array([])
s1_rat1 = np.array([])
s1_rat2 = np.array([])
s1_tns1 = np.array([])
s1_tns2 = np.array([])
s1_scl1 = np.array([])
s1_scl2 = np.array([])
s1_scl3 = np.array([])

# s2_men1 = np.array([])
# s2_men2 = np.array([])
# s2_men3 = np.array([])
# s2_spr1 = np.array([])
# s2_spr2 = np.array([])
# s2_spr3 = np.array([])
# s2_rat1 = np.array([])
# s2_rat2 = np.array([])
# s2_tns1 = np.array([])
# s2_tns2 = np.array([])
# s2_scl1 = np.array([])
# s2_scl2 = np.array([])
# s2_scl3 = np.array([])
# 
# s3_men1 = np.array([])
# s3_men2 = np.array([])
# s3_men3 = np.array([])
# s3_spr1 = np.array([])
# s3_spr2 = np.array([])
# s3_spr3 = np.array([])
# s3_rat1 = np.array([])
# s3_rat2 = np.array([])
# s3_tns1 = np.array([])
# s3_tns2 = np.array([])
# s3_scl1 = np.array([])
# s3_scl2 = np.array([])
# s3_scl3 = np.array([])
# 
# s4_men1 = np.array([])
# s4_men2 = np.array([])
# s4_men3 = np.array([])
# s4_spr1 = np.array([])
# s4_spr2 = np.array([])
# s4_spr3 = np.array([])
# s4_rat1 = np.array([])
# s4_rat2 = np.array([])
# s4_tns1 = np.array([])
# s4_tns2 = np.array([])
# s4_scl1 = np.array([])
# s4_scl2 = np.array([])
# s4_scl3 = np.array([])
# 
# s5_men1 = np.array([])
# s5_men2 = np.array([])
# s5_men3 = np.array([])
# s5_spr1 = np.array([])
# s5_spr2 = np.array([])
# s5_spr3 = np.array([])
# s5_rat1 = np.array([])
# s5_rat2 = np.array([])
# s5_tns1 = np.array([])
# s5_tns2 = np.array([])
# s5_scl1 = np.array([])
# s5_scl2 = np.array([])
# s5_scl3 = np.array([])

k = 0
for i in range(np.size(tsr)):
    s1_men1 = np.append(s1_men1,men1[k])
    s1_men2 = np.append(s1_men2,men2[k])
    s1_men3 = np.append(s1_men3,men3[k])
    s1_spr1 = np.append(s1_spr1,spr1[k])
    s1_spr2 = np.append(s1_spr2,spr2[k])
    s1_spr3 = np.append(s1_spr3,spr3[k])
    s1_spr4 = np.append(s1_spr4,spr4[k])
    s1_rat1 = np.append(s1_rat1,rat1[k])
    s1_rat2 = np.append(s1_rat2,rat2[k])
    s1_tns1 = np.append(s1_tns1,tns1[k])
    s1_tns2 = np.append(s1_tns2,tns2[k])
    s1_scl1 = np.append(s1_scl1,scl1[k])
    s1_scl2 = np.append(s1_scl2,scl2[k])
    s1_scl3 = np.append(s1_scl3,scl3[k])
    k += 1
# for i in range(np.size(tsr)):
#     s2_men1 = np.append(s2_men1,men1[k])
#     s2_men2 = np.append(s2_men2,men2[k])
#     s2_men3 = np.append(s2_men3,men3[k])
#     s2_spr1 = np.append(s2_spr1,spr1[k])
#     s2_spr2 = np.append(s2_spr2,spr2[k])
#     s2_spr3 = np.append(s2_spr3,spr3[k])
#     s2_rat1 = np.append(s2_rat1,rat1[k])
#     s2_rat2 = np.append(s2_rat2,rat2[k])
#     s2_tns1 = np.append(s2_tns1,tns1[k])
#     s2_tns2 = np.append(s2_tns2,tns2[k])
#     s2_scl1 = np.append(s2_scl1,scl1[k])
#     s2_scl2 = np.append(s2_scl2,scl2[k])
#     s2_scl3 = np.append(s2_scl3,scl3[k])
#     k += 1
# for i in range(np.size(tsr)):
#     s3_men1 = np.append(s3_men1,men1[k])
#     s3_men2 = np.append(s3_men2,men2[k])
#     s3_men3 = np.append(s3_men3,men3[k])
#     s3_spr1 = np.append(s3_spr1,spr1[k])
#     s3_spr2 = np.append(s3_spr2,spr2[k])
#     s3_spr3 = np.append(s3_spr3,spr3[k])
#     s3_rat1 = np.append(s3_rat1,rat1[k])
#     s3_rat2 = np.append(s3_rat2,rat2[k])
#     s3_tns1 = np.append(s3_tns1,tns1[k])
#     s3_tns2 = np.append(s3_tns2,tns2[k])
#     s3_scl1 = np.append(s3_scl1,scl1[k])
#     s3_scl2 = np.append(s3_scl2,scl2[k])
#     s3_scl3 = np.append(s3_scl3,scl3[k])
#     k += 1
# for i in range(np.size(tsr)):
#     s4_men1 = np.append(s4_men1,men1[k])
#     s4_men2 = np.append(s4_men2,men2[k])
#     s4_men3 = np.append(s4_men3,men3[k])
#     s4_spr1 = np.append(s4_spr1,spr1[k])
#     s4_spr2 = np.append(s4_spr2,spr2[k])
#     s4_spr3 = np.append(s4_spr3,spr3[k])
#     s4_rat1 = np.append(s4_rat1,rat1[k])
#     s4_rat2 = np.append(s4_rat2,rat2[k])
#     s4_tns1 = np.append(s4_tns1,tns1[k])
#     s4_tns2 = np.append(s4_tns2,tns2[k])
#     s4_scl1 = np.append(s4_scl1,scl1[k])
#     s4_scl2 = np.append(s4_scl2,scl2[k])
#     s4_scl3 = np.append(s4_scl3,scl3[k])
#     k += 1
# for i in range(np.size(tsr)):
#     s5_men1 = np.append(s5_men1,men1[k])
#     s5_men2 = np.append(s5_men2,men2[k])
#     s5_men3 = np.append(s5_men3,men3[k])
#     s5_spr1 = np.append(s5_spr1,spr1[k])
#     s5_spr2 = np.append(s5_spr2,spr2[k])
#     s5_spr3 = np.append(s5_spr3,spr3[k])
#     s5_rat1 = np.append(s5_rat1,rat1[k])
#     s5_rat2 = np.append(s5_rat2,rat2[k])
#     s5_tns1 = np.append(s5_tns1,tns1[k])
#     s5_tns2 = np.append(s5_tns2,tns2[k])
#     s5_scl1 = np.append(s5_scl1,scl1[k])
#     s5_scl2 = np.append(s5_scl2,scl2[k])
#     s5_scl3 = np.append(s5_scl3,scl3[k])
#     k += 1

## Writing data to csv file
# basepath = path.dirname(path.realpath(__file__))
# fdata = basepath + path.sep + 'vortdatabase_new.csv'
fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/data/velodatabase_SMG_opt11.csv'

with open(fdata,'w') as fp:
    a = csv.writer(fp)
    
    data = np.array(['TSR'])
    data = np.append(data,tsr_write.astype(np.str))
    
    for i in range(np.size(solidity_write)):
        sol = str(i+1)
        
        sollab = 'solidity'
        
        exec('solrow = np.array([sollab,str(solidity_write[i])])')
        for j in range(np.size(tsr)-1):
            solrow = np.append(solrow,'')
        data = np.vstack([data,solrow])
        
        coef = np.zeros((14,np.size(tsr)+1)).astype(np.str)
        coef[0,0] = 'men1'
        coef[1,0] = 'men2'
        coef[2,0] = 'men3'
        coef[3,0] = 'spr1'
        coef[4,0] = 'spr2'
        coef[5,0] = 'spr3'
        coef[6,0] = 'spr4'
        coef[7,0] = 'scl1'
        coef[8,0] = 'scl2'
        coef[9,0] = 'scl3'
        coef[10,0] = 'rat1'
        coef[11,0] = 'rat2'
        coef[12,0] = 'tns1'
        coef[13,0] = 'tns2'
        
        for j in range(np.size(tsr)):
            exec('coef[0,j+1] = s'+sol+'_men1[j]')
            exec('coef[1,j+1] = s'+sol+'_men2[j]')
            exec('coef[2,j+1] = s'+sol+'_men3[j]')
            exec('coef[3,j+1] = s'+sol+'_spr1[j]')
            exec('coef[4,j+1] = s'+sol+'_spr2[j]')
            exec('coef[5,j+1] = s'+sol+'_spr3[j]')
            exec('coef[6,j+1] = s'+sol+'_spr4[j]')
            exec('coef[7,j+1] = s'+sol+'_scl1[j]')
            exec('coef[8,j+1] = s'+sol+'_scl2[j]')
            exec('coef[9,j+1] = s'+sol+'_scl3[j]')
            exec('coef[10,j+1] = s'+sol+'_rat1[j]')
            exec('coef[11,j+1] = s'+sol+'_rat2[j]')
            exec('coef[12,j+1] = s'+sol+'_tns1[j]')
            exec('coef[13,j+1] = s'+sol+'_tns2[j]')
            
                
            
        data = np.vstack([data,coef])
        
    a.writerows(data)


## Vorticity Database Output
print '\n****************COPY THESE FILES****************\n'

print '\ns1_men1 = np.array(',s1_men1.tolist(),')'
print 's1_men2 = np.array(',s1_men2.tolist(),')'
print 's1_men3 = np.array(',s1_men3.tolist(),')'
print '\ns1_spr1 = np.array(',s1_spr1.tolist(),')'
print 's1_spr2 = np.array(',s1_spr2.tolist(),')'
print 's1_spr3 = np.array(',s1_spr3.tolist(),')'
print 's1_spr4 = np.array(',s1_spr4.tolist(),')'
print '\ns1_scl1 = np.array(',s1_scl1.tolist(),')'
print 's1_scl2 = np.array(',s1_scl2.tolist(),')'
print 's1_scl3 = np.array(',s1_scl3.tolist(),')'
print '\ns1_rat1 = np.array(',s1_rat1.tolist(),')'
print 's1_rat2 = np.array(',s1_rat2.tolist(),')'
print '\ns1_tns1 = np.array(',s1_tns1.tolist(),')'
print 's1_tns2 = np.array(',s1_tns2.tolist(),')'

# print '\ns2_men1 = np.array(',s2_men1.tolist(),')'
# print 's2_men2 = np.array(',s2_men2.tolist(),')'
# print 's2_men3 = np.array(',s2_men3.tolist(),')'
# print '\ns2_spr1 = np.array(',s2_spr1.tolist(),')'
# print 's2_spr2 = np.array(',s2_spr2.tolist(),')'
# print 's2_spr3 = np.array(',s2_spr3.tolist(),')'
# print '\ns2_scl1 = np.array(',s2_scl1.tolist(),')'
# print 's2_scl2 = np.array(',s2_scl2.tolist(),')'
# print 's2_scl3 = np.array(',s2_scl3.tolist(),')'
# print '\ns2_rat1 = np.array(',s2_rat1.tolist(),')'
# print 's2_rat2 = np.array(',s2_rat2.tolist(),')'
# print '\ns2_tns1 = np.array(',s2_tns1.tolist(),')'
# print 's2_tns2 = np.array(',s2_tns2.tolist(),')'
# 
# print '\ns3_men1 = np.array(',s3_men1.tolist(),')'
# print 's3_men2 = np.array(',s3_men2.tolist(),')'
# print 's3_men3 = np.array(',s3_men3.tolist(),')'
# print '\ns3_spr1 = np.array(',s3_spr1.tolist(),')'
# print 's3_spr2 = np.array(',s3_spr2.tolist(),')'
# print 's3_spr3 = np.array(',s3_spr3.tolist(),')'
# print '\ns3_scl1 = np.array(',s3_scl1.tolist(),')'
# print 's3_scl2 = np.array(',s3_scl2.tolist(),')'
# print 's3_scl3 = np.array(',s3_scl3.tolist(),')'
# print '\ns3_rat1 = np.array(',s3_rat1.tolist(),')'
# print 's3_rat2 = np.array(',s3_rat2.tolist(),')'
# print '\ns3_tns1 = np.array(',s3_tns1.tolist(),')'
# print 's3_tns2 = np.array(',s3_tns2.tolist(),')'
# 
# print '\ns4_men1 = np.array(',s4_men1.tolist(),')'
# print 's4_men2 = np.array(',s4_men2.tolist(),')'
# print 's4_men3 = np.array(',s4_men3.tolist(),')'
# print '\ns4_spr1 = np.array(',s4_spr1.tolist(),')'
# print 's4_spr2 = np.array(',s4_spr2.tolist(),')'
# print 's4_spr3 = np.array(',s4_spr3.tolist(),')'
# print '\ns4_scl1 = np.array(',s4_scl1.tolist(),')'
# print 's4_scl2 = np.array(',s4_scl2.tolist(),')'
# print 's4_scl3 = np.array(',s4_scl3.tolist(),')'
# print '\ns4_rat1 = np.array(',s4_rat1.tolist(),')'
# print 's4_rat2 = np.array(',s4_rat2.tolist(),')'
# print '\ns4_tns1 = np.array(',s4_tns1.tolist(),')'
# print 's4_tns2 = np.array(',s4_tns2.tolist(),')'
# 
# print '\ns5_men1 = np.array(',s5_men1.tolist(),')'
# print 's5_men2 = np.array(',s5_men2.tolist(),')'
# print 's5_men3 = np.array(',s5_men3.tolist(),')'
# print '\ns5_spr1 = np.array(',s5_spr1.tolist(),')'
# print 's5_spr2 = np.array(',s5_spr2.tolist(),')'
# print 's5_spr3 = np.array(',s5_spr3.tolist(),')'
# print '\ns5_scl1 = np.array(',s5_scl1.tolist(),')'
# print 's5_scl2 = np.array(',s5_scl2.tolist(),')'
# print 's5_scl3 = np.array(',s5_scl3.tolist(),')'
# print '\ns5_rat1 = np.array(',s5_rat1.tolist(),')'
# print 's5_rat2 = np.array(',s5_rat2.tolist(),')'
# print '\ns5_tns1 = np.array(',s5_tns1.tolist(),')'
# print 's5_tns2 = np.array(',s5_tns2.tolist(),')'
