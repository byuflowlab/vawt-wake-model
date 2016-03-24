import numpy as np
import matplotlib.pyplot as plt
from EMG_dist_fit import fit
import csv
from os import path

s = np.array([])
t = np.array([])

solidity = np.array(['s0.15','s0.25','s0.50','s0.75','s1.0'])
tsr = np.linspace(150,700,23)

tsr_write = tsr/100.
solidity_write = np.array([0.15,0.25,0.5,0.75,1.0])

for i in range(np.size(solidity)):
    for j in range(np.size(tsr)):
        s = np.append(s,solidity[i])
        t = np.append(t,str(int(tsr[j])))
        
s1length = np.array([210,210,205,196,185,178,170,165,160,145,140,123,115,112,108,101,101,90,85,80,78,75,70])
s2length = np.array([197,193,185,176,140,146,126,114,103,96,100,86,77,72,70,68,60,64,54,50,47,45,44])
s3length = np.array([185,150,100,95,83,76,72,63,60,49,50,41,39,36,34,33,30,31,28,30,29,28,27])
s4length = np.array([145,100,73,60,53,44,42,37,38,30,33,26,22,24,23,21,21,19,24,23,22,21,20])
s5length = np.array([78,70,52,43,37,32,29,27,26,23,20,20,23,21,20,19,19,18,18,16,16,15,14])

slength = np.array([])
slength = np.append(slength,s1length)
slength = np.append(slength,s2length)
slength = np.append(slength,s3length)
slength = np.append(slength,s4length)
slength = np.append(slength,s5length)
slength = slength*1.

loc1 = np.array([])
loc2 = np.array([])
loc3 = np.array([])
spr1 = np.array([])
spr2 = np.array([])
skw1 = np.array([])
skw2 = np.array([])
scl1 = np.array([])
scl2 = np.array([])
scl3 = np.array([])

for i in range(np.size(s)):
    _,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,_,loc1d,loc2d,loc3d,spr1d,spr2d,skw1d,skw2d,scl1d,scl2d,scl3d = fit(s[i],t[i],slength[i])
    
    loc1 = np.append(loc1,loc1d)
    loc2 = np.append(loc2,loc2d)
    loc3 = np.append(loc3,loc3d)
    spr1 = np.append(spr1,spr1d)
    spr2 = np.append(spr2,spr2d)
    skw1 = np.append(skw1,skw1d)
    skw2 = np.append(skw2,skw2d)
    scl1 = np.append(scl1,scl1d)
    scl2 = np.append(scl2,scl2d)
    scl3 = np.append(scl3,scl3d)

s1_loc1 = np.array([])
s1_loc2 = np.array([])
s1_loc3 = np.array([])
s1_spr1 = np.array([])
s1_spr2 = np.array([])
s1_skw1 = np.array([])
s1_skw2 = np.array([])
s1_scl1 = np.array([])
s1_scl2 = np.array([])
s1_scl3 = np.array([])

s2_loc1 = np.array([])
s2_loc2 = np.array([])
s2_loc3 = np.array([])
s2_spr1 = np.array([])
s2_spr2 = np.array([])
s2_skw1 = np.array([])
s2_skw2 = np.array([])
s2_scl1 = np.array([])
s2_scl2 = np.array([])
s2_scl3 = np.array([])

s3_loc1 = np.array([])
s3_loc2 = np.array([])
s3_loc3 = np.array([])
s3_spr1 = np.array([])
s3_spr2 = np.array([])
s3_skw1 = np.array([])
s3_skw2 = np.array([])
s3_scl1 = np.array([])
s3_scl2 = np.array([])
s3_scl3 = np.array([])

s4_loc1 = np.array([])
s4_loc2 = np.array([])
s4_loc3 = np.array([])
s4_spr1 = np.array([])
s4_spr2 = np.array([])
s4_skw1 = np.array([])
s4_skw2 = np.array([])
s4_scl1 = np.array([])
s4_scl2 = np.array([])
s4_scl3 = np.array([])

s5_loc1 = np.array([])
s5_loc2 = np.array([])
s5_loc3 = np.array([])
s5_spr1 = np.array([])
s5_spr2 = np.array([])
s5_skw1 = np.array([])
s5_skw2 = np.array([])
s5_scl1 = np.array([])
s5_scl2 = np.array([])
s5_scl3 = np.array([])

k = 0
for i in range(np.size(tsr)):
    s1_loc1 = np.append(s1_loc1,loc1[k])
    s1_loc2 = np.append(s1_loc2,loc2[k])
    s1_loc3 = np.append(s1_loc3,loc3[k])
    s1_spr1 = np.append(s1_spr1,spr1[k])
    s1_spr2 = np.append(s1_spr2,spr2[k])
    s1_skw1 = np.append(s1_skw1,skw1[k])
    s1_skw2 = np.append(s1_skw2,skw2[k])
    s1_scl1 = np.append(s1_scl1,scl1[k])
    s1_scl2 = np.append(s1_scl2,scl2[k])
    s1_scl3 = np.append(s1_scl3,scl3[k])
    k += 1
for i in range(np.size(tsr)):
    s2_loc1 = np.append(s2_loc1,loc1[k])
    s2_loc2 = np.append(s2_loc2,loc2[k])
    s2_loc3 = np.append(s2_loc3,loc3[k])
    s2_spr1 = np.append(s2_spr1,spr1[k])
    s2_spr2 = np.append(s2_spr2,spr2[k])
    s2_skw1 = np.append(s2_skw1,skw1[k])
    s2_skw2 = np.append(s2_skw2,skw2[k])
    s2_scl1 = np.append(s2_scl1,scl1[k])
    s2_scl2 = np.append(s2_scl2,scl2[k])
    s2_scl3 = np.append(s2_scl3,scl3[k])
    k += 1
for i in range(np.size(tsr)):
    s3_loc1 = np.append(s3_loc1,loc1[k])
    s3_loc2 = np.append(s3_loc2,loc2[k])
    s3_loc3 = np.append(s3_loc3,loc3[k])
    s3_spr1 = np.append(s3_spr1,spr1[k])
    s3_spr2 = np.append(s3_spr2,spr2[k])
    s3_skw1 = np.append(s3_skw1,skw1[k])
    s3_skw2 = np.append(s3_skw2,skw2[k])
    s3_scl1 = np.append(s3_scl1,scl1[k])
    s3_scl2 = np.append(s3_scl2,scl2[k])
    s3_scl3 = np.append(s3_scl3,scl3[k])
    k += 1
for i in range(np.size(tsr)):
    s4_loc1 = np.append(s4_loc1,loc1[k])
    s4_loc2 = np.append(s4_loc2,loc2[k])
    s4_loc3 = np.append(s4_loc3,loc3[k])
    s4_spr1 = np.append(s4_spr1,spr1[k])
    s4_spr2 = np.append(s4_spr2,spr2[k])
    s4_skw1 = np.append(s4_skw1,skw1[k])
    s4_skw2 = np.append(s4_skw2,skw2[k])
    s4_scl1 = np.append(s4_scl1,scl1[k])
    s4_scl2 = np.append(s4_scl2,scl2[k])
    s4_scl3 = np.append(s4_scl3,scl3[k])
    k += 1
for i in range(np.size(tsr)):
    s5_loc1 = np.append(s5_loc1,loc1[k])
    s5_loc2 = np.append(s5_loc2,loc2[k])
    s5_loc3 = np.append(s5_loc3,loc3[k])
    s5_spr1 = np.append(s5_spr1,spr1[k])
    s5_spr2 = np.append(s5_spr2,spr2[k])
    s5_skw1 = np.append(s5_skw1,skw1[k])
    s5_skw2 = np.append(s5_skw2,skw2[k])
    s5_scl1 = np.append(s5_scl1,scl1[k])
    s5_scl2 = np.append(s5_scl2,scl2[k])
    s5_scl3 = np.append(s5_scl3,scl3[k])
    k += 1

## Writing data to csv file
# basepath = path.dirname(path.realpath(__file__))
# fdata = basepath + path.sep + 'vortdatabase_new.csv'
fdata = 'C:\Users\TingeyPC\Documents\zStar-CCM\NACA0021\MoveForward\\Vorticity Sections\\vortdatabase_new.csv'

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
        
        coef = np.zeros((10,np.size(tsr)+1)).astype(np.str)
        coef[0,0] = 'loc1'
        coef[1,0] = 'loc2'
        coef[2,0] = 'loc3'
        coef[3,0] = 'spr1'
        coef[4,0] = 'spr2'
        coef[5,0] = 'skw1'
        coef[6,0] = 'skw2'
        coef[7,0] = 'scl1'
        coef[8,0] = 'scl2'
        coef[9,0] = 'scl3'
        
        for j in range(np.size(tsr)):
            exec('coef[0,j+1] = s'+sol+'_loc1[j]')
            exec('coef[1,j+1] = s'+sol+'_loc2[j]')
            exec('coef[2,j+1] = s'+sol+'_loc3[j]')
            exec('coef[3,j+1] = s'+sol+'_spr1[j]')
            exec('coef[4,j+1] = s'+sol+'_spr2[j]')
            exec('coef[5,j+1] = s'+sol+'_skw1[j]')
            exec('coef[6,j+1] = s'+sol+'_skw2[j]')
            exec('coef[7,j+1] = s'+sol+'_scl1[j]')
            exec('coef[8,j+1] = s'+sol+'_scl2[j]')
            exec('coef[9,j+1] = s'+sol+'_scl3[j]')
                
            
        data = np.vstack([data,coef])
        
    a.writerows(data)


## Vorticity Database Output
print '\n****************COPY THESE FILES****************\n'

print '\ns1_loc1 = np.array(',s1_loc1.tolist(),')'
print 's1_loc2 = np.array(',s1_loc2.tolist(),')'
print 's1_loc3 = np.array(',s1_loc3.tolist(),')'
print '\ns1_spr1 = np.array(',s1_spr1.tolist(),')'
print 's1_spr2 = np.array(',s1_spr2.tolist(),')'
print '\ns1_skw1 = np.array(',s1_skw1.tolist(),')'
print 's1_skw2 = np.array(',s1_skw2.tolist(),')'
print '\ns1_scl1 = np.array(',s1_scl1.tolist(),')'
print 's1_scl2 = np.array(',s1_scl2.tolist(),')'
print 's1_scl3 = np.array(',s1_scl3.tolist(),')'

print '\ns2_loc1 = np.array(',s2_loc1.tolist(),')'
print 's2_loc2 = np.array(',s2_loc2.tolist(),')'
print 's2_loc3 = np.array(',s2_loc3.tolist(),')'
print '\ns2_spr1 = np.array(',s2_spr1.tolist(),')'
print 's2_spr2 = np.array(',s2_spr2.tolist(),')'
print '\ns2_skw1 = np.array(',s2_skw1.tolist(),')'
print 's2_skw2 = np.array(',s2_skw2.tolist(),')'
print '\ns2_scl1 = np.array(',s2_scl1.tolist(),')'
print 's2_scl2 = np.array(',s2_scl2.tolist(),')'
print 's2_scl3 = np.array(',s2_scl3.tolist(),')'

print '\ns3_loc1 = np.array(',s3_loc1.tolist(),')'
print 's3_loc2 = np.array(',s3_loc2.tolist(),')'
print 's3_loc3 = np.array(',s3_loc3.tolist(),')'
print '\ns3_spr1 = np.array(',s3_spr1.tolist(),')'
print 's3_spr2 = np.array(',s3_spr2.tolist(),')'
print '\ns3_skw1 = np.array(',s3_skw1.tolist(),')'
print 's3_skw2 = np.array(',s3_skw2.tolist(),')'
print '\ns3_scl1 = np.array(',s3_scl1.tolist(),')'
print 's3_scl2 = np.array(',s3_scl2.tolist(),')'
print 's3_scl3 = np.array(',s3_scl3.tolist(),')'

print '\ns4_loc1 = np.array(',s4_loc1.tolist(),')'
print 's4_loc2 = np.array(',s4_loc2.tolist(),')'
print 's4_loc3 = np.array(',s4_loc3.tolist(),')'
print '\ns4_spr1 = np.array(',s4_spr1.tolist(),')'
print 's4_spr2 = np.array(',s4_spr2.tolist(),')'
print '\ns4_skw1 = np.array(',s4_skw1.tolist(),')'
print 's4_skw2 = np.array(',s4_skw2.tolist(),')'
print '\ns4_scl1 = np.array(',s4_scl1.tolist(),')'
print 's4_scl2 = np.array(',s4_scl2.tolist(),')'
print 's4_scl3 = np.array(',s4_scl3.tolist(),')'

print '\ns5_loc1 = np.array(',s5_loc1.tolist(),')'
print 's5_loc2 = np.array(',s5_loc2.tolist(),')'
print 's5_loc3 = np.array(',s5_loc3.tolist(),')'
print '\ns5_spr1 = np.array(',s5_spr1.tolist(),')'
print 's5_spr2 = np.array(',s5_spr2.tolist(),')'
print '\ns5_skw1 = np.array(',s5_skw1.tolist(),')'
print 's5_skw2 = np.array(',s5_skw2.tolist(),')'
print '\ns5_scl1 = np.array(',s5_scl1.tolist(),')'
print 's5_scl2 = np.array(',s5_scl2.tolist(),')'
print 's5_scl3 = np.array(',s5_scl3.tolist(),')'
