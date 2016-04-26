import numpy as np
import matplotlib.pyplot as plt
# from velocity_profile_read_skew_def_quad import fit
from velocity_optimize_SMG_surf import fit
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

men = np.array([])
sdv = np.array([])
rat = np.array([])
spr = np.array([])
scl1 = np.array([])
scl2 = np.array([])
scl3 = np.array([])

for i in range(np.size(s)):
    mend,sdvd,ratd,sprd,scl1d,scl2d,scl3d = fit(s[i],t[i],slength[i],False,'mac',4,False)
    
    men = np.append(men,mend)
    sdv = np.append(sdv,sdvd)
    rat = np.append(rat,ratd)
    spr = np.append(spr,sprd)
    scl1 = np.append(scl1,scl1d)
    scl2 = np.append(scl2,scl2d)
    scl3 = np.append(scl3,scl3d)
    print s[i],t[i]

s_men = np.array([])
s_sdv = np.array([])
s_rat = np.array([])
s_spr = np.array([])
s_scl1 = np.array([])
s_scl2 = np.array([])
s_scl3 = np.array([])

k = 0
for i in range(np.size(tsr)):
    s_men = np.append([s_men,men])
    s_sdv = np.append([s_sdv,sdv])
    s_rat = np.append([s_rat,rat])
    s_spr = np.append([s_spr,spr])
    s_scl1 = np.append([s_scl1,scl1])
    s_scl2 = np.append([s_scl2,scl2])
    s_scl3 = np.append([s_scl3,scl3])
    k += 1

## Writing data to csv file
# basepath = path.dirname(path.realpath(__file__))
# fdata = basepath + path.sep + 'vortdatabase_new.csv'
fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/data/velodatabase_SMG_surf1.csv'

with open(fdata,'w') as fp:
    a = csv.writer(fp)
    
    data = np.array(['TSR'])
    data = np.append(data,tsr_write.astype(np.str))

    sollab = 'solidity'

    solrow = np.array([sollab,str(solidity_write[i])])
    for j in range(np.size(tsr)-1):
        solrow = np.append(solrow,'')
    data = np.vstack([data,solrow])

    coef = np.zeros((7,np.size(tsr)+1)).astype(np.str)
    coef[0,0] = 'men1'
    coef[1,0] = 'sdv'
    coef[2,0] = 'rat'
    coef[3,0] = 'spr'
    coef[4,0] = 'scl1'
    coef[5,0] = 'scl2'
    coef[6,0] = 'scl3'

    for j in range(np.size(tsr)):
        coef[0,j+1] = s_men[j]
        coef[1,j+1] = s_sdv[j]
        coef[2,j+1] = s_rat[j]
        coef[3,j+1] = s_spr[j]
        coef[4,j+1] = s_scl1[j]
        coef[5,j+1] = s_scl2[j]
        coef[6,j+1] = s_scl3[j]

    data = np.vstack([data,coef])
        
    a.writerows(data)


## Vorticity Database Output
print '\n****************COPY THESE FILES****************\n'

print '\n'
print 's_men =',s_men
print 's_sdv =',s_sdv
print 's_rat =',s_rat
print 's_spr =',s_spr
print 's_scl1 =',s_scl1
print 's_scl2 =',s_scl2
print 's_scl3 =',s_scl3
