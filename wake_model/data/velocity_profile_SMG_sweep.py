import numpy as np
import matplotlib.pyplot as plt
# from velocity_profile_read_skew_def_quad import fit
from velocity_optimize_SMG_surf import fit
import csv
from sys import argv

comp = 'mac'
# comp = 'fsl'

if comp == 'fsl':
    solrun = int(argv[1])

    print 'solrun_num',solrun

    if solrun == 1:
        solidity = np.array(['s1'])#,'s2','s3','s4','s5'])
        solidity_write = np.array([0.15])
        s1length = np.array([210,210,205,196,185,178,170,165,160,145,140,123,115,112,108,101,101,90,85,80,78,75,70])
    elif solrun == 2:
        solidity = np.array(['s2'])
        solidity_write = np.array([0.25])
        s1length = np.array([197,193,185,176,140,146,126,114,103,96,100,86,77,72,70,68,60,64,54,50,47,45,44])
    elif solrun == 3:
        solidity = np.array(['s3'])
        solidity_write = np.array([0.5])
        s1length = np.array([185,150,100,95,83,76,72,63,60,49,50,41,39,36,34,33,30,31,28,30,29,28,27])
    elif solrun == 4:
        solidity = np.array(['s4'])
        solidity_write = np.array([0.75])
        s1length = np.array([145,100,73,60,53,44,42,37,38,30,33,26,22,24,23,21,21,19,24,23,22,21,20])
    elif solrun == 5:
        solidity = np.array(['s5'])
        solidity_write = np.array([1.0])
        s1length = np.array([78,70,52,43,37,32,29,27,26,23,20,20,23,21,20,19,19,18,18,16,16,15,14])
elif comp == 'mac':
    solidity = np.array(['s1'])#,'s2','s3','s4','s5'])
    solidity_write = np.array([0.15])#,0.25,0.5,0.75,1.0])
    s1length = np.array([210,210,205,196,185,178,170,165,160,145,140,123,115,112,108,101,101,90,85,80,78,75,70])
    # s1length = np.array([197,193,185,176,140,146,126,114,103,96,100,86,77,72,70,68,60,64,54,50,47,45,44])
    # s1length = np.array([185,150,100,95,83,76,72,63,60,49,50,41,39,36,34,33,30,31,28,30,29,28,27])
    # s1length = np.array([145,100,73,60,53,44,42,37,38,30,33,26,22,24,23,21,21,19,24,23,22,21,20])
    # s1length = np.array([78,70,52,43,37,32,29,27,26,23,20,20,23,21,20,19,19,18,18,16,16,15,14])

tsr = np.linspace(150,700,23)
tsr_write = tsr/100.

s = np.array([])
t = np.array([])
for i in range(np.size(solidity)):
    for j in range(np.size(tsr)):
        s = np.append(s,solidity[i])
        t = np.append(t,str(int(tsr[j])))

slength = np.array([])
slength = np.append(slength,s1length)
# slength = np.append(slength,s2length)
# slength = np.append(slength,s3length)
# slength = np.append(slength,s4length)
# slength = np.append(slength,s5length)
slength = slength*1.

men = np.array([])
sdv1 = np.array([])
sdv2 = np.array([])
sdv3 = np.array([])
sdv4 = np.array([])
rat = np.array([])
tns = np.array([])
spr1 = np.array([])
spr2 = np.array([])
spr3 = np.array([])
spr4 = np.array([])
scl1 = np.array([])
scl2 = np.array([])
scl3 = np.array([])

for i in range(np.size(s)):
    mend,sdv1d,sdv2d,sdv3d,sdv4d,ratd,tnsd,spr1d,spr2d,spr3d,spr4d,scl1d,scl2d,scl3d = fit(s[i],t[i],slength[i],False,comp,4,False)
    
    men = np.append(men,mend)
    sdv1 = np.append(sdv1,sdv1d)
    sdv2 = np.append(sdv2,sdv2d)
    sdv3 = np.append(sdv3,sdv3d)
    sdv4 = np.append(sdv4,sdv4d)
    rat = np.append(rat,ratd)
    tns = np.append(tns,tnsd)
    spr1 = np.append(spr1,spr1d)
    spr2 = np.append(spr2,spr2d)
    spr3 = np.append(spr3,spr3d)
    spr4 = np.append(spr4,spr4d)
    scl1 = np.append(scl1,scl1d)
    scl2 = np.append(scl2,scl2d)
    scl3 = np.append(scl3,scl3d)
    print s[i],t[i]

s_men = np.array([])
s_sdv1 = np.array([])
s_sdv2 = np.array([])
s_sdv3 = np.array([])
s_sdv4 = np.array([])
s_rat = np.array([])
s_tns = np.array([])
s_spr1 = np.array([])
s_spr2 = np.array([])
s_spr3 = np.array([])
s_spr4 = np.array([])
s_scl1 = np.array([])
s_scl2 = np.array([])
s_scl3 = np.array([])

k = 0
for i in range(np.size(tsr)):
    s_men = np.append(s_men,men)
    s_sdv1 = np.append(s_sdv1,sdv1)
    s_sdv2 = np.append(s_sdv2,sdv2)
    s_sdv3 = np.append(s_sdv3,sdv3)
    s_sdv4 = np.append(s_sdv4,sdv4)
    s_rat = np.append(s_rat,rat)
    s_tns = np.append(s_tns,tns)
    s_spr1 = np.append(s_spr1,spr1)
    s_spr2 = np.append(s_spr2,spr2)
    s_spr3 = np.append(s_spr3,spr3)
    s_spr4 = np.append(s_spr4,spr4)
    s_scl1 = np.append(s_scl1,scl1)
    s_scl2 = np.append(s_scl2,scl2)
    s_scl3 = np.append(s_scl3,scl3)
    k += 1

## Writing data to csv file
# basepath = path.dirname(path.realpath(__file__))
# fdata = basepath + path.sep + 'vortdatabase_new.csv'
if comp == 'mac':
    fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/data/velodatabase_SMG_surf2_1.csv'
elif comp == 'fsl':
    if solrun == 1:
        fdata = '/fslhome/ebtingey/compute/VAWTWakeModel/data/velodatabase_SMG_surf2_1.csv'
    elif solrun == 2:
        fdata = '/fslhome/ebtingey/compute/VAWTWakeModel/data/velodatabase_SMG_surf2_2.csv'
    elif solrun == 3:
        fdata = '/fslhome/ebtingey/compute/VAWTWakeModel/data/velodatabase_SMG_surf2_3.csv'
    elif solrun == 4:
        fdata = '/fslhome/ebtingey/compute/VAWTWakeModel/data/velodatabase_SMG_surf2_4.csv'
    elif solrun == 5:
        fdata = '/fslhome/ebtingey/compute/VAWTWakeModel/data/velodatabase_SMG_surf2_5.csv'

with open(fdata,'w') as fp:
    a = csv.writer(fp)
    
    data = np.array(['TSR'])
    data = np.append(data,tsr_write.astype(np.str))

    sollab = 'solidity'

    solrow = np.array([sollab,str(solidity_write[0])])
    for j in range(np.size(tsr)-1):
        solrow = np.append(solrow,'')
    data = np.vstack([data,solrow])

    coef = np.zeros((14,np.size(tsr)+1)).astype(np.str)
    coef[0,0] = 'men'
    coef[1,0] = 'sdv1'
    coef[2,0] = 'sdv2'
    coef[3,0] = 'sdv3'
    coef[4,0] = 'sdv4'
    coef[5,0] = 'rat'
    coef[6,0] = 'tns'
    coef[7,0] = 'spr1'
    coef[8,0] = 'spr2'
    coef[9,0] = 'spr3'
    coef[10,0] = 'spr4'
    coef[11,0] = 'scl1'
    coef[12,0] = 'scl2'
    coef[13,0] = 'scl3'

    for j in range(np.size(tsr)):
        coef[0,j+1] = s_men[j]
        coef[1,j+1] = s_sdv1[j]
        coef[2,j+1] = s_sdv2[j]
        coef[3,j+1] = s_sdv3[j]
        coef[4,j+1] = s_sdv4[j]
        coef[5,j+1] = s_rat[j]
        coef[6,j+1] = s_tns[j]
        coef[7,j+1] = s_spr1[j]
        coef[8,j+1] = s_spr2[j]
        coef[9,j+1] = s_spr3[j]
        coef[10,j+1] = s_spr4[j]
        coef[11,j+1] = s_scl1[j]
        coef[12,j+1] = s_scl2[j]
        coef[13,j+1] = s_scl3[j]

    data = np.vstack([data,coef])
        
    a.writerows(data)


## Vorticity Database Output
print '\n****************COPY THESE FILES****************\n'

print '\n'
print 's_men =',s_men
print 's_sdv1 =',s_sdv1
print 's_sdv2 =',s_sdv2
print 's_sdv3 =',s_sdv3
print 's_sdv4 =',s_sdv4
print 's_rat =',s_rat
print 's_tns =',s_tns
print 's_spr1 =',s_spr1
print 's_spr2 =',s_spr2
print 's_spr3 =',s_spr3
print 's_spr4 =',s_spr4
print 's_scl1 =',s_scl1
print 's_scl2 =',s_scl2
print 's_scl3 =',s_scl3
