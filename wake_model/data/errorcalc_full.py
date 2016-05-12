import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from VAWT_Wake_Model import velocity_field
from database_call import vorticity,velocity
import csv

# Enter the values desired

velf = 15.0 # free stream wind speed (m/s)
dia = 6.0  # turbine diameter (m)
rad = dia/2.
# tsr = 4.0  # tip speed ratio
B = 3. # number of blades
chord = 0.25 # chord lenth (m)
solidity = (chord*B)/(dia/2.)

# Enter the positions of the turbine and velocity calculation
xt = 0. # downstream position of turbine (m)
yt = 0. # later position of turbine (m)
x0 = 5. # downstream position of velocity calculation from turbine (m)
y0 = 0. # lateral position of velocity calculation from turbine (m)

# Choose whether CFD vorticity or velocity data will be used as the basis
cfd_data = 'vort'
cfd_data = 'velo'

s = np.array([])
t = np.array([])
sole = np.array([])
tsre = np.array([])

solidityf = np.array(['s1','s2','s3','s4','s5'])
solf = np.array([0.15,0.25,0.5,0.75,1.0])
tsrf = np.linspace(150,700,23)

tsr_write = tsrf/100.
solidity_write = np.array([0.15,0.25,0.5,0.75,1.0])

for i in range(np.size(solidityf)):
    for j in range(np.size(tsrf)):
        s = np.append(s,solidityf[i])
        t = np.append(t,str(int(tsrf[j])))
        sole = np.append(sole,solf[i])
        tsre = np.append(tsre,tsrf[j]/100.)

error = np.zeros((7,np.size(s)))
std = np.zeros((7,np.size(s)))

for k in range(np.size(s)):
    
    # if cfd_data == 'vort':
    #     wfit = s[k]+'_'+t[k]
    # elif cfd_data == 'velo':
    wfit = s[k]+'_'+t[k]+'.0'
    
    tsr = tsre[k]
    solidity = sole[k]

    if cfd_data == 'vort':
        loc,spr,skw,scl = vorticity(tsr,solidity)
        param = np.array([loc,spr,skw,scl])
    elif cfd_data == 'velo':
        men1,sdv1,rat1,tns1,spr1,scl1,tsrn1,_ = velocity(tsr-0.1249,solidity)
        men2,sdv2,rat2,tns2,spr2,scl2,tsrn2,_ = velocity(tsr+0.1249,solidity)
        if solidity >= 0.35:
            men3,sdv3,rat3,tns3,spr3,scl3,_,soln1 = velocity(tsr,solidity-0.1249)
            men4,sdv4,rat4,tns4,spr4,scl4,_,soln2 = velocity(tsr,solidity+0.1249)
        elif solidity >=0.25:
            men3,sdv3,rat3,tns3,spr3,scl3,_,soln1 = velocity(tsr,solidity-0.049)
            men4,sdv4,rat4,tns4,spr4,scl4,_,soln2 = velocity(tsr,solidity+0.1249)
        else:
            men3,sdv3,rat3,tns3,spr3,scl3,_,soln1 = velocity(tsr,solidity-0.049)
            men4,sdv4,rat4,tns4,spr4,scl4,_,soln2 = velocity(tsr,solidity+0.049)
        if tsrn1 == tsrn2:
            p = 0.
        else:
            p = (tsr-tsrn1)/(tsrn2-tsrn1)
        if soln1 == soln2:
            q = 0.
        else:
            q = (sol-soln1)/(soln2-soln1)

        param = np.array([men1,sdv1,rat1,tns1,spr1,scl1,men2,sdv2,rat2,tns2,spr2,scl2,men3,sdv3,rat3,tns3,spr3,scl3,men4,sdv4,rat4,tns4,spr4,scl4,p,q])

    

    ## Reading in STAR-CCM+ File
    


# def starccm_read(fdata,length,rad,wind,rot,tsr,s):
    fdata = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/Velocity Sections/Dia_test_'+wfit+'.csv'
    
    for i in range(30):
        name = str(i+1)
        exec('pos'+name+' = np.array([])')
        exec('velo'+name+' = np.array([])')

    f = open(fdata)

    csv_f = csv.reader(f)

    i = 0
    for row in csv_f:
        if row[0] != 'null' and i != 0:
            pos1 = np.append(pos1,float(row[0]))
        if row[2] != 'null' and i != 0:
            pos2 = np.append(pos2,float(row[2]))
        if row[4] != 'null' and i != 0:
            pos3 = np.append(pos3,float(row[4]))
        if row[6] != 'null' and i != 0:
            pos4 = np.append(pos4,float(row[6]))
        if row[8] != 'null' and i != 0:
            pos5 = np.append(pos5,float(row[8]))
        if row[10] != 'null' and i != 0:
            pos6 = np.append(pos6,float(row[10]))
        if row[12] != 'null' and i != 0:
            pos7 = np.append(pos7,float(row[12]))
        if row[14] != 'null' and i != 0:
            pos8 = np.append(pos8,float(row[14]))
        if row[16] != 'null' and i != 0:
            pos9 = np.append(pos9,float(row[16]))
        if row[18] != 'null' and i != 0:
            pos10 = np.append(pos10,float(row[18]))
        if row[20] != 'null' and i != 0:
            pos11 = np.append(pos11,float(row[20]))
        if row[22] != 'null' and i != 0:
            pos12 = np.append(pos12,float(row[22]))
        if row[24] != 'null' and i != 0:
            pos13 = np.append(pos13,float(row[24]))
        if row[26] != 'null' and i != 0:
            pos14 = np.append(pos14,float(row[26]))
        if row[28] != 'null' and i != 0:
            pos15 = np.append(pos15,float(row[28]))
        if row[30] != 'null' and i != 0:
            pos16 = np.append(pos16,float(row[30]))
        if row[32] != 'null' and i != 0:
            pos17 = np.append(pos17,float(row[32]))
        if row[34] != 'null' and i != 0:
            pos18 = np.append(pos18,float(row[34]))
        if row[36] != 'null' and i != 0:
            pos19 = np.append(pos19,float(row[36]))
        if row[38] != 'null' and i != 0:
            pos20 = np.append(pos20,float(row[38]))
        if row[40] != 'null' and i != 0:
            pos21 = np.append(pos21,float(row[40]))
        if row[42] != 'null' and i != 0:
            pos22 = np.append(pos22,float(row[42]))
        if row[44] != 'null' and i != 0:
            pos23 = np.append(pos23,float(row[44]))
        if row[46] != 'null' and i != 0:
            pos24 = np.append(pos24,float(row[46]))
        if row[48] != 'null' and i != 0:
            pos25 = np.append(pos25,float(row[48]))
        if row[50] != 'null' and i != 0:
            pos26 = np.append(pos26,float(row[50]))
        if row[52] != 'null' and i != 0:
            pos27 = np.append(pos27,float(row[52]))
        if row[54] != 'null' and i != 0:
            pos28 = np.append(pos28,float(row[54]))
        if row[56] != 'null' and i != 0:
            pos29 = np.append(pos29,float(row[56]))
        if row[58] != 'null' and i != 0:
            pos30 = np.append(pos30,float(row[58]))
            
        if row[1] != 'null' and i != 0:
            velo1 = np.append(velo1,float(row[1]))
        if row[3] != 'null' and i != 0:
            velo2 = np.append(velo2,float(row[3]))
        if row[5] != 'null' and i != 0:
            velo3 = np.append(velo3,float(row[5]))
        if row[7] != 'null' and i != 0:
            velo4 = np.append(velo4,float(row[7]))
        if row[9] != 'null' and i != 0:
            velo5 = np.append(velo5,float(row[9]))
        if row[11] != 'null' and i != 0:
            velo6 = np.append(velo6,float(row[11]))
        if row[13] != 'null' and i != 0:
            velo7 = np.append(velo7,float(row[13]))
        if row[15] != 'null' and i != 0:
            velo8 = np.append(velo8,float(row[15]))
        if row[17] != 'null' and i != 0:
            velo9 = np.append(velo9,float(row[17]))
        if row[19] != 'null' and i != 0:
            velo10 = np.append(velo10,float(row[19]))
        if row[21] != 'null' and i != 0:
            velo11 = np.append(velo11,float(row[21]))
        if row[23] != 'null' and i != 0:
            velo12 = np.append(velo12,float(row[23]))
        if row[25] != 'null' and i != 0:
            velo13 = np.append(velo13,float(row[25]))
        if row[27] != 'null' and i != 0:
            velo14 = np.append(velo14,float(row[27]))
        if row[29] != 'null' and i != 0:
            velo15 = np.append(velo15,float(row[29]))
        if row[31] != 'null' and i != 0:
            velo16 = np.append(velo16,float(row[31]))
        if row[33] != 'null' and i != 0:
            velo17 = np.append(velo17,float(row[33]))
        if row[35] != 'null' and i != 0:
            velo18 = np.append(velo18,float(row[35]))
        if row[37] != 'null' and i != 0:
            velo19 = np.append(velo19,float(row[37]))
        if row[39] != 'null' and i != 0:
            velo20 = np.append(velo20,float(row[39]))
        if row[41] != 'null' and i != 0:
            velo21 = np.append(velo21,float(row[41]))
        if row[43] != 'null' and i != 0:
            velo22 = np.append(velo22,float(row[43]))
        if row[45] != 'null' and i != 0:
            velo23 = np.append(velo23,float(row[45]))
        if row[47] != 'null' and i != 0:
            velo24 = np.append(velo24,float(row[47]))
        if row[49] != 'null' and i != 0:
            velo25 = np.append(velo25,float(row[49]))
        if row[51] != 'null' and i != 0:
            velo26 = np.append(velo26,float(row[51]))
        if row[53] != 'null' and i != 0:
            velo27 = np.append(velo27,float(row[53]))
        if row[55] != 'null' and i != 0:
            velo28 = np.append(velo28,float(row[55]))
        if row[57] != 'null' and i != 0:
            velo29 = np.append(velo29,float(row[57]))
        if row[59] != 'null' and i != 0:
            velo30 = np.append(velo30,float(row[59]))
        i += 1

    f.close()
    # print 'DATA IMPORTED'

    wind = velf
    for i in range(30):
        name = str(i+1)

        #Ordering the data numerically by the position
        exec('pos'+name+', velo'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', velo'+name+'))))')
        #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
        exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvelo'+name+' = np.copy(velo'+name+'_0)/velf')
        #Deleting wall boundary data
        exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5.*dia or pos'+name+'[j] < -5.*dia:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')
        
        # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif velo'+name+'[j] > 1.*wind or fabs(pos'+name+'[j]) > 1.2*dia:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')




    # start = length/30.
    # x = np.linspace(start,length,30)





## Error
    error2 = np.zeros_like(pos2)
    error4 = np.zeros_like(pos4)
    error6 = np.zeros_like(pos6)
    error8 = np.zeros_like(pos8)
    error10 = np.zeros_like(pos10)
    error15 = np.zeros_like(pos15)
    error20 = np.zeros_like(pos20)
    
    tot2 = np.size(pos2)
    tot4 = np.size(pos4)
    tot6 = np.size(pos6)
    tot8 = np.size(pos8)
    tot10 = np.size(pos10)
    tot15 = np.size(pos15)
    tot20 = np.size(pos20)
    
    for i in range(np.size(pos2)):
        velp = velocity_field(xt,yt,2*dia,pos2[i],velf,dia,tsr,solidity,cfd_data,param)
        error2[i] = (velp-velo2[i])/velo2[i]
        print '2D',i,'of',tot2,'TSR:',tsr,'Solidity:',solidity
    for i in range(np.size(pos4)):
        velp = velocity_field(xt,yt,4*dia,pos4[i],velf,dia,tsr,solidity,cfd_data,param)
        error4[i] = (velp-velo4[i])/velo4[i]
        print '4D',i,'of',tot4,'TSR:',tsr,'Solidity:',solidity
    for i in range(np.size(pos6)):
        velp = velocity_field(xt,yt,6*dia,pos6[i],velf,dia,tsr,solidity,cfd_data,param)
        error6[i] = (velp-velo6[i])/velo6[i]
        print '6D',i,'of',tot6,'TSR:',tsr,'Solidity:',solidity
    for i in range(np.size(pos8)):
        velp = velocity_field(xt,yt,8*dia,pos8[i],velf,dia,tsr,solidity,cfd_data,param)
        error8[i] = (velp-velo8[i])/velo8[i]
        print '8D',i,'of',tot8,'TSR:',tsr,'Solidity:',solidity
    for i in range(np.size(pos10)):
        velp = velocity_field(xt,yt,10*dia,pos10[i],velf,dia,tsr,solidity,cfd_data,param)
        error10[i] = (velp-velo10[i])/velo10[i]
        print '10D',i,'of',tot10,'TSR:',tsr,'Solidity:',solidity
    for i in range(np.size(pos15)):
        velp = velocity_field(xt,yt,15*dia,pos15[i],velf,dia,tsr,solidity,cfd_data,param)
        error15[i] = (velp-velo15[i])/velo15[i]
        print '15D',i,'of',tot15,'TSR:',tsr,'Solidity:',solidity
    for i in range(np.size(pos20)):
        velp = velocity_field(xt,yt,20*dia,pos20[i],velf,dia,tsr,solidity,cfd_data,param)
        error20[i] = (velp-velo20[i])/velo20[i]
        print '20D',i,'of',tot20,'TSR:',tsr,'Solidity:',solidity
    
    error2m = np.average(error2)
    error2std = np.std(error2)
    error4m = np.average(error4)
    error4std = np.std(error4)
    error6m = np.average(error6)
    error6std = np.std(error6)
    error8m = np.average(error8)
    error8std = np.std(error8)
    error10m = np.average(error10)
    error10std = np.std(error10)
    error15m = np.average(error15)
    error15std = np.std(error15)
    error20m = np.average(error20)
    error20std = np.std(error20)
    
    error[0,k] = error2m
    error[1,k] = error4m
    error[2,k] = error6m
    error[3,k] = error8m
    error[4,k] = error10m
    error[5,k] = error15m
    error[6,k] = error20m
    std[0,k] = error2std
    std[1,k] = error4std
    std[2,k] = error6std
    std[3,k] = error8std
    std[4,k] = error10std
    std[5,k] = error15std
    std[6,k] = error20std
    
    
    # print 'x = 2D',error2m,error2std
    # print 'x = 4D',error4m,error4std
    # print 'x = 6D',error6m,error6std
    # print 'x = 8D',error8m,error8std
    # print 'x = 10D',error10m,error10std
    # print 'x = 15D',error15m,error15std
    # print 'x = 20D',error20m,error20std
    
## Write
fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/validation/error_cfd_SMG_surf.csv'

q = 0
with open(fdata,'w') as fp:
    a = csv.writer(fp)
    
    data = np.array(['TSR'])
    data = np.append(data,tsr_write.astype(np.str))
    
    for i in range(np.size(solidity_write)):
        sol = str(i+1)
        
        sollab = 'solidity'
        
        solrow = np.array([sollab,str(solidity_write[i])])
        for j in range(np.size(tsrf)-1):
            solrow = np.append(solrow,'')
        data = np.vstack([data,solrow])
        
        errorval = np.zeros((14,np.size(tsrf)+1)).astype(np.str)
        errorval[0,0] = '2D_error'
        errorval[7,0] = '2D_std'
        errorval[1,0] = '4D_error'
        errorval[8,0] = '4D_std'
        errorval[2,0] = '6D_error'
        errorval[9,0] = '6D_std'
        errorval[3,0] = '8D_error'
        errorval[10,0] = '8D_std'
        errorval[4,0] = '10D_error'
        errorval[11,0] = '10D_std'
        errorval[5,0] = '15D_error'
        errorval[12,0] = '15D_std'
        errorval[6,0] = '20D_error'
        errorval[13,0] = '20D_std'
        
        for j in range(np.size(tsrf)):
            errorval[0,j+1] = error[0,q]
            errorval[1,j+1] = error[1,q]
            errorval[2,j+1] = error[2,q]
            errorval[3,j+1] = error[3,q]
            errorval[4,j+1] = error[4,q]
            errorval[5,j+1] = error[5,q]
            errorval[6,j+1] = error[6,q]
            errorval[7,j+1] = std[0,q]
            errorval[8,j+1] = std[1,q]
            errorval[9,j+1] = std[2,q]
            errorval[10,j+1] = std[3,q]
            errorval[11,j+1] = std[4,q]
            errorval[12,j+1] = std[5,q]
            errorval[13,j+1] = std[6,q]
            q += 1
            
        data = np.vstack([data,errorval])
    
    a.writerows(data)





