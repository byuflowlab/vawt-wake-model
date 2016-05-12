import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from database_call import vorticity,velocity
from VAWT_Wake_Model import velocity_field
import csv

# Enter the values desired

velf = 15.0 # free stream wind speed (m/s)
dia = 6.0  # turbine diameter (m)
tsr = 4.25  # tip speed ratio
B = 3. # number of blades
chord = 0.75 # chord lenth (m)
solidity = (chord*B)/(dia/2.)

# Enter the positions of the turbine and velocity calculation
xt = 0. # downstream position of turbine (m)
yt = 0. # later position of turbine (m)
x0 = 5. # downstream position of velocity calculation from turbine (m)
y0 = 0. # lateral position of velocity calculation from turbine (m)

# Choose whether CFD vorticity or velocity data will be used as the basis
cfd_data = 'vort'
cfd_data = 'velo'

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
wfit = 's4_425.0'
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
print 'Data Imported!'

wind = velf
for i in range(30):
    name = str(i+1)

    #Ordering the data numerically by the position
    exec('pos'+name+', velo'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', velo'+name+'))))')
    #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
    exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvelo'+name+' = np.copy(velo'+name+'_0)/velf')
    #Deleting wall boundary data
    exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5.*dia or pos'+name+'[j] < -5.*dia:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')

## Plotting
fs = 15 # font size for plots

# Option to plot velocity profiles
vel_slice = True
vel_slice = False # comment this out if desired on

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
    print '2D',i,'of',tot2
for i in range(np.size(pos4)):
    velp = velocity_field(xt,yt,4*dia,pos4[i],velf,dia,tsr,solidity,cfd_data,param)
    error4[i] = (velp-velo4[i])/velo4[i]
    print '4D',i,'of',tot4
for i in range(np.size(pos6)):
    velp = velocity_field(xt,yt,6*dia,pos6[i],velf,dia,tsr,solidity,cfd_data,param)
    error6[i] = (velp-velo6[i])/velo6[i]
    print '6D',i,'of',tot6
for i in range(np.size(pos8)):
    velp = velocity_field(xt,yt,8*dia,pos8[i],velf,dia,tsr,solidity,cfd_data,param)
    error8[i] = (velp-velo8[i])/velo8[i]
    print '8D',i,'of',tot8
for i in range(np.size(pos10)):
    velp = velocity_field(xt,yt,10*dia,pos10[i],velf,dia,tsr,solidity,cfd_data,param)
    error10[i] = (velp-velo10[i])/velo10[i]
    print '10D',i,'of',tot10
for i in range(np.size(pos15)):
    velp = velocity_field(xt,yt,15*dia,pos15[i],velf,dia,tsr,solidity,cfd_data,param)
    error15[i] = (velp-velo15[i])/velo15[i]
    print '15D',i,'of',tot15
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

print 'x = 2D',error2m,error2std
print 'x = 4D',error4m,error4std
print 'x = 6D',error6m,error6std
print 'x = 8D',error8m,error8std
print 'x = 10D',error10m,error10std
print 'x = 15D',error15m,error15std
print 'x = 20D',error20m,error20std



# Plotting velocity profiles
if vel_slice == True:
    leng = 100 # data points in the velocity profile
    wide = 2.0*dia # width of the profile
    
    d_lab1 = str(wide/dia) # y-axis label
    d_lab2 = str(wide/(2*dia)) # y-axis label
    
    x = np.array([2*dia,4*dia,6*dia,8*dia,10*dia,15*dia]) # plotting at 2D, 4D, 6D, 8D, 10D, and 15D (changeable)
    y = np.linspace(-wide,wide,leng)
    
    color = np.array(['b','c','g','y','r','m']) # identifying six colors to use for differentiation
    
    iterp = 0
    for i in range(int(np.size(x))):
        vel = np.array([])
        val = str(x[i]/dia)
        lab = '$x/D$ = '+val
        for j in range(int(np.size(y))):
            velp = velocity_field(xt,yt,x[i],y[j],velf,dia,tsr,solidity,cfd_data,param)
            vel = np.append(vel,velp)
            iterp += 1
            print 'Vel Slice ('+str(iterp)+' of '+str(leng*np.size(x))+')'
        plt.figure(1)
        plt.plot(vel,y,color[i],label=lab)
    
    tix = np.array([-wide,-wide/2.,0.,wide/2.,wide])
    tixl = np.array([d_lab1,d_lab2,'0.0',d_lab2,d_lab1])
    # plt.legend(loc="upper left",bbox_to_anchor=(1, 1),fontsize=fs) # legend off to the side
    plt.legend(loc=2,fontsize=fs) # legend in the plot
    plt.xlim(0.,1.2)
    plt.ylim(-wide,wide)
    plt.xticks(fontsize=fs)
    plt.yticks(tix,tixl,fontsize=fs)
    plt.xlabel(r'$u/U_\infty$', fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)

    

plt.show()


# x = 2D -0.0390882979871 0.0523644511862
# x = 4D -0.0386888996738 0.0412770950732
# x = 6D -0.0337707026147 0.0293425108821
# x = 8D -0.0331687565618 0.0272483053483
# x = 10D -0.026554754511 0.0179845075491
# x = 15D -0.0161324655252 0.0258669721121


