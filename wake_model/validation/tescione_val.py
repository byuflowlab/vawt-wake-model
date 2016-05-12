import numpy as np
import matplotlib.pyplot as plt
import csv
from VAWT_Wake_Model import velocity_field
from database_call import vorticity,velocity,velocity2
from scipy.io import loadmat
from numpy import fabs

from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

rom = True

r = 0.5
v = 1.

flip = -1.

# fdata = 'path/to/tes.csv' # adjust this for your case
fdata = '/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/Validation/tes.csv'
f = open(fdata)

csv_f = csv.reader(f)

for i in range(6):
    name = str(i+1)
    exec('pos'+name+' = np.array([])')
    exec('vel'+name+' = np.array([])')

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
        
    if row[1] != 'null' and i != 0:
        vel1 = np.append(vel1,float(row[1]))
    if row[3] != 'null' and i != 0:
        vel2 = np.append(vel2,float(row[3]))
    if row[5] != 'null' and i != 0:
        vel3 = np.append(vel3,float(row[5]))
    if row[7] != 'null' and i != 0:
        vel4 = np.append(vel4,float(row[7]))
    if row[9] != 'null' and i != 0:
        vel5 = np.append(vel5,float(row[9]))
    if row[11] != 'null' and i != 0:
        vel6 = np.append(vel6,float(row[11]))
    i += 1

f.close()



for i in range(6):
    name = str(i+1)
    
    exec('pos'+name+' = (1./(r*2.))*pos'+name+'')
    exec('vel'+name+' = (1./v)*vel'+name+'')
    #Ordering the data numerically by the position
    exec('pos'+name+', vel'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', vel'+name+'))))')
    #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
    exec('pos'+name+'_0 = np.array([])\nvel'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvel'+name+'_0 = np.append(vel'+name+'_0,vel'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvel'+name+' = np.copy(vel'+name+'_0)')

# tesdata = loadmat('/path/to/dataH_VAWT.mat') # adjust this for your case
tesdata = loadmat('/Users/ning1/Dropbox/dataH_VAWT.mat')

x15 = np.zeros(33)
x20 = np.zeros(33)
x25 = np.zeros(33)
x30 = np.zeros(33)
x35 = np.zeros(33)
x40 = np.zeros(26)
y15 = np.zeros(33)
y20 = np.zeros(33)
y25 = np.zeros(33)
y30 = np.zeros(33)
y35 = np.zeros(33)
y40 = np.zeros(26)

space = 52.4545
for i in range(33):
    x15[i] = flip*(tesdata['y'][int(i*space),1706])/1000.
    x20[i] = flip*(tesdata['y'][int(i*space),1956])/1000.
    x25[i] = flip*(tesdata['y'][int(i*space),2206])/1000.
    x30[i] = flip*(tesdata['y'][int(i*space),2456])/1000.
    x35[i] = flip*(tesdata['y'][int(i*space),2706])/1000.
    
    y15[i] = (tesdata['u_x'][int(i*space),1706])/9.3
    y20[i] = (tesdata['u_x'][int(i*space),1956])/9.3
    y25[i] = (tesdata['u_x'][int(i*space),2206])/9.3
    y30[i] = (tesdata['u_x'][int(i*space),2456])/9.3
    y35[i] = (tesdata['u_x'][int(i*space),2706])/9.3

for i in range(26):
    x40[i] = flip*(tesdata['y'][int(i*space+7.*space),2956])/1000.
    
    y40[i] = (tesdata['u_x'][int(i*space+7.*space),2956])/9.3
    
index15 = np.zeros_like(x15)
index20 = np.zeros_like(x20)
index25 = np.zeros_like(x25)
index30 = np.zeros_like(x30)
index35 = np.zeros_like(x35)
index40 = np.zeros_like(x40)

for i in range(33):
    indext15 = np.fabs(pos1-x15[i])
    index15[i] = np.argmin(indext15)
    indext20 = np.fabs(pos2-x20[i])
    index20[i] = np.argmin(indext20)
    indext25 = np.fabs(pos3-x25[i])
    index25[i] = np.argmin(indext25)
    indext30 = np.fabs(pos4-x30[i])
    index30[i] = np.argmin(indext30)
    indext35 = np.fabs(pos5-x35[i])
    index35[i] = np.argmin(indext35)
for i in range(26):
    indext40 = np.fabs(pos6-x40[i])
    index40[i] = np.argmin(indext40)

cfd15t = np.zeros(33)
cfd20t = np.zeros(33)
cfd25t = np.zeros(33)
cfd30t = np.zeros(33)
cfd35t = np.zeros(33)
cfd40t = np.zeros(26)

for i in range(33):
    cfd15t[i] = (vel1[index15[i]]-y15[i])/y15[i]
    cfd20t[i] = (vel2[index20[i]]-y20[i])/y20[i]
    cfd25t[i] = (vel3[index25[i]]-y25[i])/y25[i]
    cfd30t[i] = (vel4[index30[i]]-y30[i])/y30[i]
    cfd35t[i] = (vel5[index35[i]]-y35[i])/y35[i]
for i in range(26):
    cfd40t[i] = (vel6[index40[i]]-y40[i])/y40[i]

cfd15error = np.average(fabs(cfd15t))
cfd15errorstd = np.std(cfd15t)
cfd20error = np.average(fabs(cfd20t))
cfd20errorstd = np.std(cfd20t)
cfd25error = np.average(fabs(cfd25t))
cfd25errorstd = np.std(cfd25t)
cfd30error = np.average(fabs(cfd30t))
cfd30errorstd = np.std(cfd30t)
cfd35error = np.average(fabs(cfd35t))
cfd35errorstd = np.std(cfd35t)
cfd40error = np.average(fabs(cfd40t))
cfd40errorstd = np.std(cfd40t)
cfdoaerror = (cfd15error+cfd20error+cfd25error+cfd30error+cfd35error+cfd40error)/6.
cfdoaerrorstd = (cfd15errorstd+cfd20errorstd+cfd25errorstd+cfd30errorstd+cfd35errorstd+cfd40errorstd)/6.

    
## Plot CFD
fig1 = plt.figure(1,figsize=(12.5,6))
fig1.subplots_adjust(left=.05,right=.86,wspace=.36,hspace=.35)
plt.subplot(2,3,1)
plt.plot(x15,y15,'b.')
plt.plot(pos1,vel1,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,1.05,r'$x/D$ = 0.75')
print '1.5 cfd',(min(vel1)-min(y15))/min(y15),cfd15error,cfd15errorstd
plt.subplot(2,3,2)
plt.plot(x20,y20,'b.')
plt.plot(pos2,vel2,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,1.05,r'$x/D$ = 1.0')
print '2.0 cfd',(min(vel2)-min(y20))/min(y20),cfd20error,cfd20errorstd
plt.subplot(2,3,3)
plt.plot(x25,y25,'b.',label='PIV')
plt.plot(pos3,vel3,'r-',label='CFD')
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,1.05,r'$x/D$ = 1.25')
print '2.5 cfd',(min(vel3)-min(y25))/min(y25),cfd25error,cfd25errorstd
plt.legend(loc="upper left", bbox_to_anchor=(1,1))
plt.subplot(2,3,4)
plt.plot(x30,y30,'b.')
plt.plot(pos4,vel4,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,1.05,r'$x/D$ = 1.5')
print '3.0 cfd',(min(vel4)-min(y30))/min(y30),cfd30error,cfd30errorstd
plt.subplot(2,3,5)
plt.plot(x35,y35,'b.')
plt.plot(pos5,vel5,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,1.05,r'$x/D$ = 1.75')
print '3.5 cfd',(min(vel5)-min(y35))/min(y35),cfd35error,cfd35errorstd
plt.subplot(2,3,6)
plt.plot(x40,y40,'b.')
plt.plot(pos6,vel6,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,1.05,r'$x/D$ = 2.0')
print '4.0 cfd',(min(vel6)-min(y40))/min(y40),cfd40error,cfd40errorstd

print cfdoaerror,cfdoaerrorstd


## Plot Model
if rom == True:
    for k in range(1):
        rad = 0.5
        dia = 2*rad
        velf = 9.308422677
        sol = 0.24
        tsr = 4.5
        xt = 0.
        yt = 0.
        # Choose whether CFD vorticity or velocity data will be used as the basis
        if k == 1:
            cfd_data = 'velo'
        elif k == 0:
            cfd_data = 'velo2'
        # cfd_data = 'quad'
        
        if cfd_data == 'vort':
            loc,spr,skw,scl = vorticity(tsr,sol)
            param = np.array([loc,spr,skw,scl])
            
        elif cfd_data == 'velo':
            men1,sdv1,rat1,wdt1,spr1,scl1,tsrn1,_ = velocity(tsr-0.1249,sol)
            men2,sdv2,rat2,wdt2,spr2,scl2,tsrn2,_ = velocity(tsr+0.1249,sol)
            if sol >= 0.35:
                men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,sol-0.1249)
                men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,sol+0.1249)
            elif sol >=0.25:
                men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,sol-0.049)
                men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,sol+0.1249)
            else:
                men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,sol-0.049)
                men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,sol+0.049)
            if tsrn1 == tsrn2:
                p = 0.
            else:
                p = (tsr-tsrn1)/(tsrn2-tsrn1)
            if soln1 == soln2:
                q = 0.
            else:
                q = (sol-soln1)/(soln2-soln1)

            # import time
            # print tsrn1,tsrn2
            # print soln1,soln2
            # print p,q
            # time.sleep(10)
            param = np.array([men1,sdv1,rat1,wdt1,spr1,scl1,men2,sdv2,rat2,wdt2,spr2,scl2,men3,sdv3,rat3,wdt3,spr3,scl3,men4,sdv4,rat4,wdt4,spr4,scl4,p,q])

        elif cfd_data == 'velo2':
            men,sdv,rat,wdt,spr,scl = velocity2(tsr,sol)
            param = np.array([men,sdv,rat,wdt,spr,scl])
            # print param
            # import time
            #
            # time.sleep(10)

        
        rom15 = np.zeros(33)
        rom20 = np.zeros(33)
        rom25 = np.zeros(33)
        rom30 = np.zeros(33)
        rom35 = np.zeros(33)
        rom40 = np.zeros(26)
        rom40f = np.zeros(33)
        
        rom15t = np.zeros(33)
        rom20t = np.zeros(33)
        rom25t = np.zeros(33)
        rom30t = np.zeros(33)
        rom35t = np.zeros(33)
        rom40t = np.zeros(26)
        for i in range(33):
            rom15[i] = velocity_field(xt,yt,0.75,x15[i]*dia,velf,dia,tsr,sol,cfd_data,param)
            rom15t[i] = (rom15[i]-y15[i])/y15[i]
            rom20[i] = velocity_field(xt,yt,1.0,x20[i]*dia,velf,dia,tsr,sol,cfd_data,param)
            rom20t[i] = (rom20[i]-y20[i])/y20[i]
            rom25[i] = velocity_field(xt,yt,1.25,x25[i]*dia,velf,dia,tsr,sol,cfd_data,param)
            rom25t[i] = (rom25[i]-y25[i])/y25[i]
            rom30[i] = velocity_field(xt,yt,1.5,x30[i]*dia,velf,dia,tsr,sol,cfd_data,param)
            rom30t[i] = (rom30[i]-y30[i])/y30[i]
            rom35[i] = velocity_field(xt,yt,1.75,x35[i]*dia,velf,dia,tsr,sol,cfd_data,param)
            rom35t[i] = (rom35[i]-y35[i])/y35[i]
            rom40f[i] = velocity_field(xt,yt,2.0,x35[i]*dia,velf,dia,tsr,sol,cfd_data,param)
            print i
        for i in range(26):
            rom40[i] = velocity_field(xt,yt,4.0*3.,x40[i]*3.,9.308422677,6.,4.5,0.24,cfd_data,param)
            rom40t[i] = (rom40[i]-y40[i])/y40[i]
            print i
        
        rom15error = np.average(fabs(rom15t))
        rom15errorstd = np.std(fabs(rom15t))
        rom20error = np.average(fabs(rom20t))
        rom20errorstd = np.std(fabs(rom20t))
        rom25error = np.average(fabs(rom25t))
        rom25errorstd = np.std(fabs(rom25t))
        rom30error = np.average(fabs(rom30t))
        rom30errorstd = np.std(fabs(rom30t))
        rom35error = np.average(fabs(rom35t))
        rom35errorstd = np.std(fabs(rom35t))
        rom40error = np.average(fabs(rom40t))
        rom40errorstd = np.std(fabs(rom40t))
        oaerror = (rom15error+rom20error+rom25error+rom30error+rom35error+rom40error)/6.
        oaerrorstd = (rom15errorstd+rom20errorstd+rom25errorstd+rom30errorstd+rom35errorstd+rom40errorstd)/6.
        
        fig2 = plt.figure(2,figsize=(12.5,6))
        fig2.subplots_adjust(left=.05,right=.86,wspace=.36,hspace=.35)
        plt.subplot(2,3,1)
        if k == 0:
            plt.plot(x15,y15,'b.')
            plt.plot(pos1,vel1,'r-')
            plt.plot(x15,rom15,'m-')
        elif k == 1:
            plt.plot(x15,rom15,'g-')
        if k == 0:
            plt.xlim(-1,1)
            plt.ylim(0.1,1.2)
            plt.xlabel('$y/D$')
            plt.ylabel(r'$u/U_\infty$')
            plt.text(-0.25,1.05,r'$x/D$ = 0.75')
        print '1.5 mod',(min(rom15)-min(y15))/min(y15),rom15error,rom15errorstd
        print '1.5 cfd',(min(vel1)-min(y15))/min(y15),cfd15error,cfd15errorstd
        plt.subplot(2,3,2)
        if k == 0:
            plt.plot(x20,y20,'b.')
            plt.plot(pos2,vel2,'r-')
            plt.plot(x20,rom20,'m-')
        elif k == 1:
            plt.plot(x20,rom20,'g-')
        if k == 0:
            plt.xlim(-1,1)
            plt.ylim(0.1,1.2)
            plt.xlabel('$y/D$')
            plt.ylabel(r'$u/U_\infty$')
            plt.text(-0.25,1.05,r'$x/D$ = 1.0')
        print '2.0 mod',(min(rom20)-min(y20))/min(y20),rom20error,rom20errorstd
        print '2.0 cfd',(min(vel2)-min(y20))/min(y20),cfd20error,cfd20errorstd
        plt.subplot(2,3,3)
        if k == 0:
            plt.plot(x25,y25,'b.',label='PIV')
            plt.plot(pos3,vel3,'r-',label='CFD')
            plt.plot(x25,rom25,'m-',label='Vorticity')
        elif k == 1:
            plt.plot(x25,rom25,'g-',label='Velocity')
        if k == 0:
            plt.xlim(-1,1)
            plt.ylim(0.1,1.2)
            plt.xlabel('$y/D$')
            plt.ylabel(r'$u/U_\infty$')
            plt.text(-0.25,1.05,r'$x/D$ = 1.25')
        # elif k == 1:
        plt.legend(loc="upper left", bbox_to_anchor=(1,1))
        print '2.5 mod',(min(rom25)-min(y25))/min(y25),rom25error,rom25errorstd
        print '2.5 cfd',(min(vel3)-min(y25))/min(y25),cfd25error,cfd25errorstd
        plt.subplot(2,3,4)
        if k == 0:
            plt.plot(x30,y30,'b.')
            plt.plot(pos4,vel4,'r-')
            plt.plot(x30,rom30,'m-')
        elif k == 1:
            plt.plot(x30,rom30,'g-')
        if k == 0:
            plt.xlim(-1,1)
            plt.ylim(0.1,1.2)
            plt.xlabel('$y/D$')
            plt.ylabel(r'$u/U_\infty$')
            plt.text(-0.25,1.05,r'$x/D$ = 1.5')
        print '3.0 mod',(min(rom30)-min(y30))/min(y30),rom30error,rom30errorstd
        print '3.0 cfd',(min(vel4)-min(y30))/min(y30),cfd30error,cfd30errorstd
        plt.subplot(2,3,5)
        if k == 0:
            plt.plot(x35,y35,'b.')
            plt.plot(pos5,vel5,'r-')
            plt.plot(x35,rom35,'m-')
        elif k == 1:
            plt.plot(x35,rom35,'g-')
        if k == 0:
            plt.xlim(-1,1)
            plt.ylim(0.1,1.2)
            plt.xlabel('$y/D$')
            plt.ylabel(r'$u/U_\infty$')
            plt.text(-0.25,1.05,r'$x/D$ = 1.75')
        print '3.5 mod',(min(rom35)-min(y35))/min(y35),rom35error,rom35errorstd
        print '3.5 cfd',(min(vel5)-min(y35))/min(y35),cfd35error,cfd35errorstd
        plt.subplot(2,3,6)
        if k == 0:
            plt.plot(x40,y40,'b.')
            plt.plot(pos6,vel6,'r-')
            plt.plot(x35,rom40f,'m-')
        elif k == 1:
            plt.plot(x35,rom40f,'g-')
        if k == 0:
            plt.xlim(-1,1)
            plt.ylim(0.1,1.2)
            plt.xlabel('$y/D$')
            plt.ylabel(r'$u/U_\infty$')
            plt.text(-0.25,1.05,r'$x/D$ = 2.0')
        print '4.0 mod',(min(rom40f)-min(y40))/min(y40),rom40error,rom40errorstd
        print '4.0 cfd',(min(vel6)-min(y40))/min(y40),cfd40error,cfd40errorstd

        print oaerror,oaerrorstd
    
plt.show()

# CFD
# 1.5 cfd 0.085222697309 0.126583824064 0.145411541539
# 2.0 cfd 0.100156455326 0.128983354685 0.153469587078
# 2.5 cfd 0.157759607312 0.150857817984 0.170220943027
# 3.0 cfd 0.229735606591 0.185110325841 0.202881236867
# 3.5 cfd 0.221745704619 0.221675365171 0.241148027704
# 4.0 cfd 0.229931173793 0.229559158614 0.230378847745
# 0.173794974393 0.19058503066
# Vort
# 1.5 mod 0.132406448951 0.123814266617 0.0552752953268
# 2.0 mod 0.0899133305592 0.144327259244 0.0820073765436
# 2.5 mod 0.121791964707 0.168153726068 0.1329647045
# 3.0 mod 0.18331069938 0.199964293136 0.162818805975
# 3.5 mod 0.178599395018 0.243762270671 0.204391618207
# 4.0 mod 0.195484403221 0.227921962845 0.158527848792
# 0.18465729643 0.132664274891
# Velo- val
# 1.5 mod 0.0764360005636 0.130154973392 0.107578354084
# 2.0 mod 0.114032956588 0.162965148463 0.11913931666
# 2.5 mod 0.19796690214 0.207943768574 0.159908756532
# 3.0 mod 0.296349051996 0.256526949057 0.213454303636
# 3.5 mod 0.309388725898 0.310579829187 0.2733377453
# 4.0 mod 0.337190466884 0.230182235911 0.132586893453
# 0.216392150764 0.167667561611
# Velo- flip
# 1.5 mod 0.0768091662811 0.0653944355885 0.0726633578795
# 2.0 mod 0.114357543787 0.0787570867999 0.0793591266985
# 2.5 mod 0.198250408813 0.12999380938 0.105127979898
# 3.0 mod 0.296586892386 0.174819197659 0.15220512067
# 3.5 mod 0.309562485502 0.226642726267 0.207679397223
# 4.0 mod 0.337304393801 0.174324433215 0.139728360992
# 0.141655281485 0.126127223893


