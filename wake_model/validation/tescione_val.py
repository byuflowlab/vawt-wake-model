import numpy as np
import matplotlib.pyplot as plt
import csv
from VWM_Fortran import velocity_field
from scipy.io import loadmat

from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

rom = True

r = 0.5
v = 1.

fdata = 'path/to/tes.csv' # adjust this for your case
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


tesdata = loadmat('/path/to/dataH_VAWT.mat') # adjust this for your case


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
    x15[i] = (tesdata['y'][int(i*space),1706])/1000.
    x20[i] = (tesdata['y'][int(i*space),1956])/1000.
    x25[i] = (tesdata['y'][int(i*space),2206])/1000.
    x30[i] = (tesdata['y'][int(i*space),2456])/1000.
    x35[i] = (tesdata['y'][int(i*space),2706])/1000.
    
    y15[i] = (tesdata['u_x'][int(i*space),1706])/9.3
    y20[i] = (tesdata['u_x'][int(i*space),1956])/9.3
    y25[i] = (tesdata['u_x'][int(i*space),2206])/9.3
    y30[i] = (tesdata['u_x'][int(i*space),2456])/9.3
    y35[i] = (tesdata['u_x'][int(i*space),2706])/9.3

for i in range(26):
    x40[i] = (tesdata['y'][int(i*space+7.*space),2956])/1000.
    
    y40[i] = (tesdata['u_x'][int(i*space+7.*space),2956])/9.3
    

plt.figure(1)
plt.subplot(2,3,1)
plt.plot(x15,y15,'b.')
plt.plot(pos1,vel1,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,0.9,r'$x/D$ = 0.75')
print '1.5 cfd',(min(vel1)-min(y15))/min(y15)
plt.subplot(2,3,2)
plt.plot(x20,y20,'b.')
plt.plot(pos2,vel2,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,0.9,r'$x/D$ = 1.0')
print '2.0 cfd',(min(vel2)-min(y20))/min(y20)
plt.subplot(2,3,3)
plt.plot(x25,y25,'b.',label='PIV')
plt.plot(pos3,vel3,'r-',label='CFD')
plt.xlim(-1,1)
plt.ylim(0.1,1)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,0.9,r'$x/D$ = 1.25')
print '2.5 cfd',(min(vel3)-min(y25))/min(y25)
plt.legend(loc="upper left", bbox_to_anchor=(1,1))
plt.subplot(2,3,4)
plt.plot(x30,y30,'b.')
plt.plot(pos4,vel4,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,0.9,r'$x/D$ = 1.5')
print '3.0 cfd',(min(vel4)-min(y30))/min(y30)
plt.subplot(2,3,5)
plt.plot(x35,y35,'b.')
plt.plot(pos5,vel5,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,0.9,r'$x/D$ = 1.75')
print '3.5 cfd',(min(vel5)-min(y35))/min(y35)
plt.subplot(2,3,6)
plt.plot(x40,y40,'b.')
plt.plot(pos6,vel6,'r-')
plt.xlim(-1,1)
plt.ylim(0.1,1)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.25,0.9,r'$x/D$ = 2.0')
print '4.0 cfd',(min(vel6)-min(y40))/min(y40)


if rom == True:
    rad = 0.5
    dia = 2*rad
    velf = 9.308422677
    sol = 0.24
    tsr = 4.5
    rom15 = np.zeros(33)
    rom20 = np.zeros(33)
    rom25 = np.zeros(33)
    rom30 = np.zeros(33)
    rom35 = np.zeros(33)
    rom40 = np.zeros(33)
    for i in range(33):
        rom15[i] = velocity_field(0.75,x15[i]*dia,velf,dia,tsr,sol)
        rom20[i] = velocity_field(1.0,x20[i]*dia,velf,dia,tsr,sol)
        rom25[i] = velocity_field(1.25,x25[i]*dia,velf,dia,tsr,sol)
        rom30[i] = velocity_field(1.5,x30[i]*dia,velf,dia,tsr,sol)
        rom35[i] = velocity_field(1.75,x35[i]*dia,velf,dia,tsr,sol)
        rom40[i] = velocity_field(2.0,x35[i]*dia,velf,dia,tsr,sol)
        print i
    # for i in range(26):
    #     rom40[i] = velocity_field(4.0*3.,x40[i]*3.,9.308422677,6.,4.5,0.24)
    #     print i
    plt.figure(2)
    plt.subplot(2,3,1)
    plt.plot(x15,y15,'b.')
    plt.plot(pos1,vel1,'r-')
    plt.plot(x15,rom15,'g-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1)
    plt.xlabel('$y/D$')
    plt.ylabel(r'$u/U_\infty$')
    plt.text(-0.25,0.9,r'$x/D$ = 0.75')
    print '1.5 mod',(min(rom15)-min(y15))/min(y15)
    print '1.5 cfd',(min(vel1)-min(y15))/min(y15)
    plt.subplot(2,3,2)
    plt.plot(x20,y20,'b.')
    plt.plot(pos2,vel2,'r-')
    plt.plot(x20,rom20,'g-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1)
    plt.xlabel('$y/D$')
    plt.ylabel(r'$u/U_\infty$')
    plt.text(-0.25,0.9,r'$x/D$ = 1.0')
    print '2.0 mod',(min(rom20)-min(y20))/min(y20)
    print '2.0 cfd',(min(vel2)-min(y20))/min(y20)
    plt.subplot(2,3,3)
    plt.plot(x25,y25,'b.',label='PIV')
    plt.plot(pos3,vel3,'r-',label='CFD')
    plt.plot(x25,rom25,'g-',label='Model')
    plt.xlim(-1,1)
    plt.ylim(0.1,1)
    plt.xlabel('$y/D$')
    plt.ylabel(r'$u/U_\infty$')
    plt.text(-0.25,0.9,r'$x/D$ = 1.25')
    print '2.5 mod',(min(rom25)-min(y25))/min(y25)
    print '2.5 cfd',(min(vel3)-min(y25))/min(y25)
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.subplot(2,3,4)
    plt.plot(x30,y30,'b.')
    plt.plot(pos4,vel4,'r-')
    plt.plot(x30,rom30,'g-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1)
    plt.xlabel('$y/D$')
    plt.ylabel(r'$u/U_\infty$')
    plt.text(-0.25,0.9,r'$x/D$ = 1.5')
    print '3.0 mod',(min(rom30)-min(y30))/min(y30)
    print '3.0 cfd',(min(vel4)-min(y30))/min(y30)
    plt.subplot(2,3,5)
    plt.plot(x35,y35,'b.')
    plt.plot(pos5,vel5,'r-')
    plt.plot(x35,rom35,'g-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1)
    plt.xlabel('$y/D$')
    plt.ylabel(r'$u/U_\infty$')
    plt.text(-0.25,0.9,r'$x/D$ = 1.75')
    print '3.5 mod',(min(rom35)-min(y35))/min(y35)
    print '3.5 cfd',(min(vel5)-min(y35))/min(y35)
    plt.subplot(2,3,6)
    plt.plot(x40,y40,'b.')
    plt.plot(pos6,vel6,'r-')
    plt.plot(x35,rom40,'g-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1)
    plt.xlabel('$y/D$')
    plt.ylabel(r'$u/U_\infty$')
    plt.text(-0.25,0.9,r'$x/D$ = 2.0')
    print '4.0 mod',(min(rom40)-min(y40))/min(y40)
    print '4.0 cfd',(min(vel6)-min(y40))/min(y40)
plt.show()


# 1.5 mod 0.145925858405
# 1.5 cfd 0.0959903461365
# 2.0 mod 0.121056813459
# 2.0 cfd 0.129697198407
# 2.5 mod 0.138308060483
# 2.5 cfd 0.173569357736
# 3.0 mod 0.133010208818
# 3.0 cfd 0.17657924573
# 3.5 mod 0.146076519692
# 3.5 cfd 0.187252581349
# 4.0 mod 0.238551828422
# 4.0 cfd 0.273366784078
