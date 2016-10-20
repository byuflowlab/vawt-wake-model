import csv
import numpy as np
import matplotlib.pyplot as plt
from VAWT_Wake_Model import velocity_field

pos = np.array([])
vel = np.array([])


fdata = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/y_vel.csv'
f = open(fdata)

csv_f = csv.reader(f)

i = 0
for row in csv_f:
    if i != 0:
        pos = np.append(pos,float(row[0]))
        vel = np.append(vel,float(row[1]))
    i += 1


f.close()

delval = 3.3
delvec = np.array([])
for j in range(np.size(pos)):
    if pos[j] > -delval and pos[j] < delval:
        delvec = np.append(delvec,j)
delvec = delvec[::-1]
for j in range(np.size(delvec)):
    vel = np.delete(vel,delvec[j])
    pos = np.delete(pos,delvec[j])


print 'data'


posy = np.linspace(-50,-delval,100)
posy2 = np.linspace(delval,50,100)
posy = np.append(posy,posy2)
vely = np.zeros_like(posy)

for j in range(np.size(posy)):
    vely[j] = velocity_field(0., 0., 0., posy[j], 15., 6., 4.0, 0.25, veltype='y')*15.
    print j

plt.figure()
plt.plot(vel,pos,'b.',label='CFD')
plt.plot(vely,posy,'r.',label='Model')
plt.xlabel('Induced Y-Velocity (m/s)')
plt.ylabel('Lateral Position (m)')
plt.ylim(-50,50)
plt.legend(loc=1)
plt.show()
