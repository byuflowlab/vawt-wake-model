import numpy as np
import matplotlib.pyplot as plt
import csv
from VAWT_Wake_Model import velocity_field
from scipy.io import loadmat
from numpy import fabs
from os import path

from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

rom = True
# rom = False


r = 0.5 # radius
v = 1.0 # velocity

veltype = 'x'
errortype = 'abs'
errortype = 'rel'
epsilon = 1e-3
errortype = 'rms'

load_mat = True # load the original data from dataH_VAWT.mat if available
load_mat = False # use already imported values

# Import Star-CCM+ simulation data (from validation study)
basepath = path.join(path.dirname(path.realpath('__file__')))
fdata = basepath + path.sep + 'tes.csv'

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

if load_mat == True:
    tesdata = loadmat(basepath + path.sep + 'dataH_VAWT.mat')

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

elif load_mat == False:
    x15 = np.array( [-0.867, -0.815, -0.763, -0.71, -0.658, -0.605, -0.553, -0.5, -0.448, -0.395, -0.343, -0.291, -0.238, -0.186, -0.133, -0.081, -0.028, 0.024, 0.077, 0.129, 0.182, 0.234, 0.286, 0.339, 0.391, 0.444, 0.496, 0.549, 0.601, 0.654, 0.706, 0.759, 0.811] )
    x20 = np.array( [-0.867, -0.815, -0.763, -0.71, -0.658, -0.605, -0.553, -0.5, -0.448, -0.395, -0.343, -0.291, -0.238, -0.186, -0.133, -0.081, -0.028, 0.024, 0.077, 0.129, 0.182, 0.234, 0.286, 0.339, 0.391, 0.444, 0.496, 0.549, 0.601, 0.654, 0.706, 0.759, 0.811] )
    x25 = np.array( [-0.867, -0.815, -0.763, -0.71, -0.658, -0.605, -0.553, -0.5, -0.448, -0.395, -0.343, -0.291, -0.238, -0.186, -0.133, -0.081, -0.028, 0.024, 0.077, 0.129, 0.182, 0.234, 0.286, 0.339, 0.391, 0.444, 0.496, 0.549, 0.601, 0.654, 0.706, 0.759, 0.811] )
    x30 = np.array( [-0.867, -0.815, -0.763, -0.71, -0.658, -0.605, -0.553, -0.5, -0.448, -0.395, -0.343, -0.291, -0.238, -0.186, -0.133, -0.081, -0.028, 0.024, 0.077, 0.129, 0.182, 0.234, 0.286, 0.339, 0.391, 0.444, 0.496, 0.549, 0.601, 0.654, 0.706, 0.759, 0.811] )
    x35 = np.array( [-0.867, -0.815, -0.763, -0.71, -0.658, -0.605, -0.553, -0.5, -0.448, -0.395, -0.343, -0.291, -0.238, -0.186, -0.133, -0.081, -0.028, 0.024, 0.077, 0.129, 0.182, 0.234, 0.286, 0.339, 0.391, 0.444, 0.496, 0.549, 0.601, 0.654, 0.706, 0.759, 0.811] )
    x40 = np.array( [-0.5, -0.448, -0.395, -0.343, -0.291, -0.238, -0.186, -0.133, -0.081, -0.028, 0.024, 0.077, 0.129, 0.182, 0.234, 0.286, 0.339, 0.391, 0.444, 0.496, 0.549, 0.601, 0.654, 0.706, 0.759, 0.811] )
    y15 = np.array( [0.994109059733553, 0.9963677091332376, 1.0005562115755295, 1.0129863984422722, 1.0507201982582595, 0.9041429426683261, 0.6716433705834841, 0.5201309340176287, 0.508144821028986, 0.45366206008875143, 0.4120816289764176, 0.3841500969704306, 0.36762989356268116, 0.34878204855733175, 0.3239136841769655, 0.30087738601863206, 0.2879935890797586, 0.26790770257593155, 0.2625797365850182, 0.292431451988392, 0.32394314025359344, 0.3491215866639411, 0.3443858871235991, 0.3482695942305222, 0.3129449223429124, 0.36582747813800515, 0.4154704750687794, 0.49039503797192485, 0.6038500532436828, 0.7927363063069286, 0.9055264980137487, 0.9492042298526092, 0.9678480377419288] )
    y20 = np.array( [0.9920370542708112, 0.9937027721069257, 0.9957821547724749, 1.0042398848773346, 0.9749571783229305, 0.8545201774180555, 0.565802548086458, 0.5143332054399019, 0.489198302575, 0.44408180130876845, 0.3816702901351002, 0.35005844980465, 0.33679345047454295, 0.3230305737612392, 0.3146901353633469, 0.29915244503218225, 0.2790206166599524, 0.2464443364444221, 0.2475139147846546, 0.25357920856345817, 0.27126299099141044, 0.2987673093397647, 0.3182385501649433, 0.3243813328675722, 0.30742297967502097, 0.32253464736566645, 0.3693305499571722, 0.4191276334361715, 0.5015171898966418, 0.6228502433753057, 0.8230607176338183, 0.9264600739810046, 0.9616530515736079] )
    y25 = np.array( [0.9879572983671737, 0.9905509911272896, 0.9933676654604374, 0.9923430478507566, 0.9160865587232668, 0.727774179726256, 0.5833627736796199, 0.4815735966162955, 0.4258818988364446, 0.40697818686203785, 0.3645090556659908, 0.33624432950148214, 0.3254855613810157, 0.3103304514855841, 0.3045151031176352, 0.28401646740896264, 0.2593430020697244, 0.23659060256721978, 0.22300420317944888, 0.2244438460371643, 0.23642978330838485, 0.2568650503421204, 0.28114157843083193, 0.294601202395863, 0.3006303134268269, 0.3118773622351477, 0.3203024532857655, 0.3747931924965308, 0.41075281837916877, 0.5033971635645369, 0.6381178175981282, 0.8832861499445961, 0.9780827152012185] )
    y30 = np.array( [0.9807919873645041, 0.9799943204500705, 0.9756659438025117, 0.9597111733105987, 0.8640379300795783, 0.6756090098761603, 0.5574514345456549, 0.48935854692854497, 0.428523583438216, 0.3833992822748339, 0.3480531699708427, 0.32487761318471114, 0.3153752965437699, 0.2971902045364618, 0.2830498661729626, 0.2701094817124857, 0.2525140339109516, 0.22990689461698513, 0.21388156547631884, 0.2009225725260476, 0.2109170460152375, 0.22197598259760387, 0.24007485599064649, 0.2629099716848817, 0.28306559237296497, 0.29465097651166405, 0.29791611270965696, 0.320843380159032, 0.365988216767817, 0.435994478045424, 0.5459295715799363, 0.7790196684279612, 0.9232640197764616] )
    y35 = np.array( [0.9771079214279068, 0.9739267876288416, 0.9659755072544549, 0.9409586395095753, 0.8372410989082885, 0.6761226132078453, 0.5567979451297009, 0.46883935558555223, 0.41825105219857445, 0.3643627166731387, 0.33608496528240905, 0.31852020274820736, 0.3061917562496233, 0.28608694443477783, 0.273103042317234, 0.26439124635620204, 0.2478225127262987, 0.22900647970420313, 0.20929176925815257, 0.19435346927934738, 0.19534832002282343, 0.19777809452729966, 0.21620168910917673, 0.23300547991382745, 0.25288024387549046, 0.27049877766131897, 0.2804228982454348, 0.2916868458391746, 0.3235522394271216, 0.3723078207006349, 0.4706069281766252, 0.6090953502184843, 0.8713603615558797] )
    y40 = np.array( [0.45465337175396475, 0.3980109024658175, 0.35904507111616846, 0.32774805291916986, 0.31088720860935865, 0.29832097587839296, 0.2821950936769575, 0.26788022954331897, 0.25719920421744913, 0.24586751730014048, 0.2286008521345077, 0.20763836919642048, 0.1912065524478209, 0.1952723101174681, 0.1863567864754898, 0.19527809328409215, 0.21004969265451465, 0.22887599974659367, 0.2507107104057936, 0.26844672767918015, 0.27999051414332765, 0.2936850402917538, 0.33560687204407424, 0.41822453322242903, 0.5404777347489278, 0.7985895165740945] )

    
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
    if errortype == 'abs':
        cfd15t[i] = fabs((1. - vel1[int(index15[i])])-(1. - y15[i]))
    elif errortype == 'rel':
        if fabs(1. - y15[i]) >= epsilon:
            cfd15t[i] = fabs(((1. - vel1[int(index15[i])])-(1. - y15[i]))/(1. - y15[i]))
        else:
            cfd15t[i] = fabs((1. - vel1[int(index15[i])])-(1. - y15[i]))
    elif errortype == 'rms':
        cfd15t[i] = ((1. - vel1[int(index15[i])])-(1. - y15[i]))**2
    if errortype == 'abs':
        cfd20t[i] = fabs((1. - vel2[int(index20[i])])-(1. - y20[i]))
    elif errortype == 'rel':
        if fabs(1. - y20[i]) >= epsilon:
            cfd20t[i] = fabs(((1. - vel2[int(index20[i])])-(1. - y20[i]))/(1. - y20[i]))
        else:
            cfd20t[i] = fabs((1. - vel2[int(index20[i])])-(1. - y20[i]))
    elif errortype == 'rms':
        cfd20t[i] = ((1. - vel2[int(index20[i])])-(1. - y20[i]))**2
    if errortype == 'abs':
        cfd25t[i] = fabs((1. - vel3[int(index25[i])])-(1. - y25[i]))
    elif errortype == 'rel':
        if fabs(1. - y25[i]) >= epsilon:
            cfd25t[i] = fabs(((1. - vel3[int(index25[i])])-(1. - y25[i]))/(1. - y25[i]))
        else:
            cfd25t[i] = fabs((1. - vel3[int(index25[i])])-(1. - y25[i]))
    elif errortype == 'rms':
        cfd25t[i] = ((1. - vel3[int(index25[i])])-(1. - y25[i]))**2
    if errortype == 'abs':
        cfd30t[i] = fabs((1. - vel4[int(index30[i])])-(1. - y30[i]))
    elif errortype == 'rel':
        if fabs(1. - y30[i]) >= epsilon:
            cfd30t[i] = fabs(((1. - vel4[int(index30[i])])-(1. - y30[i]))/(1. - y30[i]))
        else:
            cfd30t[i] = fabs((1. - vel4[int(index30[i])])-(1. - y30[i]))
    elif errortype == 'rms':
        cfd30t[i] = ((1. - vel4[int(index30[i])])-(1. - y30[i]))**2
    if errortype == 'abs':
        cfd35t[i] = fabs((1. - vel5[int(index35[i])])-(1. - y35[i]))
    elif errortype == 'rel':
        if fabs(1. - y35[i]) >= epsilon:
            cfd35t[i] = fabs(((1. - vel5[int(index35[i])])-(1. - y35[i]))/(1. - y35[i]))
        else:
            cfd35t[i] = fabs((1. - vel5[int(index35[i])])-(1. - y35[i]))
    elif errortype == 'rms':
        cfd35t[i] = ((1. - vel5[int(index35[i])])-(1. - y35[i]))**2
for i in range(26):
    if errortype == 'abs':
        cfd40t[i] = fabs((1. - vel6[int(index40[i])])-(1. - y40[i]))
    elif errortype == 'rel':
        if fabs(1. - y40[i]) >= epsilon:
            cfd40t[i] = fabs(((1. - vel6[int(index40[i])])-(1. - y40[i]))/(1. - y40[i]))
        else:
            cfd40t[i] = fabs((1. - vel6[int(index40[i])])-(1. - y40[i]))
    elif errortype == 'rms':
        cfd40t[i] = ((1. - vel6[int(index40[i])])-(1. - y40[i]))**2



if errortype == 'abs' or errortype == 'rel':
    cfd15error = np.average(cfd15t)
    cfd15errorstd = np.std(cfd15t)
    cfd20error = np.average(cfd20t)
    cfd20errorstd = np.std(cfd20t)
    cfd25error = np.average(cfd25t)
    cfd25errorstd = np.std(cfd25t)
    cfd30error = np.average(cfd30t)
    cfd30errorstd = np.std(cfd30t)
    cfd35error = np.average(cfd35t)
    cfd35errorstd = np.std(cfd35t)
    cfd40error = np.average(cfd40t)
    cfd40errorstd = np.std(cfd40t)
elif errortype == 'rms':
    cfd15error = np.sqrt(np.average(cfd15t))
    cfd15errorstd = 1.
    cfd20error = np.sqrt(np.average(cfd20t))
    cfd20errorstd = 1.
    cfd25error = np.sqrt(np.average(cfd25t))
    cfd25errorstd = 1.
    cfd30error = np.sqrt(np.average(cfd30t))
    cfd30errorstd = 1.
    cfd35error = np.sqrt(np.average(cfd35t))
    cfd35errorstd = 1.
    cfd40error = np.sqrt(np.average(cfd40t))
    cfd40errorstd = 1.

cfdoaerror = (cfd15error+cfd20error+cfd25error+cfd30error+cfd35error+cfd40error)/6.
cfdoaerrorstd = (cfd15errorstd+cfd20errorstd+cfd25errorstd+cfd30errorstd+cfd35errorstd+cfd40errorstd)/6.

    
## Plot CFD
fs = 20

fig1 = plt.figure(1,figsize=(12.5,6))
fig1.subplots_adjust(left=.07,bottom=.12,right=.84,wspace=.36,hspace=.5)
plt.subplot(2,3,1)
plt.plot(x15,y15,'g.')
plt.plot(pos1,vel1,'b--',linewidth=2)
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$',fontsize=fs)
plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.text(-0.45,1.05,r'$x/D$ = 0.75',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
print '1.5 cfd',cfd15error,cfd15errorstd
plt.subplot(2,3,2)
plt.plot(x20,y20,'g.')
plt.plot(pos2,vel2,'b--',linewidth=2)
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$',fontsize=fs)
plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.text(-0.45,1.05,r'$x/D$ = 1.0',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
print '2.0 cfd',cfd20error,cfd20errorstd
plt.subplot(2,3,3)
plt.plot(x25,y25,'g.',label='PIV')
plt.plot(pos3,vel3,'b--',linewidth=2,label='CFD')
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$',fontsize=fs)
plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.text(-0.45,1.05,r'$x/D$ = 1.25',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
print '2.5 cfd',cfd25error,cfd25errorstd
plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
plt.subplot(2,3,4)
plt.plot(x30,y30,'g.')
plt.plot(pos4,vel4,'b--',linewidth=2)
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$',fontsize=fs)
plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.text(-0.45,1.05,r'$x/D$ = 1.5',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
print '3.0 cfd',cfd30error,cfd30errorstd
plt.subplot(2,3,5)
plt.plot(x35,y35,'g.')
plt.plot(pos5,vel5,'b--',linewidth=2)
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$',fontsize=fs)
plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.text(-0.45,1.05,r'$x/D$ = 1.75',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
print '3.5 cfd',cfd35error,cfd35errorstd
plt.subplot(2,3,6)
plt.plot(x40,y40,'g.')
plt.plot(pos6,vel6,'b--',linewidth=2)
plt.xlim(-1,1)
plt.ylim(0.1,1.2)
plt.xlabel('$y/D$',fontsize=fs)
plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.text(-0.45,1.05,r'$x/D$ = 2.0',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
print '4.0 cfd',cfd40error,cfd40errorstd


print '\nOverall Error:',cfdoaerror
print 'Overall Stand Dev:',cfdoaerrorstd,'\n'



## Plot Model
if rom == True:
    rad = 0.5
    dia = 2*rad
    velf = 9.308422677
    sol = 0.24
    tsr = 4.5
    rot = tsr*velf/rad
    chord = 0.06
    B = 2
    xt = 0.
    yt = 0.

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
        rom15[i] = velocity_field(xt,yt,0.75*dia,x15[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
        rom20[i] = velocity_field(xt,yt,1.0*dia,x20[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
        rom25[i] = velocity_field(xt,yt,1.25*dia,x25[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
        rom30[i] = velocity_field(xt,yt,1.5*dia,x30[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
        rom35[i] = velocity_field(xt,yt,1.75*dia,x35[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
        rom40f[i] = velocity_field(xt,yt,2.0*dia,x35[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
        print i
    for i in range(26):
        rom40[i] = velocity_field(xt,yt,2.0*dia,x40[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
        print i

    for i in range(33):
        if errortype == 'abs':
            rom15t[i] = fabs((1. - rom15[i])-(1. - y15[i]))
        elif errortype == 'rel':
            if fabs(1. - y15[i]) >= epsilon:
                rom15t[i] = fabs(((1. - rom15[i])-(1. - y15[i]))/(1. - y15[i]))
            else:
                rom15t[i] = fabs((1. - rom15[i])-(1. - y15[i]))
        elif errortype == 'rms':
            rom15t[i] = ((1. - rom15[i])-(1. - y15[i]))**2.
        if errortype == 'abs':
            rom20t[i] = fabs((1. - rom20[i])-(1. - y20[i]))
        elif errortype == 'rel':
            if fabs(1. - y20[i]) >= epsilon:
                rom20t[i] = fabs(((1. - rom20[i])-(1. - y20[i]))/(1. - y20[i]))
            else:
                rom20t[i] = fabs((1. - rom20[i])-(1. - y20[i]))
        elif errortype == 'rms':
            rom20t[i] = ((1. - rom20[i])-(1. - y20[i]))**2.
        if errortype == 'abs':
            rom25t[i] = fabs((1. - rom25[i])-(1. - y25[i]))
        elif errortype == 'rel':
            if fabs(1. - y25[i]) >= epsilon:
                rom25t[i] = fabs(((1. - rom25[i])-(1. - y25[i]))/(1. - y25[i]))
            else:
                rom25t[i] = fabs((1. - rom25[i])-(1. - y25[i]))
        elif errortype == 'rms':
            rom25t[i] = ((1. - rom25[i])-(1. - y25[i]))**2.
        if errortype == 'abs':
            rom30t[i] = fabs((1. - rom30[i])-(1. - y30[i]))
        elif errortype == 'rel':
            if fabs(1. - y30[i]) >= epsilon:
                rom30t[i] = fabs(((1. - rom30[i])-(1. - y30[i]))/(1. - y30[i]))
            else:
                rom30t[i] = fabs((1. - rom30[i])-(1. - y30[i]))
        elif errortype == 'rms':
            rom30t[i] = ((1. - rom30[i])-(1. - y30[i]))**2.
        if errortype == 'abs':
            rom35t[i] = fabs((1. - rom35[i])-(1. - y35[i]))
        elif errortype == 'rel':
            if fabs(1. - y35[i]) >= epsilon:
                rom35t[i] = fabs(((1. - rom35[i])-(1. - y35[i]))/(1. - y35[i]))
            else:
                rom35t[i] = fabs((1. - rom35[i])-(1. - y35[i]))
        elif errortype == 'rms':
            rom35t[i] = ((1. - rom35[i])-(1. - y35[i]))**2.
    for i in range(26):
        if errortype == 'abs':
            rom40t[i] = fabs((1. - rom40[i])-(1. - y40[i]))
        elif errortype == 'rel':
            if fabs(1. - y40[i]) >= epsilon:
                rom40t[i] = fabs(((1. - rom40[i])-(1. - y40[i]))/(1. - y40[i]))
            else:
                rom40t[i] = fabs((1. - rom40[i])-(1. - y40[i]))
        elif errortype == 'rms':
            rom40t[i] = ((1. - rom40[i])-(1. - y40[i]))**2.


    if errortype == 'abs' or errortype == 'rel':
        rom15error = np.average(rom15t)
        rom15errorstd = np.std(rom15t)
        rom20error = np.average(rom20t)
        rom20errorstd = np.std(rom20t)
        rom25error = np.average(rom25t)
        rom25errorstd = np.std(rom25t)
        rom30error = np.average(rom30t)
        rom30errorstd = np.std(rom30t)
        rom35error = np.average(rom35t)
        rom35errorstd = np.std(rom35t)
        rom40error = np.average(rom40t)
        rom40errorstd = np.std(rom40t)
    elif errortype == 'rms':
        rom15error = np.sqrt(np.average(rom15t))
        rom15errorstd = 1.
        rom20error = np.sqrt(np.average(rom20t))
        rom20errorstd = 1.
        rom25error = np.sqrt(np.average(rom25t))
        rom25errorstd = 1.
        rom30error = np.sqrt(np.average(rom30t))
        rom30errorstd = 1.
        rom35error = np.sqrt(np.average(rom35t))
        rom35errorstd = 1.
        rom40error = np.sqrt(np.average(rom40t))
        rom40errorstd = 1.
    oaerror = (rom15error+rom20error+rom25error+rom30error+rom35error+rom40error)/6.
    oaerrorstd = (rom15errorstd+rom20errorstd+rom25errorstd+rom30errorstd+rom35errorstd+rom40errorstd)/6.

    rom15f = velocity_field(xt,yt,0.75*dia,-1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom20f = velocity_field(xt,yt,1.0*dia,-1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom25f = velocity_field(xt,yt,1.25*dia,-1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom30f = velocity_field(xt,yt,1.5*dia,-1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom35f = velocity_field(xt,yt,1.75*dia,-1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom40ff = velocity_field(xt,yt,2.0*dia,-1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom15l = velocity_field(xt,yt,0.75*dia,1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom20l = velocity_field(xt,yt,1.0*dia,1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom25l = velocity_field(xt,yt,1.25*dia,1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom30l = velocity_field(xt,yt,1.5*dia,1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom35l = velocity_field(xt,yt,1.75*dia,1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    rom40fl = velocity_field(xt,yt,2.0*dia,1.*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)

    rom15 = np.insert(rom15,0,rom15f)
    rom20 = np.insert(rom20,0,rom20f)
    rom25 = np.insert(rom25,0,rom25f)
    rom30 = np.insert(rom30,0,rom30f)
    rom35 = np.insert(rom35,0,rom35f)
    rom40f = np.insert(rom40f,0,rom40ff)
    rom15 = np.append(rom15,rom15l)
    rom20 = np.append(rom20,rom20l)
    rom25 = np.append(rom25,rom25l)
    rom30 = np.append(rom30,rom30l)
    rom35 = np.append(rom35,rom35l)
    rom40f = np.append(rom40f,rom40fl)

    x15p = np.copy(x15)
    x20p = np.copy(x20)
    x25p = np.copy(x25)
    x30p = np.copy(x30)
    x35p = np.copy(x35)
    x40p = np.copy(x35)

    x15p = np.insert(x15p,0,-1.)
    x20p = np.insert(x20p,0,-1.)
    x25p = np.insert(x25p,0,-1.)
    x30p = np.insert(x30p,0,-1.)
    x35p = np.insert(x35p,0,-1.)
    x40p = np.insert(x40p,0,-1.)
    x15p = np.append(x15p,1.)
    x20p = np.append(x20p,1.)
    x25p = np.append(x25p,1.)
    x30p = np.append(x30p,1.)
    x35p = np.append(x35p,1.)
    x40p = np.append(x40p,1.)



    fig2 = plt.figure(2,figsize=(12.5,6))
    fig2.subplots_adjust(left=.07,bottom=.12,right=.84,wspace=.36,hspace=.50)
    plt.subplot(2,3,1)
    plt.plot(x15,y15,'g.')
    plt.plot(pos1,vel1,'b--',linewidth=2)
    plt.plot(x15p,rom15,'r-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1.2)
    plt.xlabel('$y/D$',fontsize=fs)
    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
    plt.text(-0.45,1.05,r'$x/D$ = 0.75',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    print '1.5 mod',(min(rom15)-min(y15))/min(y15),rom15error,rom15errorstd
    print '1.5 cfd',(min(vel1)-min(y15))/min(y15),cfd15error,cfd15errorstd

    plt.subplot(2,3,2)
    plt.plot(x20,y20,'g.')
    plt.plot(pos2,vel2,'b--',linewidth=2)
    plt.plot(x20p,rom20,'r-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1.2)
    plt.xlabel('$y/D$',fontsize=fs)
    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
    plt.text(-0.45,1.05,r'$x/D$ = 1.0',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    print '2.0 mod',(min(rom20)-min(y20))/min(y20),rom20error,rom20errorstd
    print '2.0 cfd',(min(vel2)-min(y20))/min(y20),cfd20error,cfd20errorstd

    plt.subplot(2,3,3)
    plt.plot(x25,y25,'g.',label='PIV')
    plt.plot(pos3,vel3,'b--',linewidth=2,label='CFD')
    plt.plot(x25p,rom25,'r-',label='Model')
    plt.xlim(-1,1)
    plt.ylim(0.1,1.2)
    plt.xlabel('$y/D$',fontsize=fs)
    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
    plt.text(-0.45,1.05,r'$x/D$ = 1.25',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    # elif k == 1:
    plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
    print '2.5 mod',(min(rom25)-min(y25))/min(y25),rom25error,rom25errorstd
    print '2.5 cfd',(min(vel3)-min(y25))/min(y25),cfd25error,cfd25errorstd

    plt.subplot(2,3,4)
    plt.plot(x30,y30,'g.')
    plt.plot(pos4,vel4,'b--',linewidth=2)
    plt.plot(x30p,rom30,'r-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1.2)
    plt.xlabel('$y/D$',fontsize=fs)
    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
    plt.text(-0.45,1.05,r'$x/D$ = 1.5',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    print '3.0 mod',(min(rom30)-min(y30))/min(y30),rom30error,rom30errorstd
    print '3.0 cfd',(min(vel4)-min(y30))/min(y30),cfd30error,cfd30errorstd
    plt.subplot(2,3,5)
    plt.plot(x35,y35,'g.')
    plt.plot(pos5,vel5,'b--',linewidth=2)
    plt.plot(x35p,rom35,'r-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1.2)
    plt.xlabel('$y/D$',fontsize=fs)
    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
    plt.text(-0.45,1.05,r'$x/D$ = 1.75',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    print '3.5 mod',(min(rom35)-min(y35))/min(y35),rom35error,rom35errorstd
    print '3.5 cfd',(min(vel5)-min(y35))/min(y35),cfd35error,cfd35errorstd
    plt.subplot(2,3,6)
    plt.plot(x40,y40,'g.')
    plt.plot(pos6,vel6,'b--',linewidth=2)
    plt.plot(x40p,rom40f,'r-')
    plt.xlim(-1,1)
    plt.ylim(0.1,1.2)
    plt.xlabel('$y/D$',fontsize=fs)
    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
    plt.text(-0.45,1.05,r'$x/D$ = 2.0',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    print '4.0 mod',(min(rom40f)-min(y40))/min(y40),rom40error,rom40errorstd
    print '4.0 cfd',(min(vel6)-min(y40))/min(y40),cfd40error,cfd40errorstd


    if errortype == 'abs':
        print '\n Average Absolute Error'
    elif errortype == 'rel':
        print '\n Average Relative Error'
    elif errortype == 'rms':
        print '\n Root Mean Squared Error'
    print '\nAverage CFD Error:',cfdoaerror
    print 'Average CFD Stand Dev:',cfdoaerrorstd
    print '\nAverage ROM Error:',oaerror
    print 'Average ROM Stand Dev:',oaerrorstd

    low15 = y15.min()
    low20 = y20.min()
    low25 = y25.min()
    low30 = y30.min()
    low35 = y35.min()
    low40 = y40.min()
    low15cfd = vel1.min()
    low20cfd = vel2.min()
    low25cfd = vel3.min()
    low30cfd = vel4.min()
    low35cfd = vel5.min()
    low40cfd = vel6.min()
    low15rom = rom15.min()
    low20rom = rom20.min()
    low25rom = rom25.min()
    low30rom = rom30.min()
    low35rom = rom35.min()
    low40rom = rom40.min()

    if errortype == 'abs' or errortype == 'rms':
        mdfcfd15 = fabs((1.-low15cfd)-(1.-low15))
        mdfcfd20 = fabs((1.-low20cfd)-(1.-low20))
        mdfcfd25 = fabs((1.-low25cfd)-(1.-low25))
        mdfcfd30 = fabs((1.-low30cfd)-(1.-low30))
        mdfcfd35 = fabs((1.-low35cfd)-(1.-low35))
        mdfcfd40 = fabs((1.-low40cfd)-(1.-low40))

        mdfrom15 = fabs((1.-low15rom)-(1.-low15))
        mdfrom20 = fabs((1.-low20rom)-(1.-low20))
        mdfrom25 = fabs((1.-low25rom)-(1.-low25))
        mdfrom30 = fabs((1.-low30rom)-(1.-low30))
        mdfrom35 = fabs((1.-low35rom)-(1.-low35))
        mdfrom40 = fabs((1.-low40rom)-(1.-low40))

    elif errortype == 'rel':
        mdfcfd15 = fabs(((1.-low15cfd)-(1.-low15))/(1.-low15))
        mdfcfd20 = fabs(((1.-low20cfd)-(1.-low20))/(1.-low20))
        mdfcfd25 = fabs(((1.-low25cfd)-(1.-low25))/(1.-low25))
        mdfcfd30 = fabs(((1.-low30cfd)-(1.-low30))/(1.-low30))
        mdfcfd35 = fabs(((1.-low35cfd)-(1.-low35))/(1.-low35))
        mdfcfd40 = fabs(((1.-low40cfd)-(1.-low40))/(1.-low40))

        mdfrom15 = fabs(((1.-low15rom)-(1.-low15))/(1.-low15))
        mdfrom20 = fabs(((1.-low20rom)-(1.-low20))/(1.-low20))
        mdfrom25 = fabs(((1.-low25rom)-(1.-low25))/(1.-low25))
        mdfrom30 = fabs(((1.-low30rom)-(1.-low30))/(1.-low30))
        mdfrom35 = fabs(((1.-low35rom)-(1.-low35))/(1.-low35))
        mdfrom40 = fabs(((1.-low40rom)-(1.-low40))/(1.-low40))

    print 'Maximum Deficit Error (CFD):',mdfcfd15,mdfcfd20,mdfcfd25,mdfcfd30,mdfcfd35,mdfcfd40,'\n\t',max(mdfcfd15,mdfcfd20,mdfcfd25,mdfcfd30,mdfcfd35,mdfcfd40)
    print 'Maximum Deficit Error (ROM):',mdfrom15,mdfrom20,mdfrom25,mdfrom30,mdfrom35,mdfrom40,'\n\t',max(mdfrom15,mdfrom20,mdfrom25,mdfrom30,mdfrom35,mdfrom40)

    
plt.show()

