import csv
from os import path

import numpy as np
from sklearn.cross_validation import train_test_split
from scipy.optimize import curve_fit
from sys import argv

import _vawtwake

def starccm_read_vel(fdata,dia,windd,length):
    start = length/30.
    lendat = np.linspace(start,length,30)/dia

    for j in range(np.size(fdata)):
        for k in range(30):
            name = str(k+1)
            exec('pos'+name+' = np.array([])')
            exec('vel'+name+' = np.array([])')
        wind = windd[j]

        f = open(fdata[j])

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
            if row[13] != 'null' and i != 0:
                vel7 = np.append(vel7,float(row[13]))
            if row[15] != 'null' and i != 0:
                vel8 = np.append(vel8,float(row[15]))
            if row[17] != 'null' and i != 0:
                vel9 = np.append(vel9,float(row[17]))
            if row[19] != 'null' and i != 0:
                vel10 = np.append(vel10,float(row[19]))
            if row[21] != 'null' and i != 0:
                vel11 = np.append(vel11,float(row[21]))
            if row[23] != 'null' and i != 0:
                vel12 = np.append(vel12,float(row[23]))
            if row[25] != 'null' and i != 0:
                vel13 = np.append(vel13,float(row[25]))
            if row[27] != 'null' and i != 0:
                vel14 = np.append(vel14,float(row[27]))
            if row[29] != 'null' and i != 0:
                vel15 = np.append(vel15,float(row[29]))
            if row[31] != 'null' and i != 0:
                vel16 = np.append(vel16,float(row[31]))
            if row[33] != 'null' and i != 0:
                vel17 = np.append(vel17,float(row[33]))
            if row[35] != 'null' and i != 0:
                vel18 = np.append(vel18,float(row[35]))
            if row[37] != 'null' and i != 0:
                vel19 = np.append(vel19,float(row[37]))
            if row[39] != 'null' and i != 0:
                vel20 = np.append(vel20,float(row[39]))
            if row[41] != 'null' and i != 0:
                vel21 = np.append(vel21,float(row[41]))
            if row[43] != 'null' and i != 0:
                vel22 = np.append(vel22,float(row[43]))
            if row[45] != 'null' and i != 0:
                vel23 = np.append(vel23,float(row[45]))
            if row[47] != 'null' and i != 0:
                vel24 = np.append(vel24,float(row[47]))
            if row[49] != 'null' and i != 0:
                vel25 = np.append(vel25,float(row[49]))
            if row[51] != 'null' and i != 0:
                vel26 = np.append(vel26,float(row[51]))
            if row[53] != 'null' and i != 0:
                vel27 = np.append(vel27,float(row[53]))
            if row[55] != 'null' and i != 0:
                vel28 = np.append(vel28,float(row[55]))
            if row[57] != 'null' and i != 0:
                vel29 = np.append(vel29,float(row[57]))
            if row[59] != 'null' and i != 0:
                vel30 = np.append(vel30,float(row[59]))
            i += 1

        f.close()

        posdn = np.array([])
        poslt = np.array([])
        veld = np.array([])
        for i in range(30):
            name = str(i+1)
            ind = str(i)

            #Ordering the data numerically by the position
            exec('pos'+name+', vel'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', vel'+name+'))))')

            #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
            exec('pos'+name+'_0 = np.array([])\nvel'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvel'+name+'_0 = np.append(vel'+name+'_0,vel'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)/dia\nvel'+name+' = np.copy(vel'+name+'_0)/wind')

            #Deleting wall boundary data
            exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] > 5. or pos'+name+'[j] < -5.:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvel'+name+' = np.delete(vel'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')

            # #Deleting values greater than free stream wind speed
            # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif vel'+name+'[j] > 1. or fabs(pos'+name+'[j]) > 1.2:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvel'+name+' = np.delete(vel'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')

            # Creating same-sized arrays for downstream position, lateral position, and velocity
            exec('lensize = np.size(pos'+name+')')
            exec('posdn = np.append(posdn,np.ones(lensize)*lendat['+ind+'])\nposlt = np.append(poslt,pos'+name+')\nveld = np.append(veld,vel'+name+')')

        #Deleting values inside turbine region (and shortly after)
        delvec = np.array([])
        for j in range(np.size(posdn)):
            if posdn[j] < 0.58 or poslt[j] > 1.08 or poslt[j] < -1.08:
                delvec = np.append(delvec,j)
        delvec = delvec[::-1]
        for k in range(np.size(delvec)):
            posdn = np.delete(posdn,delvec[k])
            poslt = np.delete(poslt,delvec[k])
            veld = np.delete(veld,delvec[k])

    return posdn,poslt,veld


def starccm_read_vort(fdata,dia,rotd,length):
    start = length/30.
    lendat = np.linspace(start,length,30)/dia

    for j in range(np.size(fdata)):
        for k in range(30):
            name = str(k+1)
            exec('pos'+name+' = np.array([])')
            exec('vort'+name+' = np.array([])')
        rot = rotd[j]

        f = open(fdata[j])

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
                vort1 = np.append(vort1,float(row[1]))
            if row[3] != 'null' and i != 0:
                vort2 = np.append(vort2,float(row[3]))
            if row[5] != 'null' and i != 0:
                vort3 = np.append(vort3,float(row[5]))
            if row[7] != 'null' and i != 0:
                vort4 = np.append(vort4,float(row[7]))
            if row[9] != 'null' and i != 0:
                vort5 = np.append(vort5,float(row[9]))
            if row[11] != 'null' and i != 0:
                vort6 = np.append(vort6,float(row[11]))
            if row[13] != 'null' and i != 0:
                vort7 = np.append(vort7,float(row[13]))
            if row[15] != 'null' and i != 0:
                vort8 = np.append(vort8,float(row[15]))
            if row[17] != 'null' and i != 0:
                vort9 = np.append(vort9,float(row[17]))
            if row[19] != 'null' and i != 0:
                vort10 = np.append(vort10,float(row[19]))
            if row[21] != 'null' and i != 0:
                vort11 = np.append(vort11,float(row[21]))
            if row[23] != 'null' and i != 0:
                vort12 = np.append(vort12,float(row[23]))
            if row[25] != 'null' and i != 0:
                vort13 = np.append(vort13,float(row[25]))
            if row[27] != 'null' and i != 0:
                vort14 = np.append(vort14,float(row[27]))
            if row[29] != 'null' and i != 0:
                vort15 = np.append(vort15,float(row[29]))
            if row[31] != 'null' and i != 0:
                vort16 = np.append(vort16,float(row[31]))
            if row[33] != 'null' and i != 0:
                vort17 = np.append(vort17,float(row[33]))
            if row[35] != 'null' and i != 0:
                vort18 = np.append(vort18,float(row[35]))
            if row[37] != 'null' and i != 0:
                vort19 = np.append(vort19,float(row[37]))
            if row[39] != 'null' and i != 0:
                vort20 = np.append(vort20,float(row[39]))
            if row[41] != 'null' and i != 0:
                vort21 = np.append(vort21,float(row[41]))
            if row[43] != 'null' and i != 0:
                vort22 = np.append(vort22,float(row[43]))
            if row[45] != 'null' and i != 0:
                vort23 = np.append(vort23,float(row[45]))
            if row[47] != 'null' and i != 0:
                vort24 = np.append(vort24,float(row[47]))
            if row[49] != 'null' and i != 0:
                vort25 = np.append(vort25,float(row[49]))
            if row[51] != 'null' and i != 0:
                vort26 = np.append(vort26,float(row[51]))
            if row[53] != 'null' and i != 0:
                vort27 = np.append(vort27,float(row[53]))
            if row[55] != 'null' and i != 0:
                vort28 = np.append(vort28,float(row[55]))
            if row[57] != 'null' and i != 0:
                vort29 = np.append(vort29,float(row[57]))
            if row[59] != 'null' and i != 0:
                vort30 = np.append(vort30,float(row[59]))
            i += 1

        f.close()

        posdn = np.array([])
        poslt = np.array([])
        vortd = np.array([])
        for i in range(30):
            name = str(i+1)
            ind = str(i)

            # Ordering the data numerically by the position
            exec('pos'+name+', vort'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', vort'+name+'))))')

            # STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
            exec('pos'+name+'_0 = np.array([])\nvort'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvort'+name+'_0 = np.append(vort'+name+'_0,vort'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)/dia\nvort'+name+' = np.copy(vort'+name+'_0)/rot')

            # Creating same-sized arrays for downstream position, lateral position, and vorticity strength
            exec('lensize = np.size(pos'+name+')')
            exec('posdn = np.append(posdn,np.ones(lensize)*lendat['+ind+'])\nposlt = np.append(poslt,pos'+name+')\nvortd = np.append(vortd,vort'+name+')')

        # Deleting values inside turbine region (and shortly after)
        delvec = np.array([])
        for j in range(np.size(posdn)):
            if posdn[j] < 0.58 or poslt[j] > 1.08 or poslt[j] < -1.08:
                delvec = np.append(delvec,j)
        delvec = delvec[::-1]
        for k in range(np.size(delvec)):
            posdn = np.delete(posdn,delvec[k])
            poslt = np.delete(poslt,delvec[k])
            vortd = np.delete(vortd,delvec[k])

    return posdn,poslt,vortd


def externaldata_read_vel(fdata,dia,windd,lengthpoints):
    ndata = np.size(lengthpoints)
    lendat = lengthpoints/dia
    nodat = 'null'

    for j in range(np.size(fdata)):
        for k in range(ndata):
            name = str(k+1)
            exec('pos'+name+' = np.array([])')
            exec('vel'+name+' = np.array([])')
        wind = windd[j]

        f = open(fdata[j])

        csv_f = csv.reader(f)

        i = 0
        for row in csv_f:
            for p in range(ndata):
                name = str(p+1)
                exec('if row[p*2] != nodat and i != 0:\n\tpos'+name+' = np.append(pos'+name+',float(row[p*2]))')
                exec('if row[p*2+1] != nodat and i != 0:\n\tvel'+name+' = np.append(vel'+name+',float(row[p*2+1]))')
            i += 1

        f.close()

        posdn = np.array([])
        poslt = np.array([])
        veld = np.array([])
        for i in range(ndata):
            name = str(i+1)
            ind = str(i)

            #Ordering the data numerically by the position
            exec('pos'+name+', vel'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', vel'+name+'))))')

            #Eliminating repeated data
            exec('pos'+name+'_0 = np.array([])\nvel'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvel'+name+'_0 = np.append(vel'+name+'_0,vel'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)/dia\nvel'+name+' = np.copy(vel'+name+'_0)/wind')

            # Creating same-sized arrays for downstream position, lateral position, and velocity
            exec('lensize = np.size(pos'+name+')')
            exec('posdn = np.append(posdn,np.ones(lensize)*lendat['+ind+'])\nposlt = np.append(poslt,pos'+name+')\nveld = np.append(veld,vel'+name+')')

        #Deleting values inside turbine region (and shortly after)
        delvec = np.array([])
        for j in range(np.size(posdn)):
            if posdn[j] < 0.58 or poslt[j] > 1.08 or poslt[j] < -1.08:
                delvec = np.append(delvec,j)
        delvec = delvec[::-1]
        for k in range(np.size(delvec)):
            posdn = np.delete(posdn,delvec[k])
            poslt = np.delete(poslt,delvec[k])
            veld = np.delete(veld,delvec[k])

    return posdn,poslt,veld



def sheet(pos,param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18,param19,param20,param21,param22,param23,param24,param25,param26,param27,param28,param29,param30,param31,param32,param33,param34,param35,param36,param37,param38,param39,param40,param41,param42,param43,param44,param45,param46,param47,param48,param49,param50,param51,param52,param53,param54,param55,param56,param57,param58,param59,param60,param61,param62,param63,param64,param65,param66,param67,param68,param69,param70,param71,param72,param73,param74,param75,param76,param77,param78,param79,param80,param81,param82,param83,param84,param85,param86,param87,param88,param89,param90,param91,param92,param93,param94,param95,param96,param97,param98,param99,param100):
    global xttr
    global ystr
    global fit_type

    coef0 = np.array([])
    coef1 = np.array([])
    coef2 = np.array([])
    coef3 = np.array([])
    coef4 = np.array([])
    coef5 = np.array([])
    coef6 = np.array([])
    coef7 = np.array([])
    coef8 = np.array([])
    coef9 = np.array([])

    k = 1
    for i in range(10):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',param'+str(k)+')')
            k += 1

    posdn = pos[0,:]
    poslt = pos[1,:]

    # Ensuring the correct order of fit for each coefficient
    coef0[9] = 0. # loc1; TSR: 3; Solidity: 2
    coef3[6] = 0. # spr1; TSR: 2; Solidity: 3
    coef4[9] = 0. # spr2; TSR: 3; Solidity: 2
    coef5[6] = 0. # skw1; TSR: 2; Solidity: 3
    coef7[6] = 0. # scl1; TSR: 2; Solidity: 3

    if fit_type == 'vel':
        model = _vawtwake.sheet_vel(xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.,1.,220,200,1)
    elif fit_type == 'vort':
        model = _vawtwake.sheet_vort(xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.)

    return model

# Main File
if __name__ == "__main__":

    global xttr
    global ystr
    global posdntr_opt
    global vortdtr_opt
    global fit_type

    # Settings for purpose
    # cv_on = 0 # cross-validation
    cv_on = 1 # retraining with all data

    # STAR-CCM+ data sets to read in
    read_data = 4

    fit_type = 'vel'
    fit_type = 'vort'

    add_data = True
    add_data = False

    if cv_on == 0:
        print '********Cross Validating Wake Model********'
    elif cv_on == 1:
        print '********Retraining Wake Model********'
    if fit_type == 'vel':
        print '\nReading in '+str(read_data)+' STAR-CCM+ velocity data sets'
    elif fit_type == 'vort':
        print '\nReading in '+str(read_data)+' STAR-CCM+ vorticity data sets'


    # STAR-CCM+ wake data lengths
    s1length = np.array([210.,210.,205.,196.,185.,178.,170.,165.,160.,145.,140.,123.,115.,112.,108.,101.,101.,90.,85.,80.,78.,75.,70.])
    s2length = np.array([197.,193.,185.,176.,140.,146.,126.,114.,103.,96.,100.,86.,77.,72.,70.,68.,60.,64.,54.,50.,47.,45.,44.])
    s3length = np.array([185.,150.,100.,95.,83.,76.,72.,63.,60.,49.,50.,41.,39.,36.,34.,33.,30.,31.,28.,30.,29.,28.,27.])
    s4length = np.array([145.,100.,73.,60.,53.,44.,42.,37.,38.,30.,33.,26.,22.,24.,23.,21.,21.,19.,24.,23.,22.,21.,20.])
    s5length = np.array([78.,70.,52.,43.,37.,32.,29.,27.,26.,23.,20.,20.,23.,21.,20.,19.,19.,18.,18.,16.,16.,15.,14.])

    tsr = np.linspace(150,700,23)
    solidity = np.array(['0.15','0.25','0.50','0.75','1.0'])

    tsr_cv = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])
    solidity_cv = np.array([1,2,3,4,5])

    tsrr_cv = tsr/100.
    sol_cv = np.array([0.15,0.25,0.5,0.75,1.0])

    s = np.array([])
    t = np.array([])
    scv = np.array([])
    tcv = np.array([])
    xtcv = np.array([])
    yscv = np.array([])
    for i in range(np.size(solidity)):
        for j in range(np.size(tsr)):
            s = np.append(s,solidity[i])
            t = np.append(t,tsr[j])
            scv = np.append(scv,solidity_cv[i])
            tcv = np.append(tcv,tsr_cv[j])
            xtcv = np.append(xtcv,tsrr_cv[j])
            yscv = np.append(yscv,sol_cv[i])

    slength = np.array([])
    slength = np.append(slength,s1length)
    slength = np.append(slength,s2length)
    slength = np.append(slength,s3length)
    slength = np.append(slength,s4length)
    slength = np.append(slength,s5length)
    slength = slength

    q = 0
    print '\nSTAR-CCM+ data import complete at:'
    for i in range(5):
        for j in range(23):
            sname = str(i+1)
            tname = str(j+1)

            wind1 = 15.
            wind2 = 14.
            wind3 = 12.
            wind4 = 16.

            rad = 3.
            dia = rad*2.
            tsrf = t[q]/100.
            rot1 = tsrf*wind1/rad
            rot2 = tsrf*wind2/rad
            rot3 = tsrf*wind3/rad
            rot4 = tsrf*wind4/rad

            wfit1 = 'w'+str(int(wind1))+'_s'+s[q]+'_t'+str(int(t[q]))
            wfit2 = 'w'+str(int(wind2))+'_s'+s[q]+'_t'+str(int(t[q]))
            wfit3 = 'w'+str(int(wind3))+'_s'+s[q]+'_t'+str(int(t[q]))
            wfit4 = 'w'+str(int(wind4))+'_s'+s[q]+'_t'+str(int(t[q]))

            basepath = path.join(path.dirname(path.realpath('__file__')), 'data')


            if fit_type == 'vel':
                fdata1 = basepath + path.sep + 'Figshare/VelocityData/'+wfit1+'.csv'
                fdata2 = basepath + path.sep + 'Figshare/VelocityData/'+wfit2+'.csv'
                fdata3 = basepath + path.sep + 'Figshare/VelocityData/'+wfit3+'.csv'
                fdata4 = basepath + path.sep + 'Figshare/VelocityData/'+wfit4+'.csv'
            elif fit_type == 'vort':
                fdata1 = basepath + path.sep + 'Figshare/VorticityData/'+wfit1+'.csv'
                fdata2 = basepath + path.sep + 'Figshare/VorticityData/'+wfit2+'.csv'
                fdata3 = basepath + path.sep + 'Figshare/VorticityData/'+wfit3+'.csv'
                fdata4 = basepath + path.sep + 'Figshare/VorticityData/'+wfit4+'.csv'

            if fit_type == 'vel':
                if read_data == 1:
                    posdn,poslt,datad = starccm_read_vel(np.array([fdata1]),dia,np.array([wind1]),slength[q])
                if read_data == 2:
                    posdn,poslt,datad = starccm_read_vel(np.array([fdata1,fdata2]),dia,np.array([wind1,wind2]),slength[q])
                if read_data == 3:
                    posdn,poslt,datad = starccm_read_vel(np.array([fdata1,fdata2,fdata3]),dia,np.array([wind1,wind2,wind3]),slength[q])
                if read_data == 4:
                    posdn,poslt,datad = starccm_read_vel(np.array([fdata1,fdata2,fdata3,fdata4]),dia,np.array([wind1,wind2,wind3,wind4]),slength[q])
            elif fit_type == 'vort':
                if read_data == 1:
                    posdn,poslt,datad = starccm_read_vort(np.array([fdata1]),dia,np.array([rot1]),slength[q])
                if read_data == 2:
                    posdn,poslt,datad = starccm_read_vort(np.array([fdata1,fdata2]),dia,np.array([rot1,rot2]),slength[q])
                if read_data == 3:
                    posdn,poslt,datad = starccm_read_vort(np.array([fdata1,fdata2,fdata3]),dia,np.array([rot1,rot2,rot3]),slength[q])
                if read_data == 4:
                    posdn,poslt,datad = starccm_read_vort(np.array([fdata1,fdata2,fdata3,fdata4]),dia,np.array([rot1,rot2,rot3,rot4]),slength[q])

            exec('posdns'+sname+'t'+tname+' = posdn')
            exec('poslts'+sname+'t'+tname+' = poslt')
            exec('datads'+sname+'t'+tname+' = datad')

            print 'TSR: '+str(t[q]/100.)+'/solidity: '+s[q]
            q += 1

    print '\nAll STAR-CCM+ data imported'

    if add_data == True:
        #########################################################################################################
        print '\nReading in additional data sets'

        filename = 'NameOfDataFile'
        fdata_add = basepath + path.sep + 'additional_data/'+filename+'.csv'
        dia = 1. # specify the turbine diameter used (m)
        windd = np.array([1.]) # specify the free stream wind speed (m/s)
        lengthpoints = np.array([1.,2.,3.]) # specify the downstream positions of the data (m)

        posdn,poslt,datad = externaldata_read_vel(fdata_add,dia,windd,lengthpoints)

        #########################################################################################################
    elif add_data == False:
        print '\nNo additional data to read in'

    # Cross-validation setup
    if cv_on == 0:
        cvtest = 0.3 # 70-30 split
        scvtr,scvts,tcvtr,tcvts,xttr,xtts,ystr,ysts = train_test_split(scv,tcv,xtcv,yscv,test_size=cvtest)

        # Training Set
        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdn'+name+'tr = posdns'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('poslt'+name+'tr = poslts'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('datad'+name+'tr = datads'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('xt'+name+'tr = xttr[i]')
            exec('ys'+name+'tr = ystr[i]')

        # Validation Set
        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdn'+name+'ts = posdns'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('poslt'+name+'ts = poslts'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('datad'+name+'ts = datads'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('xt'+name+'ts = xtts[i]')
            exec('ys'+name+'ts = ysts[i]')

        # Combining training set into single array for calculations
        posdntr = np.array([])
        poslttr = np.array([])
        datadtr = np.array([])
        xttr = np.array([])
        ystr = np.array([])
        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdntr = np.append(posdntr,posdn'+name+'tr)')
            exec('poslttr = np.append(poslttr,poslt'+name+'tr)')
            exec('datadtr = np.append(datadtr,datad'+name+'tr)')
            exec('xttr = np.append(xttr,np.ones_like(posdn'+name+'tr)*xt'+name+'tr)')
            exec('ystr = np.append(ystr,np.ones_like(posdn'+name+'tr)*ys'+name+'tr)')
        posdntr = np.vstack([posdntr,poslttr])

        # Combining validation set into single array for calculations
        posdnts = np.array([])
        posltts = np.array([])
        datadts = np.array([])
        xtts = np.array([])
        ysts = np.array([])
        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdnts = np.append(posdnts,posdn'+name+'ts)')
            exec('posltts = np.append(posltts,poslt'+name+'ts)')
            exec('datadts = np.append(datadts,datad'+name+'ts)')
            exec('xtts = np.append(xtts,np.ones_like(posdn'+name+'ts)*xt'+name+'ts)')
            exec('ysts = np.append(ysts,np.ones_like(posdn'+name+'ts)*ys'+name+'ts)')
        posdnts = np.vstack([posdnts,posltts])

    # Retraining setup (with all data)
    elif cv_on == 1:
        cvtest = 0.0 # no validation set
        scvtr,scvts,tcvtr,tcvts,xttr,xtts,ystr,ysts = train_test_split(scv,tcv,xtcv,yscv,test_size=cvtest)

        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdn'+name+'tr = posdns'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('poslt'+name+'tr = poslts'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('datad'+name+'tr = datads'+str(int(scvtr[i]))+'t'+str(int(tcvtr[i])))
            exec('xt'+name+'tr = xttr[i]')
            exec('ys'+name+'tr = ystr[i]')

        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdn'+name+'ts = posdns'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('poslt'+name+'ts = poslts'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('datad'+name+'ts = datads'+str(int(scvts[i]))+'t'+str(int(tcvts[i])))
            exec('xt'+name+'ts = xtts[i]')
            exec('ys'+name+'ts = ysts[i]')

        posdntr = np.array([])
        poslttr = np.array([])
        datadtr = np.array([])
        xttr = np.array([])
        ystr = np.array([])
        for i in range(np.size(scvtr)):
            name = str(i+1)
            exec('posdntr = np.append(posdntr,posdn'+name+'tr)')
            exec('poslttr = np.append(poslttr,poslt'+name+'tr)')
            exec('datadtr = np.append(datadtr,datad'+name+'tr)')
            exec('xttr = np.append(xttr,np.ones_like(posdn'+name+'tr)*xt'+name+'tr)')
            exec('ystr = np.append(ystr,np.ones_like(posdn'+name+'tr)*ys'+name+'tr)')
        posdntr = np.vstack([posdntr,poslttr])

        posdnts = np.array([])
        posltts = np.array([])
        datadts = np.array([])
        xtts = np.array([])
        ysts = np.array([])
        for i in range(np.size(scvts)):
            name = str(i+1)
            exec('posdnts = np.append(posdnts,posdn'+name+'ts)')
            exec('posltts = np.append(posltts,poslt'+name+'ts)')
            exec('datadts = np.append(datadts,datad'+name+'ts)')
            exec('xtts = np.append(xtts,np.ones_like(posdn'+name+'ts)*xt'+name+'ts)')
            exec('ysts = np.append(ysts,np.ones_like(posdn'+name+'ts)*ys'+name+'ts)')
        posdnts = np.vstack([posdnts,posltts])

    # Original starting coefficients indicated from Matlab calculations (see Cross Validation Results.xlsx)
    loc10 = np.array([0.08285,-0.05305,-0.1564,0.01036,0.06302,0.08747,-0.0006704,-0.005499,-0.036,0.0])
    loc20 = np.array([-0.2383,0.1042,0.5448,-0.02112,0.05874,-0.8001,0.001587,-0.003793,-0.03121,0.4588])
    loc30 = np.array([0.2494,0.1772,0.312,-0.02348,-0.1658,0.2102,0.0009984,0.01201,0.03067,-0.2078])
    spr10 = np.array([0.04542,-0.01692,-0.2282,0.00144,0.06475,0.1965,0.0,-0.004795,-0.02548,-0.0642])
    spr20 = np.array([-0.07023,0.013,0.1328,-0.0002634,-0.0397,-0.0595,-0.00007737,0.002569,0.0138,0.0])
    skw10 = np.array([-1.663,1.484,0.5703,-0.1898,-5.12,18.05,0.0,0.7082,-0.07986,-12.5])
    skw20 = np.array([-25.8,9.066,51.12,-1.202,-18.56,-52.59,0.05548,0.6758,11.3,14.65])
    scl10 = np.array([-0.217,0.07754,1.49,-0.005692,-0.3125,-1.04,0.0,0.01769,0.1067,0.2558])
    scl20 = np.array([1.97,-1.605,-8.175,0.2362,7.28,-1.009,-0.004933,-0.5333,-2.771,4.26])
    scl30 = np.array([75.5,-13.76,-160.5,1.002,15.25,153.9,-0.02266,-0.4127,-6.459,-53.5])

    # Published results used as starting coefficients
    loc10 = np.array([0.0025703809856661534, -0.0007386258659065129, 0.004595508188667984, 0.000380123563204793, -0.0005090098755683027, 0.005744581813281894, -4.103393770815313e-05, -0.0014146918534486358, -0.013975958482495927, 0.0])
    loc20 = np.array([-0.5047504670963536, 0.23477391362058556, 0.8414256436198028, -0.04252528528617351, -0.06962875967504166, -0.6566907653208429, 0.002839318332370807, 0.00571803958194812, 0.0070744372783060295, 0.22805286438890995])
    loc30 = np.array([0.2878345841026334, 0.11512552658662782, 0.7303949879914625, -0.007035517839387948, -0.18284850673545897, -0.5241921153256568, -0.0003704899921255296, 0.010972527139685873, 0.04380801537377295, 0.1724129349605399])
    spr10 = np.array([0.08234816067475287, -0.03530687906626052, -0.3662863944976986, 0.003240141344532779, 0.12172015102204112, 0.2993048183466721, 0.0, -0.009253185586804007, -0.057469126406649716, -0.07257633583877886])
    spr20 = np.array([-0.07083579909945328, 0.016182024377569406, 0.1985436342461859, 0.0017738254727425816, -0.09111094817943823, -0.06561408122153217, -0.0005115133402638633, 0.009434288536679505, 0.022392136905926813, 0.0])
    skw10 = np.array([-1.6712830849073221, 1.5625053380692426, -6.180392756736983, -0.20407668040293722, -4.6476103643607685, 29.380064536220306, 0.0, 0.7502978877582536, -0.16358232641365608, -19.937609244085568])
    skw20 = np.array([-3.423561091777921, -9.228795430171687, 86.95722105482042, 2.772872601988039, -11.968168333741515, -150.61261090270446, -0.24715316589674527, 0.5283723108899993, 4.537286811245538, 82.50581844010263])
    scl10 = np.array([-0.19815381951708524, 0.08438758133540872, 1.2650146439483734, -0.007606115512168328, -0.2747023984740461, -0.8844640101378567, 0.0, 0.01870057580949183, 0.0699898278743648, 0.2794360008051127])
    scl20 = np.array([2.3932787625531815, -2.020874419612962, -8.938221963838357, 0.576323845480877, 2.8782448498416944, 16.598492450314534, -0.04746016700352029, -0.197101203594028, -1.3860007472886064, -8.289767128060362])
    scl30 = np.array([104.40501489600803, -29.942999569370276, -174.42008279158216, 3.708514822202037, 25.14336546356742, 132.35546551746415, -0.16479555172343271, -1.351556690339512, -6.721810844025761, -40.39565289044579])

    param0 = np.array([])
    param0 = np.append(param0,loc10)
    param0 = np.append(param0,loc20)
    param0 = np.append(param0,loc30)
    param0 = np.append(param0,spr10)
    param0 = np.append(param0,spr20)
    param0 = np.append(param0,skw10)
    param0 = np.append(param0,skw20)
    param0 = np.append(param0,scl10)
    param0 = np.append(param0,scl20)
    param0 = np.append(param0,scl30)

    print '\nStarting curve fitting'

    res,_ = curve_fit(sheet,posdntr,datadtr,p0=param0,maxfev=2820000)

    coef0 = np.array([])
    coef1 = np.array([])
    coef2 = np.array([])
    coef3 = np.array([])
    coef4 = np.array([])
    coef5 = np.array([])
    coef6 = np.array([])
    coef7 = np.array([])
    coef8 = np.array([])
    coef9 = np.array([])

    k = 0
    for i in range(10):
        for j in range(10):
            exec('coef'+str(i)+'= np.append(coef'+str(i)+',res[k])')
            k += 1

    # Output the results of cross-validation
    if cv_on == 0:
        if fit_type == 'vel':
            datac = _vawtwake.sheet_vel(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.,1.,220,200,1)
            dataf = _vawtwake.sheet_vel(xtts,ysts,posdnts[0,:],posdnts[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.,1.,220,200,1)
        elif fit_type == 'vort':
            datac = _vawtwake.sheet_vort(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.)
            dataf = _vawtwake.sheet_vort(xtts,ysts,posdnts[0,:],posdnts[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.)

        error_curve = np.sum((datac-datadtr)**2)
        error_cv = np.sum((dataf-datadts)**2)

        print '\nCross Validation Error:',error_cv
        print 'Curve Fit Error:',error_curve

        print '\ncoef0 = np.array(',coef0.tolist(),')'
        print 'coef1 = np.array(',coef1.tolist(),')'
        print 'coef2 = np.array(',coef2.tolist(),')'
        print 'coef3 = np.array(',coef3.tolist(),')'
        print 'coef4 = np.array(',coef4.tolist(),')'
        print 'coef5 = np.array(',coef5.tolist(),')'
        print 'coef6 = np.array(',coef6.tolist(),')'
        print 'coef7 = np.array(',coef7.tolist(),')'
        print 'coef8 = np.array(',coef8.tolist(),')'
        print 'coef9 = np.array(',coef9.tolist(),')'


        print '\nparam = np.array([',res[0],',',res[1],',',res[2],',',res[3],',',res[4],',',res[5],',',res[6],',',res[7],',',res[8],',',res[9],',\t#loc1'
        print res[10],',',res[11],',',res[12],',',res[13],',',res[14],',',res[15],',',res[16],',',res[17],',',res[18],',',res[19],',\t#loc2'
        print res[20],',',res[21],',',res[22],',',res[23],',',res[24],',',res[25],',',res[26],',',res[27],',',res[28],',',res[29],',\t#loc3'
        print res[30],',',res[31],',',res[32],',',res[33],',',res[34],',',res[35],',',res[36],',',res[37],',',res[38],',',res[39],',\t#spr1'
        print res[40],',',res[41],',',res[42],',',res[43],',',res[44],',',res[45],',',res[46],',',res[47],',',res[48],',',res[49],',\t#spr2'
        print res[50],',',res[51],',',res[52],',',res[53],',',res[54],',',res[55],',',res[56],',',res[57],',',res[58],',',res[59],',\t#skw1'
        print res[60],',',res[61],',',res[62],',',res[63],',',res[64],',',res[65],',',res[66],',',res[67],',',res[68],',',res[69],',\t#skw2'
        print res[70],',',res[71],',',res[72],',',res[73],',',res[74],',',res[75],',',res[76],',',res[77],',',res[78],',',res[79],',\t#scl1'
        print res[80],',',res[81],',',res[82],',',res[83],',',res[84],',',res[85],',',res[86],',',res[87],',',res[88],',',res[89],',\t#scl2'
        print res[90],',',res[91],',',res[92],',',res[93],',',res[94],',',res[95],',',res[96],',',res[97],',',res[98],',',res[99],'])\t#scl3'

    # Output the results of the retraining with all the data
    if cv_on == 1:
        if fit_type == 'vel':
            datac = _vawtwake.sheet_vel(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.,1.,220,200,1)
        elif fit_type == 'vort':
            datac = _vawtwake.sheet_vort(xttr,ystr,posdntr[0,:],posdntr[1,:],coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,1.)

        error_curve = np.sum((datac-datadtr)**2)

        print '\nCurve Fit Error:',error_curve

        print '\ncoef0 = np.array(',coef0.tolist(),')'
        print 'coef1 = np.array(',coef1.tolist(),')'
        print 'coef2 = np.array(',coef2.tolist(),')'
        print 'coef3 = np.array(',coef3.tolist(),')'
        print 'coef4 = np.array(',coef4.tolist(),')'
        print 'coef5 = np.array(',coef5.tolist(),')'
        print 'coef6 = np.array(',coef6.tolist(),')'
        print 'coef7 = np.array(',coef7.tolist(),')'
        print 'coef8 = np.array(',coef8.tolist(),')'
        print 'coef9 = np.array(',coef9.tolist(),')'

        print '\nloc1 = np.array(',coef0.tolist(),')'
        print 'loc2 = np.array(',coef1.tolist(),')'
        print 'loc3 = np.array(',coef2.tolist(),')'
        print 'spr1 = np.array(',coef3.tolist(),')'
        print 'spr2 = np.array(',coef4.tolist(),')'
        print 'skw1 = np.array(',coef5.tolist(),')'
        print 'skw2 = np.array(',coef6.tolist(),')'
        print 'scl1 = np.array(',coef7.tolist(),')'
        print 'scl2 = np.array(',coef8.tolist(),')'
        print 'scl3 = np.array(',coef9.tolist(),')'

        print '\nparam = np.array([',res[0],',',res[1],',',res[2],',',res[3],',',res[4],',',res[5],',',res[6],',',res[7],',',res[8],',',res[9],',\t#loc1'
        print res[10],',',res[11],',',res[12],',',res[13],',',res[14],',',res[15],',',res[16],',',res[17],',',res[18],',',res[19],',\t#loc2'
        print res[20],',',res[21],',',res[22],',',res[23],',',res[24],',',res[25],',',res[26],',',res[27],',',res[28],',',res[29],',\t#loc3'
        print res[30],',',res[31],',',res[32],',',res[33],',',res[34],',',res[35],',',res[36],',',res[37],',',res[38],',',res[39],',\t#spr1'
        print res[40],',',res[41],',',res[42],',',res[43],',',res[44],',',res[45],',',res[46],',',res[47],',',res[48],',',res[49],',\t#spr2'
        print res[50],',',res[51],',',res[52],',',res[53],',',res[54],',',res[55],',',res[56],',',res[57],',',res[58],',',res[59],',\t#skw1'
        print res[60],',',res[61],',',res[62],',',res[63],',',res[64],',',res[65],',',res[66],',',res[67],',',res[68],',',res[69],',\t#skw2'
        print res[70],',',res[71],',',res[72],',',res[73],',',res[74],',',res[75],',',res[76],',',res[77],',',res[78],',',res[79],',\t#scl1'
        print res[80],',',res[81],',',res[82],',',res[83],',',res[84],',',res[85],',',res[86],',',res[87],',',res[88],',',res[89],',\t#scl2'
        print res[90],',',res[91],',',res[92],',',res[93],',',res[94],',',res[95],',',res[96],',',res[97],',',res[98],',',res[99],'])\t#scl3'
