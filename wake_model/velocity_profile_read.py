import csv
import numpy as np
import matplotlib.pyplot as plt


def velprof(wfit,freestream,norm):

    # wfit = 's2_400'
    # freestream = 15.
    
    
    #Initializing matrices
    for k in range(30):
        name = str(k+1)
        exec('pos'+name+' = np.array([])')
        exec('velo'+name+' = np.array([])')
    
    
    #Files read in
    # fdata = 'C:\\Users\\TingeyPC\\Documents\\zStar-CCM\NACA0021\\MoveForward\\Velocity Sections\\'+wfit+'.csv'
    # C:\Users\TingeyPC\Documents\zStar-CCM\NACA0021\MoveForward\Velocity Sections
    fdata = '/Users/ning1/Documents/FLOW Lab/STAR-CCM+/NACA0021/MoveForward/Velocity Sections/'+wfit+'.csv'
    
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
    
    if norm == True:
        #Calculating normalized averages
        vel1 = velo1/freestream
        vel2 = velo2/freestream
        vel3 = velo3/freestream
        vel4 = velo4/freestream
        vel5 = velo5/freestream
        vel6 = velo6/freestream
        vel7 = velo7/freestream
        vel8 = velo8/freestream
        vel9 = velo9/freestream
        vel10 = velo10/freestream
        vel15 = velo15/freestream
        vel20 = velo20/freestream
        vel25 = velo25/freestream
        vel30 = velo30/freestream
    
    return vel1,vel2,vel3,vel4,vel5,vel6,vel7,vel8,vel9,vel10,vel15,vel20,vel25,vel30,pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9,pos10,pos15,pos20,pos25,pos30
    
    
## Main File
if __name__ == "__main__":
    
    wfit = 's2_400'
    freestream = 15.
    norm = True
    # norm = False
    
    vel1,vel2,vel3,vel4,vel5,vel6,vel7,vel8,vel9,vel10,vel15,vel20,vel25,vel30,pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9,pos10,pos15,pos20,pos25,pos30 = velprof(wfit,freestream,norm)
    
    
    #Plotting the results
    plt.figure()
    plt.plot(vel2,pos2,'m.',markersize=5,label='x=2D')
    plt.plot(vel4,pos4,'b.',markersize=5,label='x=4D')
    plt.plot(vel6,pos6,'c.',markersize=5,label='x=6D')
    plt.plot(vel8,pos8,'g.',markersize=5,label='x=8D')
    plt.plot(vel10,pos10,'y.',markersize=5,label='x=10D')
    plt.plot(vel15,pos15,'r.',markersize=5,label='x=15D')
    plt.xlabel('Velocity/Free Stream')
    plt.ylabel('Position (m)')
    plt.legend(loc=2)
    plt.title('Velocity Profile: TSR '+wfit[3]+'.'+wfit[4]+wfit[5])
    plt.axis([0.0,1.5,-51.0,51.0])
    plt.show()
    
    
