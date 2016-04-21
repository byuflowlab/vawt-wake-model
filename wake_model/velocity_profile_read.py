import csv
import numpy as np
import matplotlib.pyplot as plt


def velprof(wfit,freestream,norm):

    # wfit = 's2_400'
    # freestream = 15.
    
    
    #Initializing matrices
    pos1 = np.array([])
    pos2 = np.array([])
    pos3 = np.array([])
    pos4 = np.array([])
    pos5 = np.array([])
    pos6 = np.array([])
    pos7 = np.array([])
    pos8 = np.array([])
    pos9 = np.array([])
    pos10 = np.array([])
    pos15 = np.array([])
    pos20 = np.array([])
    pos25 = np.array([])
    pos30 = np.array([])
    vel1 = np.array([])
    vel2 = np.array([])
    vel3 = np.array([])
    vel4 = np.array([])
    vel5 = np.array([])
    vel6 = np.array([])
    vel7 = np.array([])
    vel8 = np.array([])
    vel9 = np.array([])
    vel10 = np.array([])
    vel15 = np.array([])
    vel20 = np.array([])
    vel25 = np.array([])
    vel30 = np.array([])
    
    
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
            pos15 = np.append(pos15,float(row[20]))
        if row[22] != 'null' and i != 0:
            pos20 = np.append(pos20,float(row[22]))
        if row[24] != 'null' and i != 0:
            pos25 = np.append(pos25,float(row[24]))
        if row[26] != 'null' and i != 0:
            pos30 = np.append(pos30,float(row[26]))
        
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
            vel15 = np.append(vel15,float(row[21]))
        if row[23] != 'null' and i != 0:
            vel20 = np.append(vel20,float(row[23]))
        if row[25] != 'null' and i != 0:
            vel25 = np.append(vel25,float(row[25]))
        if row[27] != 'null' and i != 0:
            vel30 = np.append(vel30,float(row[27]))
        i += 1
    
    f.close()
    
    if norm == True:
        #Calculating normalized averages
        vel1 = vel1/freestream
        vel2 = vel2/freestream
        vel3 = vel3/freestream
        vel4 = vel4/freestream
        vel5 = vel5/freestream
        vel6 = vel6/freestream
        vel7 = vel7/freestream
        vel8 = vel8/freestream
        vel9 = vel9/freestream
        vel10 = vel10/freestream
        vel15 = vel15/freestream
        vel20 = vel20/freestream
        vel25 = vel25/freestream
        vel30 = vel30/freestream
    
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
    
    
