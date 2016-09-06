import csv
from os import path
import numpy as np
from numpy import fabs
from scipy.interpolate import RectBivariateSpline,griddata


def vorticity(tsr, solidity):
    """
    Using EMG distribution parameters to define the vorticity strength and shape
    
    Parameters
    ----------
    tsr : float
        tip-speed ratio
    solidity : float
        turbine solidity
    
    Returns
    ----------
    loc : array
        array of the location parameter (3 values)
    spr : array
        array of the spread parameter (2 values)
    skw : array
        array of the skew parameter (2 values)
    scl : array
        array of the scale parameter (3 values)
    """
    
    # Reading in csv file (vorticity database)
    basepath = path.join(path.dirname(path.realpath(__file__)), 'data')
    fdata = basepath + path.sep + 'vortdatabase.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw,0)
            vortdat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d, float(row[1]))  # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            vortdat = np.vstack([vortdat, raw])  # adding entry to vorticity database array
        i += 1
    f.close()
    
    vortdat = np.delete(vortdat, 0, axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float)  # converting tip-speed ratio entries into floats
    vortdat = vortdat.astype(np.float)  # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_loc1 = vortdat[i*10]\ns'+sol+'_loc2 = vortdat[i*10+1]\ns'+sol+'_loc3 = vortdat[i*10+2]\ns'+sol+'_spr1 = vortdat[i*10+3]\ns'+sol+'_spr2 = vortdat[i*10+4]\ns'+sol+'_skw1 = vortdat[i*10+5]\ns'+sol+'_skw2 = vortdat[i*10+6]\ns'+sol+'_scl1 = vortdat[i*10+7]\ns'+sol+'_scl2 = vortdat[i*10+8]\ns'+sol+'_scl3 = vortdat[i*10+9]\n')
    
    # BIVARIATE SPLINE INTERPOLATION
    
    iz = np.size(sol_d)
    jz = np.size(tsr_d)
    
    # Initializing rectangular matrices
    z_loc1 = np.zeros((iz, jz))
    z_loc2 = np.zeros((iz, jz))
    z_loc3 = np.zeros((iz, jz))
    z_spr1 = np.zeros((iz, jz))
    z_spr2 = np.zeros((iz, jz))
    z_skw1 = np.zeros((iz, jz))
    z_skw2 = np.zeros((iz, jz))
    z_scl1 = np.zeros((iz, jz))
    z_scl2 = np.zeros((iz, jz))
    z_scl3 = np.zeros((iz, jz))
    
    # Transferring raw data into rectangular matrices
    for i in range(iz):
        for j in range(jz):
            sol = str(i+1)
            exec('z_loc1[i, j] = s'+sol+'_loc1[j]')
            exec('z_loc2[i, j] = s'+sol+'_loc2[j]')
            exec('z_loc3[i, j] = s'+sol+'_loc3[j]')
            exec('z_spr1[i, j] = s'+sol+'_spr1[j]')
            exec('z_spr2[i, j] = s'+sol+'_spr2[j]')
            exec('z_skw1[i, j] = s'+sol+'_skw1[j]')
            exec('z_skw2[i, j] = s'+sol+'_skw2[j]')
            exec('z_scl1[i, j] = s'+sol+'_scl1[j]')
            exec('z_scl2[i, j] = s'+sol+'_scl2[j]')
            exec('z_scl3[i, j] = s'+sol+'_scl3[j]')
    
    # Creating a rectangular bivariate spline of the parameter data
    s_loc1 = RectBivariateSpline(sol_d, tsr_d, z_loc1)
    s_loc2 = RectBivariateSpline(sol_d, tsr_d, z_loc2)
    s_loc3 = RectBivariateSpline(sol_d, tsr_d, z_loc3)
    s_spr1 = RectBivariateSpline(sol_d, tsr_d, z_spr1)
    s_spr2 = RectBivariateSpline(sol_d, tsr_d, z_spr2)
    s_skw1 = RectBivariateSpline(sol_d, tsr_d, z_skw1)
    s_skw2 = RectBivariateSpline(sol_d, tsr_d, z_skw2)
    s_scl1 = RectBivariateSpline(sol_d, tsr_d, z_scl1)
    s_scl2 = RectBivariateSpline(sol_d, tsr_d, z_scl2)
    s_scl3 = RectBivariateSpline(sol_d, tsr_d, z_scl3)
    
    # Selecting the specific parameters to use for TSR and solidity
    loc1 = s_loc1(solidity, tsr)
    loc2 = s_loc2(solidity, tsr)
    loc3 = s_loc3(solidity, tsr)
    spr1 = s_spr1(solidity, tsr)
    spr2 = s_spr2(solidity, tsr)
    skw1 = s_skw1(solidity, tsr)
    skw2 = s_skw2(solidity, tsr)
    scl1 = s_scl1(solidity, tsr)
    scl2 = s_scl2(solidity, tsr)
    scl3 = s_scl3(solidity, tsr)
    
    # Creating arrays of the parameters
    loc = np.array([loc1[0, 0], loc2[0, 0], loc3[0, 0]])
    spr = np.array([spr1[0, 0], spr2[0, 0]])
    skw = np.array([skw1[0, 0], skw2[0, 0]])
    scl = np.array([scl1[0, 0], scl2[0, 0], scl3[0, 0]])
    
    return loc, spr, skw, scl


def velocity(tsr, solidity):
    """
    Using SMG distribution parameters to define the velocity strength and shape
    
    Parameters
    ----------
    tsr : float
        tip-speed ratio
    solidity : float
        turbine solidity
    
    Returns
    ----------
    men : array
        array of the mean parameter (1 values)
    sdv : array
        array of the standard deviation parameter (4 values)
    rat : array
        array of the rate parameter (1 values)
    wdt : array
        array of the translation parameter (1 values)
    spr : array
        array of the spread parameter (4 values)
    scl : array
        array of the scale parameter (3 values)
    """
    # Reading in csv file (vorticity database)
    basepath = path.join(path.dirname(path.realpath(__file__)), 'data')
    # fdata = basepath + path.sep + 'velodatabase_SMG_surf2_edit.csv'
    fdata = basepath + path.sep + 'velodatabase_SMG_surf4.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw, 0)
            velodat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d, float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw, 0)
            velodat = np.vstack([velodat, raw]) # adding entry to vorticity database array
        i += 1
    f.close()

    velodat = np.delete(velodat, 0, axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    velodat = velodat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each SMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_men = velodat[i*14]\ns'+sol+'_sdv1 = velodat[i*14+1]\ns'+sol+'_sdv2 = velodat[i*14+2]\ns'+sol+'_sdv3 = velodat[i*14+3]\ns'+sol+'_sdv4 = velodat[i*14+4]\ns'+sol+'_rat = velodat[i*14+5]\ns'+sol+'_wdt = velodat[i*14+6]\ns'+sol+'_spr1 = velodat[i*14+7]\ns'+sol+'_spr2 = velodat[i*14+8]\ns'+sol+'_spr3 = velodat[i*14+9]\ns'+sol+'_spr4 = velodat[i*14+10]\ns'+sol+'_scl1 = velodat[i*14+11]\ns'+sol+'_scl2 = velodat[i*14+12]\ns'+sol+'_scl3 = velodat[i*14+13]\n')

    # NEAREST ND INTERPOLATION

    iz = np.size(sol_d)
    jz = np.size(tsr_d)

    # Transferring parameter arrays to data matrix
    dataset = np.zeros((14*iz,jz))

    for i in range(iz):
        sol = str(i+1)
        exec('dataset[i*14] = s'+sol+'_men')
        exec('dataset[i*14+1] = s'+sol+'_sdv1')
        exec('dataset[i*14+2] = s'+sol+'_sdv2')
        exec('dataset[i*14+3] = s'+sol+'_sdv3')
        exec('dataset[i*14+4] = s'+sol+'_sdv4')
        exec('dataset[i*14+5] = s'+sol+'_rat')
        exec('dataset[i*14+6] = s'+sol+'_wdt')
        exec('dataset[i*14+7] = s'+sol+'_spr1')
        exec('dataset[i*14+8] = s'+sol+'_spr2')
        exec('dataset[i*14+9] = s'+sol+'_spr3')
        exec('dataset[i*14+10] = s'+sol+'_spr4')
        exec('dataset[i*14+11] = s'+sol+'_scl1')
        exec('dataset[i*14+12] = s'+sol+'_scl2')
        exec('dataset[i*14+13] = s'+sol+'_scl3')

    # Transferring data matrix to TSR, soldity, and parameter arrays
    points1 = np.array([])
    points2 = np.array([])
    points3 = np.array([])
    values = np.array([])
    for i in range(jz):
        for j in range(iz):
            for k in range(14):
                points1 = np.append(points1,tsr_d[i])
                points2 = np.append(points2,sol_d[j])
                points3 = np.append(points3,k+1)
                values = np.append(values,dataset[j*14+k,i])

    # Identifying nearest ND point to CFD data set
    met = 'nearest'
    men = griddata((points1,points2,points3),values,(tsr,solidity,1),method=met)
    sdv1 = griddata((points1,points2,points3),values,(tsr,solidity,2),method=met)
    sdv2 = griddata((points1,points2,points3),values,(tsr,solidity,3),method=met)
    sdv3 = griddata((points1,points2,points3),values,(tsr,solidity,4),method=met)
    sdv4 = griddata((points1,points2,points3),values,(tsr,solidity,5),method=met)
    rat = griddata((points1,points2,points3),values,(tsr,solidity,6),method=met)
    wdt = griddata((points1,points2,points3),values,(tsr,solidity,7),method=met)
    spr1 = griddata((points1,points2,points3),values,(tsr,solidity,8),method=met)
    spr2 = griddata((points1,points2,points3),values,(tsr,solidity,9),method=met)
    spr3 = griddata((points1,points2,points3),values,(tsr,solidity,10),method=met)
    spr4 = griddata((points1,points2,points3),values,(tsr,solidity,11),method=met)
    scl1 = griddata((points1,points2,points3),values,(tsr,solidity,12),method=met)
    scl2 = griddata((points1,points2,points3),values,(tsr,solidity,13),method=met)
    scl3 = griddata((points1,points2,points3),values,(tsr,solidity,14),method=met)

    # Identifying nearest TSR and solidity to given variables
    pointst = np.copy(points1)/tsr - 1.
    pointss = np.copy(points2)/solidity - 1.
    ti = np.argmin(fabs(pointst))
    si = np.argmin(fabs(pointss))
    tsrn = points1[ti]
    soln = points2[si]

    # Creating arrays of the parameters
    men = np.array([men])
    sdv = np.array([sdv1, sdv2, sdv3, sdv4])
    rat = np.array([rat])
    wdt = np.array([wdt])
    spr = np.array([spr1, spr2, spr3, spr4])
    scl = np.array([scl1, scl2, scl3])

    return men,sdv,rat,wdt,spr,scl,tsrn,soln


def overlay(xt,ys,coef):
    a = coef[0]
    b = coef[1]
    c = coef[2]
    d = coef[3]
    e = coef[4]
    f = coef[5]
    g = coef[6]
    h = coef[7]
    i = coef[8]
    j = coef[9]

    return a + b*xt + c*ys + d*xt**2 + e*xt*ys + f*ys**2 + g*xt**3 + h*xt**2*ys + i*xt*ys**2 + j*ys**3

def velocity2(tsr,solidity):
    coef0 = np.array( [69.31120214165799, -14.54226846455916, -182.9484314845766, 0.7757500007593103, 17.077128640459204, 205.7405780507834, 0.0, 0.6952004161258284, -8.656278871516317, -80.30651611713891] )
    coef1 = np.array( [20.796591097382365, -6.244181130387176, -8.587693125209208, 0.8142615941603905, 1.8506138334190578, 0.0, -0.03363107177764374, 0.12223019591467849, 0.0, 0.0] )
    coef2 = np.array( [-0.015294073403020714, 0.020473995751279686, 0.13789844492064895, -0.0052638657413030055, -0.11469080246898106, -0.19707130831185676, 0.00028907835977146215, 0.020435994570160664, 0.12862824866464892, 0.0] )
    coef3 = np.array( [15.938353497212011, -7.420923601227339, -17.81286649008919, 1.2062417608423517, 3.955459751836299, 7.241485498819938, -0.06694101282588388, -0.21226646110737413, -0.7478822202018863, 0.0] )
    coef4 = np.array( [0.2619104899556741, -0.073703248730009, -0.18110278053024495, 0.0, 0.08236981420749905, -0.07033060400290583, 0.0, 0.0, 0.0, 0.0] )
    coef5 = np.array( [1.2772126399139734, -0.28793513736493115, -0.9675440002861441, 0.022232873411874163, 0.16702810476509866, 0.21867217108386947, -1.8776159079600002e-05, -0.004624278672205385, 0.03751786152363546, -0.060347089708739324] )
    coef6 = np.array( [0.07310683619403949, -0.0731078636275017, -0.37932110788073947, 0.017633578369773174, 0.2769037348921203, 0.3208692982952059, -0.0011170884942625948, -0.015307930271091856, -0.08194905828120999, 0.0] )
    coef7 = np.array( [83.62824233391731, -20.14985073380763, -119.3133697945351, 2.186060358249593, 14.877953353208454, 82.12647543800952, -0.08428765888226167, -0.7381718005397198, -3.2357192642279067, -24.477809412902253] )


    coef0 = np.array( [96.78787569603121, -25.958433685317107, -264.69170623326727, 1.87464555787568, 51.87562056487128, 239.2539339088309, 0.0, -2.3683807739615252, -21.77566460245587, -72.16807143827123] )
    coef1 = np.array( [-1.5607308407855849, 3.781810622756087, 3.8589765773806635, -0.7828322617478133, 1.2851127289130426, -0.9780850425505745, 0.053806398646338315, -0.004216458910203603, -1.9029882519009809, 1.1216497213945278] )
    coef2 = np.array( [0.04474980994165066, -0.01781338292690174, -0.2928285095380119, 0.0022556558140634383, 0.05240119175487355, 0.42156882995131023, -0.0001656666942084075, 0.004380334411217462, 0.0017623529889020642, -0.22543995985088924] )
    coef3 = np.array( [0.44456289023552864, -0.15060882084694283, -0.5839444598651449, 0.017653548310866428, 0.13144967434649582, 0.25020074597304015, -0.0006767773859921913, -0.007996991123069315, -0.02796902183745419, -0.031964165183268466] )
    coef4 = np.array( [-2.7247813990305123, 4.724192699399499, 10.36801276589125, -1.2166148524625373, -7.9280796594054115, 0.0, 0.08822977083653992, 0.677747763097038, 0.0, 0.0] )
    coef5 = np.array( [0.7506898031843403, -0.06444761119603926, -0.837579404198149, -0.009267404487624946, 0.19446504777032012, 0.09424648444431204, 0.0013092828164535946, -0.004935762159895988, -0.03629240478785554, 0.20951216903837414] )
    coef6 = np.array( [0.14354129036499289, -0.0713909071063587, -0.6243972531747024, 0.011973080714083267, 0.30933176461872885, 0.6699155214984721, -0.0004766541652853193, -0.018026755908751073, -0.08464726862880464, -0.20627751665213928] )
    coef7 = np.array( [41.295014624752916, -4.660859316203044, -39.90052730458264, 0.19616479682310406, 1.7128905012718616, 13.61901360510989, 0.0, 0.0, 0.0, 0.0] )


    coef0 = np.array( [132.6722372610478, -34.93465942596633, -380.4394325386311, 2.4459067170844673, 68.60880512259655, 374.49565996499416, 0.0, -2.631298457155278, -32.41512006231126, -122.91288763921081] )
    coef1 = np.array( [3.2841496359202185, 2.5258185980166776, -19.32475137654605, -0.6635605862665878, 4.39276994866165, 31.457505087283618, 0.05362178673494203, -0.1291054894957141, -3.8167603944580875, -13.111805505269034] )
    coef2 = np.array( [0.16041249254968282, -0.03660737924109727, -1.3324697946077309, 0.0006952312900736275, 0.26518562270588875, 1.8519148426513374, 0.00010378689556150557, -0.00454873844845595, -0.14497362526941746, -0.7847581114747488] )
    coef3 = np.array( [0.5129822083425051, -0.20140014304589243, -0.25138976836980925, 0.02413866231871476, 0.13084332587949254, -0.2927441073541505, -0.0008954524910658225, -0.009407338557251597, -0.011865207041416847, 0.1874098605524512] )
    coef4 = np.array( [-2.260178381037694, 4.148684620712865, 6.231819496294992, -1.1073679523291522, -5.638070571202504, 0.0, 0.07174203727980469, 0.656402421601931, 0.0, 0.0] )
    coef5 = np.array( [1.0333574828116667, -0.18231961505300875, -1.3077494182082738, 0.007527942821932585, 0.271193222135317, 0.6347854281130164, 0.00046651250333707664, -0.00603448454892895, -0.08287667731283102, 0.008341552496059878] )
    coef6 = np.array( [0.06105811054700276, -0.0441101309198774, -0.2765502943108469, 0.009849233044716597, 0.23043914347423136, 0.27613892118545386, -0.0005264583562888652, -0.012695069285973236, -0.05490202811295547, -0.05645025865265032] )
    coef7 = np.array( [39.95802886172504, -4.656227200127483, -37.89948079419315, 0.22772918100271503, 1.337720465366513, 14.171027502911437, 0.0, 0.0, 0.0, 0.0] )

    coef0 = np.array( [110.36030047488595, -27.985029641159635, -317.9413884301561, 1.878631951546971, 56.098972269155766, 312.4289628464063, 0.0, -2.0059698573093683, -27.32200171079657, -100.17261139540484] )
    coef1 = np.array( [1.1862973344821859, 3.2913660782457277, -13.343558670695955, -0.6994710704769485, 2.755577535790919, 24.41590484609918, 0.04455378755559768, 0.07032158254873834, -3.7862576950238216, -8.720769026464755] )
    coef2 = np.array( [0.1255751578337351, -0.023873145621841157, -1.1731644110837747, 0.0, 0.21525685130206756, 1.6926888990773763, 0.0, 0.0, -0.13006998640696732, -0.7203847208910938] )
    coef3 = np.array( [0.604604570567267, -0.23669313806774606, -0.34190928626245876, 0.028525238398903852, 0.1626346287030023, -0.3025727427349951, -0.0010750324435612887, -0.011521887892891916, -0.018045156540234433, 0.21032697191920527] )
    coef4 = np.array( [-2.1579735618906395, 3.8231906624406182, 6.203244641496877, -0.9977766045457122, -5.480607527455324, 0.0, 0.06455330726056291, 0.6226831089712495, 0.0, 0.0] )
    coef5 = np.array( [0.9956179768896047, -0.1851658632946506, -1.1667623179247888, 0.010377359343377409, 0.2551128576510571, 0.5398026250286322, 0.0002202373562505837, -0.006404037687646303, -0.07055287805082666, 0.0] )
    coef6 = np.array( [0.042251377261527216, -0.03314023685276692, -0.2132999223735087, 0.007497085100377211, 0.21850016639856837, 0.16869202696085084, -0.0003310698528521071, -0.012647244265627924, -0.04544937719634106, 0.0] )
    coef7 = np.array( [39.34005307348499, -4.563877850653355, -36.29775209491077, 0.215754042224762, 1.3565208067720742, 12.914792852294092, 0.0, 0.0, 0.0, 0.0] )

    spr1 = overlay(tsr,solidity,coef0)
    pow1 = overlay(tsr,solidity,coef1)
    pow2 = overlay(tsr,solidity,coef2)
    spr2 = overlay(tsr,solidity,coef3)
    skw = overlay(tsr,solidity,coef4)
    scl1 = overlay(tsr,solidity,coef5)
    scl2 = overlay(tsr,solidity,coef6)
    scl3 = overlay(tsr,solidity,coef7)

    return spr1,pow1,pow2,spr2,skw,scl1,scl2,scl3