from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp,sqrt,pi,fabs
from VAWT_Wake_Model import velocity_field
from database_call import velocity
import csv
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

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



coef0 = np.array( [0.0025703809856661534, -0.0007386258659065129, 0.004595508188667984, 0.000380123563204793, -0.0005090098755683027, 0.005744581813281894, -4.103393770815313e-05, -0.0014146918534486358, -0.013975958482495927, 0.0] )
coef1 = np.array( [-0.5047504670963536, 0.23477391362058556, 0.8414256436198028, -0.04252528528617351, -0.06962875967504166, -0.6566907653208429, 0.002839318332370807, 0.00571803958194812, 0.0070744372783060295, 0.22805286438890995] )
coef2 = np.array( [0.2878345841026334, 0.11512552658662782, 0.7303949879914625, -0.007035517839387948, -0.18284850673545897, -0.5241921153256568, -0.0003704899921255296, 0.010972527139685873, 0.04380801537377295, 0.1724129349605399] )
coef3 = np.array( [0.08234816067475287, -0.03530687906626052, -0.3662863944976986, 0.003240141344532779, 0.12172015102204112, 0.2993048183466721, 0.0, -0.009253185586804007, -0.057469126406649716, -0.07257633583877886] )
coef4 = np.array( [-0.07083579909945328, 0.016182024377569406, 0.1985436342461859, 0.0017738254727425816, -0.09111094817943823, -0.06561408122153217, -0.0005115133402638633, 0.009434288536679505, 0.022392136905926813, 0.0] )
coef5 = np.array( [-1.6712830849073221, 1.5625053380692426, -6.180392756736983, -0.20407668040293722, -4.6476103643607685, 29.380064536220306, 0.0, 0.7502978877582536, -0.16358232641365608, -19.937609244085568] )
coef6 = np.array( [-3.423561091777921, -9.228795430171687, 86.95722105482042, 2.772872601988039, -11.968168333741515, -150.61261090270446, -0.24715316589674527, 0.5283723108899993, 4.537286811245538, 82.50581844010263] )
coef7 = np.array( [-0.19815381951708524, 0.08438758133540872, 1.2650146439483734, -0.007606115512168328, -0.2747023984740461, -0.8844640101378567, 0.0, 0.01870057580949183, 0.0699898278743648, 0.2794360008051127] )
coef8 = np.array( [2.3932787625531815, -2.020874419612962, -8.938221963838357, 0.576323845480877, 2.8782448498416944, 16.598492450314534, -0.04746016700352029, -0.197101203594028, -1.3860007472886064, -8.289767128060362] )
coef9 = np.array( [104.40501489600803, -29.942999569370276, -174.42008279158216, 3.708514822202037, 25.14336546356742, 132.35546551746415, -0.16479555172343271, -1.351556690339512, -6.721810844025761, -40.39565289044579] )

# coef0 = np.array( [0.07709044464625317, -0.04740029165688312, -0.16991048271258233, 0.010185705350301196, 0.06884032150568876, 0.10299017833087143, -0.0007122943215019527, -0.008308822666077508, -0.03435156111077655, 0.0] )
# coef1 = np.array( [-0.5003717932946304, 0.22610110226599275, 1.3085708132849527, -0.04811681050703861, -0.15715350741292458, -1.2771809638703127, 0.003860117612383805, 0.0128260317514925, 0.07229701122051288, 0.38505299579079794] )
# coef2 = np.array( [0.18447471588496925, 0.22491446110851412, 0.352882721970498, -0.022220685661704473, -0.26163657606654983, 0.35858243211253543, -0.0002849994734101372, 0.02441732290603862, 0.0017870851288199138, -0.15083295080322173] )
# coef3 = np.array( [0.07916577329587404, -0.03144757083755946, -0.340832651135764, 0.002728282634468151, 0.10713756080183831, 0.25582398671769246, 0.0, -0.007915289453241086, -0.05011301604845094, -0.044729906209416866] )
# coef4 = np.array( [-0.08562759918040584, 0.014952261574604832, 0.1043146634725098, 0.0027479619455810732, -0.03298914821394256, 0.014119501371758644, -0.0006302605846889573, 0.004795277656533713, -0.0139277398842217, 0.0] )
# coef5 = np.array( [-1.655460343720484, 1.7367941271405778, -7.238562586343561, -0.23492155695145236, -4.3949689533986, 32.148612999217015, 0.0, 0.7643482585274446, -0.02223360684222774, -22.812878863992015] )
# coef6 = np.array( [-7.889402871093269, -7.480765809247869, 81.43536734767329, 2.5259032797143357, -12.91777745115153, -130.9655557181711, -0.2273087970688141, 0.46713593707455003, 3.5290597349608985, 75.13101144264893] )
# coef7 = np.array( [-0.281860241228644, 0.15337593599273622, 1.3874156400015434, -0.013349346689954193, -0.4499765244970901, -0.6531816298972863, 0.0, 0.030394272428236082, 0.13639777182324817, 0.06777399286119161] )
# coef8 = np.array( [3.107813107183466, -1.9460558282949345, -15.241180687382021, 0.16298098950801992, 8.55044346091066, 12.09032038574375, 2.9217309255614916e-05, -0.6535225876034279, -3.0377325606705443, -3.523374903798156] )
# coef9 = np.array( [12531.99583129866, -40.94972927053573, -53958.89757682066, 5.918882881361156, 28.034407867118492, 74594.37812720957, -0.28546799520575034, -2.0428378590399507, -4.530858713114041, -33138.3719592977] )



coefarray = np.array([coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9])

for k in range(10):
    n = 50
    x = np.linspace(1.5,7.,n)
    y = np.linspace(0.15,1.0,n)


    fs = 15
    fig = plt.figure(k+1,figsize=(10,6))
    fig.subplots_adjust(bottom=0,left=0,right=1,top=1)
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(x, y)
    G = np.zeros((n,n))
    H = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            G[i,j] = overlay(X[i,j],Y[i,j],coefarray[k])
            # H[i,j] = overlay(X[i,j],Y[i,j],coef4)

    surf = ax.plot_surface(X, Y, G, rstride=1, cstride=1, color='g', linewidth=0, antialiased=True, alpha=0.5, label='SMG Surface')#cmap=cm.parula
    # surf = ax.plot_surface(X, Y, H*100., rstride=1, cstride=1, color='b', linewidth=0, antialiased=True, alpha=0.5, label='SMG Surface')#cmap=cm.parula
    surf = ax.plot_surface(X, Y, np.ones_like(G)*0., rstride=1, cstride=1, color='r', linewidth=0, antialiased=True, alpha=0.5, label='SMG Surface')#cmap=cm.parula
    # ax.azim = -160
    # ax.elev = 35
    ax.set_xlabel('TSR',fontsize=fs)
    ax.set_ylabel('Solidity',fontsize=fs)
    ax.set_zlabel('Value',fontsize=fs)
    # ax.set_xlim(0,20)
    # ax.set_ylim(-2,2)
    # ax.set_zlim(0,1.2)
    # ax.set_xticklabels(np.array([0,5,10,15,20]),fontsize=fs)
    # ax.set_yticklabels(np.array([-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0]),fontsize=fs)
    # ax.set_zticklabels(np.array([0.0,0.2,0.4,0.6,0.8,1.0]),fontsize=fs)



plt.show()