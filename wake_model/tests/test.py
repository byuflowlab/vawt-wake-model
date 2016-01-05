# Unit test for VAWT_Wake_Model using PIV and wind tunnel experiments as references

import unittest
import numpy as np
from scipy.io import loadmat
from VWM_Fortran import velocity_field

class Testwakemodel(unittest.TestCase):
    def test_PIV(self):

        # Importing PIV test data
        tesdata = loadmat('dataH_VAWT.mat')
        
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
        for i in range(np.size(x15)):
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
        
        for i in range(np.size(x40)):
            x40[i] = (tesdata['y'][int(i*space+7.*space),2956])/1000.
            y40[i] = (tesdata['u_x'][int(i*space+7.*space),2956])/9.3
        
        # Wake model calculation using PIV data specifications
        rad = 0.5
        dia = 2*rad
        velf = 9.308422677
        sol = 0.24
        tsr = 4.5
        rom15 = np.zeros_like(x15)
        rom20 = np.zeros_like(x20)
        rom25 = np.zeros_like(x25)
        rom30 = np.zeros_like(x30)
        rom35 = np.zeros_like(x35)
        rom40 = np.zeros_like(x40)
        for i in range(np.size(x15)):
            rom15[i] = velocity_field(0.75,x15[i]*dia,velf,dia,tsr,sol)
            rom20[i] = velocity_field(1.0,x20[i]*dia,velf,dia,tsr,sol)
            rom25[i] = velocity_field(1.25,x25[i]*dia,velf,dia,tsr,sol)
            rom30[i] = velocity_field(1.5,x30[i]*dia,velf,dia,tsr,sol)
            rom35[i] = velocity_field(1.75,x35[i]*dia,velf,dia,tsr,sol)
        for i in range(np.size(x40)):
            rom40[i] = velocity_field(2.0,x40[i]*dia,velf,dia,tsr,sol)
            
            
        idx1 = np.array([4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
        idx2 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25])
        np.testing.assert_allclose(rom15[idx1],y15[idx1],atol=0.4)
        np.testing.assert_allclose(rom20[idx1],y20[idx1],atol=0.4)
        np.testing.assert_allclose(rom25[idx1],y25[idx1],atol=0.4)
        np.testing.assert_allclose(rom30[idx1],y30[idx1],atol=0.4)
        np.testing.assert_allclose(rom35[idx1],y35[idx1],atol=0.4)
        np.testing.assert_allclose(rom40[idx2],y40[idx2],atol=0.6)
        
        
        
        # Wind tunnel test data
        x15wt = np.linspace(-1.7,1.7,21)
        y15wt = np.array([0.98,0.98,0.985,0.99,0.995,1.005,0.96,0.73,0.5,0.44,0.54,0.51,0.66,0.9,1.01,1.0,0.99,0.985,0.98,0.98,0.97])
        
        # Wake model calculation using wind tunnel data specifications
        rad = 0.515
        dia = 2*rad
        velf = 16.0
        sol = 0.5
        tsr = 1.6
        rom15wt = np.zeros_like(x15wt)
        for i in range(np.size(rom15wt)):
            rom15wt[i] = velocity_field(1.5*dia,x15wt[i]*dia,velf,dia,tsr,sol)
        
        np.testing.assert_allclose(rom15wt,y15wt,atol=0.4)
        
if __name__ == '__main__':
    unittest.main(exit=False)
        