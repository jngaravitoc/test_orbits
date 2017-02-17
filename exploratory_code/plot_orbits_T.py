import matplotlib.pyplot as plt
import numpy as np
from soda import *

def galactocentic(xyzMW, xyzSat):
    xyz_sat_gal = np.array([xyzSat[:,0]-xyzMW[:,0], \
                            xyzSat[:,1]-xyzMW[:,1],\
                            xyzSat[:,2]-xyzMW[:,2]]).T
    return xyz_sat_gal




font = {'size':18, 'family':'serif'}
plt.matplotlib.rc('font', **font)


data = np.loadtxt('test.txt')

time = data[:,0]
posNGCx = data[:,1]
posNGCy = data[:,2]
posNGCz = data[:,3]

velNGCx = data[:,4]
velNGCy = data[:,5]
velNGCz = data[:,6]

posSagx = data[:,7]
posSagy = data[:,8]
posSagz = data[:,9]

velSagx = data[:,10]
velSagy = data[:,11]
velSagz = data[:,12]


#plt.subplot(1, 2, 2)
#plt.plot(time,np.sqrt(relNGCSAG[:,0]**2+relNGCSAG[:,1]**2+relNGCSAG[:,2]**2),\
#         c='green', label='LMC{}'.format(str(i+1)), lw=(6-i)/2.)
#plt.legend(fontsize=13)
#plt.xlim(-5.4,0)
#plt.ylabel('$R_{gal} (Kpc)$')
#plt.xlabel('$Time (Gyrs)$')
#plt.title('$Massari+16$')

#plt.subplot(1, 2, 1)
plt.plot(time, (posSagx**2.0 + posSagy**2.0 + \
         posSagz**2.0)**0.5, label='$Sgr(+LMC)$', \
         c='b', lw=1.5)

plt.plot(time, (posNGCx**2.0 + posNGCy**2.0 + \
         posNGCz**2.0)**0.5, label='$NGC2419(+Sgr+LMC)$',\
         c='darkorange', lw=1.5)
#plt.plot(time, (posLMC[:,0]**2.0 + posLMC[:,1]**2.0 +\
#         posLMC[:,2]**2.0)**0.5, alpha=0.5, c='k',\
#         label='$LMC$', lw=(6-i)/2.)

plt.ylim(0, 160)
plt.legend(loc='best', ncol=2, fontsize=14)
plt.ylabel('$R_{gal} (Kpc)$')
plt.xlabel('$Time (Gyrs)$')
plt.show()
