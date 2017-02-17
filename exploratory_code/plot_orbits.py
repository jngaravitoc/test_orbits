import matplotlib.pyplot as plt
import numpy as np
from soda import *

font = {'size':18, 'family':'serif'}
plt.matplotlib.rc('font', **font)

time = 4

# GC 6D coordinates

satellite_model = ['hernquist', 1E10, 9.8]
satellite_model_sgr2 = ['NFW', 1E10, 44, 8]

sgr_pos = [16.1, 2.35, -6.12]
sgr_vel = [242.5, 5.6, 228.1]

pos_host = [0,0,0]
vel_host = [0,0,0]

pos_NGC2419 = [-87.43, -0.51, 37.31]
vel_NGC2419 = [16.55, 48.46, -31.33]

pos_NGC2419_2 = [-79.1-8.3, -0.5, 37.4-0.014]
vel_NGC2419_2 = [-32.6+11.1, -177.2+240.24, -119.3+7.25]
host_model = ['NFW', 1E12, 261, 9.86]
disk_params = [6.5E10, 3.0, 0.53]
bulge_params = [1E10, 0.7]


def galactocentic(xyzMW, xyzSat):
    """
    Transforming to galactocentric coordinates
    """
    xyz_sat_gal = np.array([xyzSat[:,0]-xyzMW[:,0], \
                           xyzSat[:,1]-xyzMW[:,1],\
                           xyzSat[:,2]-xyzMW[:,2]]).T
    return xyz_sat_gal


def NGC_Sgr_LMC_orbit(LMC_mod, satellite_model_sgr2):
    t5, posLMC, velLMC, posMWsLMC, velMWsLMC, posNGCsLMC, velNGCsLMC,\
    posSagLMC, velSagLMC = leapfrog.integrate_sat(4, pos_host,\
    vel_host,host_model, disk_params, bulge_params, pos_sat=[-1,-41, -28],\
    vel_sat=[-57,-226, 221], satellite_model=LMC_mod, pos_p=pos_NGC2419_2,\
    vel_p=vel_NGC2419_2, dt=0.001, alpha=[0, 0.3], pos_sat2=sgr_pos,\
    vel_sat2=sgr_vel, satellite_model2=satellite_model_sgr2)

    posSag_GLMC = galactocentic(posMWsLMC, posSagLMC)
    posNGC_GLMC = galactocentic(posMWsLMC, posNGCsLMC)
    posLMC_GLMC = galactocentic(posMWsLMC, posLMC)
    rel_d_sag_NGCLMC = galactocentic(posSagLMC, posNGCsLMC)

    return t5, posNGC_GLMC, posSag_GLMC, posLMC_GLMC, rel_d_sag_NGCLMC


LMC_Ms = [3E10]#, 1E11,  2.5E11]
LMC_as = [3.0]#, 12.7,  25.2]

fig = plt.figure(figsize=(14, 5))
for i in range(len(LMC_Ms)):
    LMC_mod1 = ['hernquist', LMC_Ms[i], LMC_as[i]]
    time, posNGC, posSag, posLMC, relNGCSAG = NGC_Sgr_LMC_orbit(LMC_mod1, satellite_model_sgr2)
    plt.subplot(1, 2, 2)
    plt.plot(time,np.sqrt(relNGCSAG[:,0]**2+relNGCSAG[:,1]**2+relNGCSAG[:,2]**2),\
             c='green', label='LMC{}'.format(str(i+1)), lw=(6-i)/2.)
    plt.legend(fontsize=13)
    plt.xlim(-5.4,0)
    plt.ylabel('$R_{gal} (Kpc)$')
    plt.xlabel('$Time (Gyrs)$')
    plt.title('$Massari+16$')

    plt.subplot(1, 2, 1)
    plt.plot(time, (posSag[:,0]**2.0 + posSag[:,1]**2.0 + \
             posSag[:,2]**2.0)**0.5, label='$Sgr(+LMC)$', \
             c='b', lw=(6-i)/2.)

    plt.plot(time, (posNGC[:,0]**2.0 + posNGC[:,1]**2.0 + \
             posNGC[:,2]**2.0)**0.5, label='$NGC2419(+Sgr+LMC)$',\
             c='darkorange', lw=(6-i)/2.)
    plt.plot(time, (posLMC[:,0]**2.0 + posLMC[:,1]**2.0 +\
             posLMC[:,2]**2.0)**0.5, alpha=0.5, c='k',\
             label='$LMC$', lw=(6-i)/2.)

    plt.ylim(0, 160)
    #plt.legend(loc='best', ncol=2, fontsize=14)
    plt.ylabel('$R_{gal} (Kpc)$')
    plt.xlabel('$Time (Gyrs)$')

#plt.savefig('gal_orbits_all_LMCs_massari.pdf', dpi=300, bbox_inches='tight')
