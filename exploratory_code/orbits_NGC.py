import numpy as np
import soda
import sys




def extract(dct, namespace=None):
    # function that extracts variables from kwargs
    # from:
    # http://stackoverflow.com/questions/4357851/creating-or-assigning-variables-from-a-dictionary-in-python

    if not namespace: namespace = globals()
    namespace.update(dct)

def galactocentic(xyzMW, xyzSat):
    """
    Transforming to galactocentric coordinates
    """
    xyz_sat_gal = np.array([xyzSat[:,0]-xyzMW[:,0], \
                           xyzSat[:,1]-xyzMW[:,1],\
                           xyzSat[:,2]-xyzMW[:,2]]).T
    return xyz_sat_gal


def NGC_Sgr_LMC_orbit(time, Host_model, disk_params,  Sag_model, LMC_model,\
    NGC_IC_pos, NGC_IC_vel, Sag_IC_pos, Sag_IC_vel, LMC_IC_pos,\
    LMC_IC_vel):

    t5, posLMC, velLMC, posMWsLMC, velMWsLMC, posNGCsLMC, velNGCsLMC,\
    posSagLMC, velSagLMC = soda.leapfrog.integrate_sat(time, pos_host,\
    vel_host, Host_model, disk_params, bulge_params,alpha=[0, 0.3],\
    dt=0.001, pos_sat=LMC_IC_pos,\
    vel_sat=LMC_IC_vel, satellite_model=LMC_model, pos_p=NGC_IC_pos,\
    vel_p=NGC_IC_vel, pos_sat2=Sag_IC_pos,\
    vel_sat2=Sag_IC_vel, \
    satellite_model2 = Sag_model)

    posSag_GLMC = galactocentic(posMWsLMC, posSagLMC)
    posNGC_GLMC = galactocentic(posMWsLMC, posNGCsLMC)
    posLMC_GLMC = galactocentic(posMWsLMC, posLMC)
    rel_d_sag_NGCLMC = galactocentic(posSagLMC, posNGCsLMC)
    velSag_GLMC = galactocentic(velMWsLMC, velSagLMC)
    velNGC_GLMC = galactocentic(velMWsLMC, velNGCsLMC)
    velLMC_GLMC = galactocentic(velMWsLMC, velLMC)

    del(Sag_model)

    return t5, posNGC_GLMC, velNGC_GLMC, posSag_GLMC, velSag_GLMC,\
           posLMC_GLMC, velLMC_GLMC, rel_d_sag_NGCLMC


def NGC_Sgr_orbit(time, Host_model, disk_params, Sag_model, NGC_IC_pos, NGC_IC_vel, Sag_IC_pos, Sag_IC_vel):
    t2, posSag2, velSag2, posMWs2, velMWs2, posNGCs2, velNGCs2 =\
    soda.leapfrog.integrate_sat(time, pos_host, vel_host, Host_model,\
    disk_params, bulge_params, pos_p=NGC_IC_pos, vel_p=NGC_IC_vel,\
    dt=0.001,  alpha=[0, 0.3], pos_sat=Sag_IC_pos, vel_sat=Sag_IC_vel,\
    satellite_model=satellite_model_sgr2)

    posSag_G = galactocentic(posMWs2, posSag2)
    posNGC_G = galactocentic(posMWs2, posNGCs2)

    velSag_G = galactocentic(velMWs2, velSag2)
    velNGC_G = galactocentic(velMWs2, velNGCs2)

    rel_d_sag_NGC = galactocentic(posSag2, posNGCs2)

    return t2, posNGC_G, velNGC_G, posSag_G, velSag_G, rel_d_sag_NGC


def NGC_orbit(time, Host_model, disk_params, NGC_IC_pos, NGC_IC_vel):
    point_p_model = ['NFW', 0.0, 0, 0]
    t2, posNGCG, velNGCG, posMWs2, velMWs2 =\
    soda.leapfrog.integrate_sat(time, pos_host, vel_host, Host_model,\
    disk_params, bulge_params, dt=0.001, alpha=[0,0],\
    pos_sat=NGC_IC_pos, vel_sat=NGC_IC_vel, satellite_model=point_p_model)

    return t2, posNGCG, velNGCG, posMWs2, velMWs2


def output_file(filename, time, posNGC, velNGC, posSag, velSag, rel, **kwargs):
    extract(kwargs)
    if 'posLMC' in kwargs:

        f = open(filename, 'w')
        f.write('# time, posNGC, velNGC, posSag, velSag, posLMC, velLMC, relNGCSag \n')
        for i in range(len(time)):
            f.write(('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n')\
            %(time[i], posNGC[i,0], posNGC[i,1], posNGC[i,2],\
            velNGC[i,0], velNGC[i,1], velNGC[i,2], posSag[i,0],\
            posSag[i,1], posSag[i,2], velSag[i,0], velSag[i,1],\
            velSag[i,2], posLMC[i,0], posLMC[i,1], posLMC[i,2], \
            velLMC[i,0], velLMC[i,1], velLMC[i,2], rel[i,0],\
            rel[i,1], rel[i,2]))
        f.close()
    else:
        f = open(filename, 'w')
        f.write('# time, posNGC, velNGC, posSag, velSag \n')
        for i in range(len(time)):
            f.write(('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n')\
            %(time[i], posNGC[i,0], posNGC[i,1], posNGC[i,2],\
            velNGC[i,0], velNGC[i,1], velNGC[i,2], posSag[i,0],\
            posSag[i,1], posSag[i,2], velSag[i,0], velSag[i,1],\
            velSag[i,2], rel[i,0], rel[i,1], rel[i,2]))
        f.close()
if __name__ == "__main__":

    T = 4

    pos_host = [0,0,0]
    vel_host = [0,0,0]

    host_model_l = ['NFW', 1E12, 261, 9.86]
    host_model_h = ['NFW', 1.5E12, 299, 9.56]

    bulge_params = [1E10, 0.7]

    disk_params_h = [6.5E10, 3.5, 0.53]
    disk_params_l = [5.5E10, 3.5, 0.53]


    #disk_params_l = [5E10, 3.0, 0.53] Massari
    #host_model_T = ['NFW_T', 1.45E12, 295.26, 20, 0.8, 1]


    LMC_Ms = [3E10, 1E11,  2.5E11]
    LMC_as = [3.0, 12.7,  25.2]


    sgr_pos = np.array([16.1, 2.35, -6.12])
    sgr_vel = np.array([242.5, 5.6, 228.1])


    ics = np.loadtxt('mc2.ngc2419')

    pos = [[i,j,k] for i,j,k in zip(ics[:,2], ics[:,3], ics[:,4])]
    vel = [[i,j,k] for i,j,k in zip(ics[:,5], ics[:,6], ics[:,7])]

    assert(len(pos)==10000)
    assert(len(vel)==10000)

    #pos = [-87.43, -0.51, 37.31]#NGC2419
    #vel = [16.55, 48.46, -31.33]#NGC2419

    #pos_NGC2419_2 = [-79.1-8.3, -0.5, 37.4-0.014]
    #vel_NGC2419_2 = [-32.6+11.1, -177.2+240.24, -119.3+7.25]


    LMC_model =  ['hernquist', LMC_Ms[0], LMC_as[0]]

    lmc_pos = [-1, -41, -28]
    lmc_vel = [-57, -226, 221]

    satellite_model_sgr2 = ['NFW', 1E10, 44, 8]

    if sys.argv[1] == 'MW':
        print('orbtis around the MW')

    elif sys.argv[1] == 'MWlNGC':
        for i in range(len(pos)):
            time, posNGC, velNGC, posMW, velMW = \
            NGC_orbit(T, host_model_l, disk_params_h, pos[i], vel[i])

            output_file('MWlSgr_ICs{}.txt'.format(str(i)), time, posNGC,\
            velNGC, posMW, velMW, np.ones([len(time), 3]))

    elif sys.argv[1] == 'MWlSgr':
        for i in range(len(pos)):
            time, posNGC, velNGC, posSag, velSag, relNGCSAG\
            = NGC_Sgr_orbit(T, host_model_l, disk_params_h,\
            satellite_model_sgr2, pos[i], vel[i], sgr_pos, sgr_vel)

            output_file('MWlSgr_ICs{}.txt'.format(str(i)), time, posNGC,\
            velNGC, posSag, velSag, relNGCSAG)


    elif sys.argv[1] == 'MWhSgr':
        for i in range(len(pos)):
            time, posNGC, velNGC, posSag, velSag, relNGCSAG\
            = NGC_Sgr_orbit(T, host_model_h, disk_params_l,\
            satellite_model_sgr2, pos[i], vel[i], sgr_pos, sgr_vel)

            output_file('MWhSgr_ICs{}.txt'.format(str(i)), time, posNGC, velNGC, posSag, velSag, relNGCSAG)


    elif sys.argv[1] == 'MWlSgrLMC':
        for i in range(len(pos)):
            time, posNGC, velNGC, posSag, velSag, posLMC, velLMC, \
            relNGCSAG = NGC_Sgr_LMC_orbit(T, host_model_l, disk_params_h,\
            satellite_model_sgr2, LMC_model, pos[i], vel[i], sgr_pos, sgr_vel,\
            lmc_pos, lmc_vel)

            output_file('MWlLMC0Sgr_ICs{}.txt'.format(str(i)), time, posNGC, velNGC, posSag, velSag, relNGCSAG, posLMC = posLMC, velLMC = velLMC)

    elif sys.argv[1] == 'MWhSgrLMC':

        for i in range(len(pos)):
            time, posNGC, velNGC, posSag, velSag, posLMC, velLMC, \
            relNGCSAG = NGC_Sgr_LMC_orbit(T, host_model_h, disk_params_l,\
            satellite_model_sgr2, LMC_model, pos[i], vel[i], sgr_pos, sgr_vel,\
            lmc_pos, lmc_vel)

            output_file('MWhLMC0Sgr_ICs{}.txt'.format(str(i)), time, posNGC, velNGC, posSag, velSag, relNGCSAG, posLMC = posLMC, velLMC = velLMC)
