import numpy as np
import matplotlib.pyplot as plt
import sys

def load_orbit(filename, lmc):
    data = np.loadtxt(filename)

    t = data[:,0]
    x_ngc = data[:,1]
    y_ngc = data[:,2]
    z_ngc = data[:,3]

    vx_ngc = data[:,4]
    vy_ngc = data[:,5]
    vz_ngc = data[:,6]


    x_sag = data[:,7]
    y_sag = data[:,8]
    z_sag = data[:,9]

    vx_sag = data[:,10]
    vy_sag = data[:,11]
    vz_sag = data[:,12]


    if lmc==1:
        x_rel = data[:,19]
        y_rel = data[:,20]
        z_rel = data[:,21]

    if lmc==0:
        x_rel = data[:,13]
        y_rel = data[:,14]
        z_rel = data[:,15]

    r_ngc = np.array([x_ngc, y_ngc, z_ngc]).T
    v_ngc = np.array([vx_ngc, vy_ngc, vz_ngc]).T

    r_sag = np.array([x_sag, y_sag, z_sag]).T
    v_sag = np.array([vx_sag, vy_sag, vz_sag]).T

    d_rel = np.sqrt(x_rel**2 + y_rel**2 + z_rel**2)

    return t, r_ngc,  v_ngc, r_sag, v_sag , d_rel



def angular_m(r, v):
    """
    Compute the angular momentum of the particle
    at the present time.
    to-do:

    Guarantee that the given position and velocity
    is the present one and not the ICs.
    """

    L = np.cross(r[0],v[0])
    return L/np.linalg.norm(L)


def angles(t, r1, v1, r2, v2):
    L_NGC = angular_m(r1, v1)
    L_sag = angular_m(r2, v2)
    orb_dot = np.zeros(len(t))
    norm_NGC = np.linalg.norm(L_NGC)
    norm_sag = np.linalg.norm(L_sag)
    orb_dot = np.dot(L_NGC/norm_NGC, L_sag/norm_sag)
    theta =  np.arccos(orb_dot)*180/np.pi
    return np.mean(theta)


def min_dist(t, drel):
    """"
    Computes the minimum relative distance and the time at where that
    happens.
    """

    index = np.argmin(drel)
    return t[index], drel[index]


def tidal_radius_sag(possag):
    r = np.sqrt(possag[:,0]**2 + possag[:,1]**2 + possag[:,2]**2)
    m = 0.22072889
    b = -0.46641561
    return r*m+b


def dif_potentials(r, r_tan):
    """
    compute the potentials of the MW and Sgr at a perpedincular radius
    to the galactocentric vector of Sgr.

    To-Do:

    Generalize to any MW potential!

    """

    mw_r = (r**2+r_tan**2)**0.5
    mw_pot = soda.profiles.pot_NFW(9.86, mw_r, 0, 0, 1E12)
    sgr_pot = soda.profiles.pot_NFW(8, r_tan, 0, 0, 1E10)
    return (np.abs(mw_pot) - np.abs(sgr_pot))

def encounters(t_r, possag, drel, time):
    clo_enc = np.where(drel<t_r)[0]
    return time[clo_enc], possag[clo_enc,0], possag[clo_enc,1], possag[clo_enc,2]


def proper_motions(posngc, x_NGC, y_NGC, z_NGC, pmw, pmn):
    """
    Finding the proper motions from the positions of Sgr
    """

    index = np.where((posngc[0][0] == x_NGC) & (posngc[0][1] == y_NGC) & (posngc[0][2] == z_NGC))

    if len(index[0])>1:
        print(len(index[0]), index[0])
        print('more than 1 proper motion found!')
        return(pmw[index[0][0]], pmn[index[0][0]])
    else:
        return (pmw[index[0]], pmn[index[0]])




def orbit_properties(t, R_sag):
    """
    Function that computes the pericenters and apocenters of orbits

    input: galactocentric distance


    output:
    ------

    N_peris: Number of percienters.
    t_peri: Time at which the percienters occure.
    r_peri: Radius of pericenters.
    N_apos: Number of apocenters.
    t_apo: Time at which the apocenters occure.
    r_apo: Radius of apocenters.

    """

    t_peri = []
    r_peri = []

    t_apo = []
    r_apo = []
    for i in range(2,len(R_sag)-2,2):
        # The condition bellow can be accomplish by two nearby points 
        # therefore doing the increase by more than 1 is better.
        # the lower the resoution in time the decrease the number.
        if ((R_sag[i]<R_sag[i+2]) & (R_sag[i]<R_sag[i-2])):
            #
            r_peri.append(R_sag[i])
            t_peri.append(t[i])


        if ((R_sag[i]>R_sag[i+2]) & (R_sag[i]>R_sag[i-2])):
            #
            r_apo.append(R_sag[i])
            t_apo.append(t[i])

    N_peris = len(t_peri)
    N_apos = len(t_apo) 
    period_peri = np.abs(t_peri[1]) - np.abs(t_peri[0])
    period_apo= np.abs(t_apo[1]) - np.abs(t_apo[0])


    return t_peri[0], r_peri[0], t_apo[0], r_apo[0], period_peri, period_apo

    #return np.mean(r_peri), np.mean(r_apo)


if __name__ == "__main__": 

    #N_files = 10000

    lmc=int(sys.argv[1])
    path = sys.argv[2]
    file_name = sys.argv[3]
    out_name = sys.argv[4]

    with open(path+file_name) as f:
       name = f.read().splitlines()

    N_files = len(name)

    Theta = np.zeros(N_files)
    r_rel_mins = np.zeros(N_files)
    t_min_rel = np.zeros(N_files)
    r_peri = np.zeros(N_files)
    t_peri = np.zeros(N_files)
    r_apo = np.zeros(N_files)
    t_apo = np.zeros(N_files)
    p_apo = np.zeros(N_files)
    p_peri = np.zeros(N_files)



    ppms = np.loadtxt('mc2.ngc2419')

    pmw = ppms[:,0]
    pmn = ppms[:,1]
    xNGC = ppms[:,2]
    yNGC = ppms[:,3]
    zNGC = ppms[:,4]


    #plt.figure(figsize=(10,5))

    with open(path+file_name) as f:
       name = f.read().splitlines()
    f.close()

    f = open(out_name,'w')
    f.write('#t_peri, r_perim t_apo, r_apo, p_peri, p_apo, t_min_rel, r_rel_min, theta \n')
    for i in range(len(name)):
       t, posNGC, velNGC, posSag, velSag, d_rel = load_orbit(path+name[i], lmc)

       t_r_sag = tidal_radius_sag(posSag)
       t_enc, x_enc, y_enc, z_enc = encounters(t_r_sag, posSag, d_rel, t)
       #_enc = (x_enc**2 + y_enc**2 + z_enc**2)**0.5
       if len(t_enc>0):
           enc = 1
       else: 
           enc = 0
       pmw_here, pmn_here = proper_motions(posNGC, xNGC, yNGC, zNGC, pmw, pmn)
       Theta[i] = angles(t, posNGC, velNGC, posSag, velSag)
       t_min_rel[i] , r_rel_mins[i] = min_dist(t, d_rel)
       NGC_r_G = np.sqrt(posNGC[:,0]**2 + posNGC[:,1]**2 + posNGC[:,2]**2)
       t_peri[i], r_peri[i], t_apo[i], r_apo[i], p_peri[i], p_apo[i] = orbit_properties(t, NGC_r_G)


       f.write(("%f %f %f %f %f %f %f %f %f %d %f %f \n")%(t_peri[i],\
                r_peri[i], t_apo[i], r_apo[i], p_peri[i], p_apo[i],\
                t_min_rel[i], r_rel_mins[i], Theta[i], enc, pmw_here,\
                pmn_here))

    f.close()
    """
    font = {'size':18, 'family':'serif'}
    plt.matplotlib.rc('font', **font)

    plt.figure(figsize=(14, 14))
    plt.subplot(2, 2, 1)
    plt.scatter(t_min_rel, r_rel_mins, c='k', s=1, alpha=0.5)
    plt.ylabel('$Relative\ Distance [kpc]$')
    plt.xlabel(r'$t[Gyrs]$')
    # tidal radius
    plt.axhline(7.36)

    plt.subplot(2, 2, 2)
    h = plt.hist(Theta, color='darkorange', rwidth=0.9, normed=True, bins=15)
    plt.xlabel(r'$\theta [^{\circ}]$')


    plt.subplot(2, 2, 3)
    h2 = plt.hist(peris, color='purple', rwidth=0.9, normed=True, bins=15)
    plt.xlabel('$r_{peri}[kpc]$')

    plt.subplot(2, 2, 4)
    h3 = plt.hist(apos, color='k', rwidth=0.9, normed=True, bins=15)
    plt.xlabel('$r_{apo}[kpc]$')
    plt.savefig('model1_analysis.pdf', bbox_inches='tight', dpi=300)
    """
