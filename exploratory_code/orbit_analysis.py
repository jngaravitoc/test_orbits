import numpy as np
import matplotlib.pyplot as plt


def load_orbit(filename):
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

    x_rel = data[:,19]
    y_rel = data[:,20]
    z_rel = data[:,21]

    r_ngc = np.array([x_ngc, y_ngc, z_ngc]).T
    v_ngc = np.array([vx_ngc, vy_ngc, vz_ngc]).T

    r_sag = np.array([x_sag, y_sag, z_sag]).T
    v_sag = np.array([vx_sag, vy_sag, vz_sag]).T

    d_rel = np.sqrt(x_rel**2 + y_rel**2 + z_rel**2)

    return t, r_ngc,  v_ngc, r_sag, v_sag , d_rel



def angular_m(r, v):
    L = np.cross(r,v)
    return L/np.linalg.norm(L)

def angles(t, r1, v1, r2, v2):
    L_NGC = angular_m(r1, v1)
    L_sag = angular_m(r2, v2)
    orb_dot = np.zeros(len(t))
    for i in range(len(t)):
        norm_NGC = np.linalg.norm(L_NGC[i])
        norm_sag = np.linalg.norm(L_sag[i])
        orb_dot[i] = np.dot(L_NGC[i]/norm_NGC, L_sag[i]/norm_sag)
    theta =  np.arccos(orb_dot)*180/np.pi
    return np.mean(theta)


def min_dist(t, drel):
    index = np.argmin(drel)
    return t[index], drel[index]

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

    #return N_peris, t_peri, np.mean(r_peri), N_apos, t_apo, np.mean(r_apo)

    return np.mean(r_peri), np.mean(r_apo)


if __name__ == "__main__": 

    N_files = 10000

    Theta = np.zeros(N_files)
    r_rel_mins = np.zeros(N_files)
    t_min_rel = np.zeros(N_files)
    peris = np.zeros(N_files)
    apos = np.zeros(N_files)

    plt.figure(figsize=(10,5))
    for i in range(0,N_files):
       t, posNGC, velNGC, posSag, velSag, d_rel = load_orbit('../orbits/MWlLMC4Sgr_ICs{}.txt'.format(str(i)))
       Theta[i-N_files] = angles(t, posNGC, velNGC, posSag, velSag)
       t_min_rel[i-N_files] , r_rel_mins[i-N_files] = min_dist(t, d_rel)
       NGC_r_G = np.sqrt(posNGC[:,0]**2 + posNGC[:,1]**2 + posNGC[:,2]**2)
       peris[i-N_files] , apos[i-N_files] = orbit_properties(t, NGC_r_G)

    font = {'size':18, 'family':'serif'}
    plt.matplotlib.rc('font', **font)

    plt.figure(figsize=(14, 14))
    plt.subplot(2, 2, 1)
    plt.scatter(t_min_rel, r_rel_mins, c='k', s=1, alpha=0.5)
    plt.ylabel('$Relative\ Distance [kpc]$')
    plt.xlabel(r'$t[Gyrs]$')


    plt.subplot(2, 2, 2)
    h = plt.hist(Theta, color='darkorange', rwidth=0.95, normed=True)
    plt.xlabel(r'$\theta [^{\circ}]$')


    plt.subplot(2, 2, 3)
    h2 = plt.hist(peris, color='purple', rwidth=0.95, normed=True)
    plt.xlabel('$r_{peri}[kpc]$')

    plt.subplot(2, 2, 4)
    h3 = plt.hist(apos, color='k', rwidth=0.95, normed=True)
    plt.xlabel('$r_{apo}[kpc]$')
    plt.savefig('model1_analysis.pdf', bbox_inches='tight', dpi=300)
