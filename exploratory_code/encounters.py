import soda
import numpy as np

def load_data(filename):
    data = np.loadtxt(filename)

    t = data[:,0]
    x_sag = data[:,7]
    y_sag = data[:,8]
    z_sag = data[:,9]

    x_rel = data[:,19]
    y_rel = data[:,20]
    z_rel = data[:,21]

    return t, x_sag, y_sag, z_sag, x_rel, y_rel, z_rel



"""
def tidal_radius_sag(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    t_r_100 = np.logspace(0.001, 2, 1000)
    t_r = np.zeros(len(x))
    for i in range(len(r)):
        pot_r_100 = dif_potentials(r[i], t_r_100)
        max_pot = np.argmax(pot_r_100)
        t_r[i] = t_r_100[max_pot]
    return t_r
"""

def tidal_radius_sag(x, y, z):
    """
    Approximate tidal radius of sgr at the tangential point
    of the ...
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    m = 0.22072889
    b = -0.46641561
    return r*m+ b

def dif_potentials(r, r_tan):
    mw_r = (r**2 + r_tan**2)**0.5
    mw_pot = soda.profiles.pot_NFW(9.86, mw_r, 0, 0, 1E12)
    sgr_pot = soda.profiles.pot_NFW(8, r_tan, 0, 0, 1E10)
    #print('MW pot:', np.abs(mw_pot))
    #print('Sgr pot:', np.abs(sgr_pot))
    return (np.abs(mw_pot) - np.abs(sgr_pot))


def encounters(t_r, x_sag, y_sag, z_sag, x_rel, y_rel, z_rel, time):
    clo_enc = np.where((x_rel**2 + y_rel**2 + z_rel**2)**0.5<t_r)[0]
    return time[clo_enc], x_sag[clo_enc], y_sag[clo_enc], z_sag[clo_enc]


if __name__ == "__main__":



    with open("../orbits/files.txt") as f:
        files = list(f.read().split()) 
    file_path = '../orbits/'

    f = open('encounters_saglmc.txt', 'w')
    f.write('# tenc, xenc, yenc, zenc, renc  \n')
    fraction = 0
    for name in files:
        t, x_sag, y_sag, z_sag, x_rel, y_rel, z_rel = load_data(file_path + name)
        t_r_sag = tidal_radius_sag(x_sag, y_sag, z_sag)
        t_enc, x_enc, y_enc, z_enc = encounters(t_r_sag, x_sag, y_sag, z_sag, x_rel, y_rel, z_rel, t)
        r_enc = (x_enc**2 + y_enc**2 + z_enc**2)**0.5
        print(r_enc)
        for j in range(len(r_enc)):
            f.write(('%f %f %f %f %f %s \n')%(t_enc[j], x_enc[j],\
                                              y_enc[j], z_enc[j],\
                                              r_enc[j], name))
        if len(r_enc)>0:
            fraction+=1

    f.close()
    print(fraction)
