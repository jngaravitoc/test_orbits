#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np



font = {'size':18, 'family':'serif'}
plt.matplotlib.rc('font', **font)

def read_file(filename):
    data = np.loadtxt(filename)
    t = data[:,0]
    xNGC = data[:,1]
    yNGC = data[:,2]
    zNGC = data[:,3]

    xSag = data[:,7]
    ySag = data[:,8]
    zSag = data[:,9]

    ncols = np.shape(data)[1]

    if ncols==22:
        xLMC = data[:,13]
        yLMC = data[:,14]
        zLMC = data[:,15]
        x_rel = data[:,19]
        y_rel = data[:,20]
        z_rel = data[:,21]

        return t, xNGC, yNGC, zNGC, xSag, ySag, zSag, xLMC, yLMC, zLMC, x_rel, y_rel, z_rel

    if ncols==16:
        x_rel = data[:,13]
        y_rel = data[:,14]
        z_rel = data[:,15]

        return t, xNGC, yNGC, zNGC, xSag, ySag, zSag, x_rel, y_rel, z_rel

def plot_orbits(filename):

    font = {'size':18, 'family':'serif'}
    plt.matplotlib.rc('font', **font)

    data = read_file(filename)
    ncols = np.shape(data)[0]
    print(ncols)
    if ncols == 13:
        t = data[0]
        xNGC = data[1]
        yNGC = data[2]
        zNGC = data[3]
        xSag = data[4]
        ySag = data[5]
        zSag = data[6]
        xLMC = data[7]
        yLMC = data[8]
        zLMC = data[9]
        x_rel = data[10]
        y_rel = data[11]
        z_rel = data[12]

        fig = plt.figure(figsize=(14, 5))
        plt.subplot(1, 2, 2)
        plt.plot(t, np.sqrt(x_rel**2+y_rel**2+z_rel**2),\
                 c='green', label='LMC', lw=2.)
        plt.legend(fontsize=13)
        plt.xlim(-5.4,0)
        plt.ylabel('$R_{gal} (Kpc)$')
        plt.xlabel('$Time (Gyrs)$')

        plt.subplot(1, 2, 1)
        plt.plot(t, (xSag**2.0 + ySag**2.0 + \
                 zSag**2.0)**0.5, label='$Sgr(+LMC)$', \
                 c='b', lw=2.)

        plt.plot(t, (xNGC**2.0 + yNGC**2.0 + \
                 zNGC**2.0)**0.5, label='$NGC2419(+Sgr+LMC)$',\
                 c='darkorange', lw=2.)
        plt.plot(t, (xLMC**2.0 + yLMC**2.0 +\
                 zLMC**2.0)**0.5, alpha=0.5, c='k',\
                 label='$LMC$', lw=2.)

        plt.ylim(0, 160)

        plt.ylabel('$R_{gal} (Kpc)$')
        plt.xlabel('$Time (Gyrs)$')
        return fig
        #plt.savefig('gal_orbits_all_LMCs_massari.pdf', dpi=300, bbox_inches='tight')
