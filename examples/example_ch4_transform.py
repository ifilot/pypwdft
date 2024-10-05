# -*- coding: utf-8 -*-

from pypwdft import SystemBuilder, PyPWDFT, PeriodicSystem
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.image import AxesImage
import pickle
import os

def main():
    sz = 10
    npts = 32
    s = SystemBuilder().from_name('CH4', sz, npts)
    
    if os.path.exists('ch4.pickle'):
        with open('ch4.pickle', 'rb') as f:
            res = pickle.load(f)
    else:
        res = PyPWDFT(s).scf(verbose=True)
        with open('ch4.pickle', 'wb') as f:
            pickle.dump(res, f)
       
    fig, ax = plt.subplots(3, 5, dpi=144, figsize=(16,8))
    im = np.zeros((3,5), dtype=AxesImage)
    for i in range(5):
        limit = np.max(np.abs(res['orbc_rs'][i,npts//2,:,:].real))
        im[0][i] = ax[0,i].imshow(res['orbc_rs'][i,npts//2,:,:].real, extent=(0,sz,0,sz), 
                   interpolation='bicubic', cmap='PiYG',
                   vmin=-limit, vmax=limit)
        limit = np.max(np.abs(res['orbc_rs'][i,npts//2,:,:].imag))
        im[1][i] = ax[1,i].imshow(res['orbc_rs'][i,npts//2,:,:].imag, extent=(0,sz,0,sz), 
                   interpolation='bicubic', cmap='PiYG',
                   vmin=-limit, vmax=limit)
        im[2][i] = ax[2,i].imshow((res['orbc_rs'][i,npts//2,:,:] * res['orbc_rs'][i,npts//2,:,:].conjugate()).real, 
                   extent=(0,sz,0,sz), interpolation='bicubic')
    
        for j in range(0,3):
            ax[j,i].set_xlabel(r'$x$ [a.u.]')
            ax[j,i].set_ylabel(r'$y$ [a.u.]')
            
            divider = make_axes_locatable(ax[j,i])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im[j][i], cax=cax, orientation='vertical')
            
        ax[0,i].set_title(r'$\mathbb{R}\;[\psi_{%i}]$' % (i+1))
        ax[1,i].set_title(r'$\mathbb{I}\;[\psi_{%i}]$' % (i+1))
        ax[2,i].set_title(r'$\rho_{%i}$' % (i+1))

    plt.tight_layout()
    
    ### TRANSFORMATION
    
    # calculate overlap matrix prior to transformation
    S = calculate_overlap_matrix(res['orbc_rs'], sz, npts)
    print('S = ', S)
    
    # calculate kinetic energies prior to transformation
    print('Kinetic energies:')
    for i in range(5):
        print(calculate_kinetic_energy(res['orbc_fft'][i], sz, npts).real)
    
    # perform transformation
    print('\nPerforming Transformation\n')
    for i in range(5):
        res['orbc_rs'][i] = np.sign(res['orbc_rs'][i].real)*np.abs(res['orbc_rs'][i])
        
    # calculate overlap matrix after transformation
    S = calculate_overlap_matrix(res['orbc_rs'], sz, npts)
    print('S = ', S)
    
    # calculate kinetic energies after to transformation
    print('Kinetic energies:')
    Ct = np.sqrt(sz**3) / npts**3
    for i in range(5):
        print(calculate_kinetic_energy(np.fft.fftn(res['orbc_rs'][i]) * Ct, sz, npts).real)
    
    # reproduce plots
    fig, ax = plt.subplots(3, 5, dpi=144, figsize=(16,8))
    im = np.zeros((3,5), dtype=AxesImage)
    for i in range(5):
        limit = np.max(np.abs(res['orbc_rs'][i,npts//2,:,:].real))
        im[0][i] = ax[0,i].imshow(res['orbc_rs'][i,npts//2,:,:].real, extent=(0,sz,0,sz), 
                   interpolation='bicubic', cmap='PiYG',
                   vmin=-limit, vmax=limit)
        limit = 1e-16
        im[1][i] = ax[1,i].imshow(res['orbc_rs'][i,npts//2,:,:].imag, extent=(0,sz,0,sz), 
                   interpolation='bicubic', cmap='PiYG',
                   vmin=-limit, vmax=limit)
        im[2][i] = ax[2,i].imshow((res['orbc_rs'][i,npts//2,:,:] * res['orbc_rs'][i,npts//2,:,:].conjugate()).real, 
                   extent=(0,sz,0,sz), interpolation='bicubic')
    
        for j in range(0,3):
            ax[j,i].set_xlabel(r'$x$ [a.u.]')
            ax[j,i].set_ylabel(r'$y$ [a.u.]')
            
            divider = make_axes_locatable(ax[j,i])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im[j][i], cax=cax, orientation='vertical')
            
        ax[0,i].set_title(r'$\mathbb{R}\;[\psi_{%i}]$' % (i+1))
        ax[1,i].set_title(r'$\mathbb{I}\;[\psi_{%i}]$' % (i+1))
        ax[2,i].set_title(r'$\rho_{%i}$' % (i+1))

    plt.tight_layout()
    
def calculate_overlap_matrix(orbc, sz, npts):
    """
    Calculate the overlap matrix in real-space
    """
    N = len(orbc)
    S = np.zeros((N,N))
    dV = (sz / npts)**3
    for i in range(N):
        for j in range(N):
            S[i,j] = (np.sum(orbc[i].conjugate() * orbc[j]) * dV).real
            
    return S

def calculate_kinetic_energy(orbc_fft, sz, npts):
    """
    Calculate the kinetic energy of a molecular orbital as represented by a
    set of plane-wave coefficients
    """
    s = PeriodicSystem(sz=sz, npts=npts)
    
    return 0.5 * np.einsum('ijk,ijk,ijk', orbc_fft.conjugate(), s.get_pw_k2(), orbc_fft)

if __name__ == '__main__':
    main()