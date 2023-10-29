#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 10:32:26 2023

@author: f003fh5
"""



from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

hdu_list = fits.open('specgrid.fits')
logwave = hdu_list['LOGWAVE'].data
flux = hdu_list['FLUX'].data

plt.plot(logwave, flux[0, :])
plt.ylabel('flux [$10^{−17}$ erg s$^{−1}$ cm$^{−2}$ A$^{-1}$]', fontsize = 16)
plt.xlabel('wavelength [$A$]', fontsize = 16)

flux_sum = np.sum(flux, axis = 1)
flux_normalized = flux/np.tile(flux_sum, (np.shape(flux)[1], 1)).T
means_normalized = np.mean(flux_normalized, axis=1)
flux_normalized_0_mean = flux_normalized-np.tile(means_normalized, (np.shape(flux)[1], 1)).T

def sorted_eigs(r, return_eigvalues = False):
    """
    Calculate the eigenvectors and eigenvalues of the correlation matrix of r
    -----------------------------------------------------
    """
    corr=r.T@r
    eigs=np.linalg.eig(corr) #calculate eigenvectors and values of original 
    arg=np.argsort(eigs[0])[::-1] #get indices for sorted eigenvalues
    eigvec=eigs[1][:,arg] #sort eigenvectors
    eig = eigs[0][arg] # sort eigenvalues
    if return_eigvalues == True:
        return eig, eigvec
    else:
        return eigvec
#Corr
    
r = flux_normalized_0_mean 
r_subset = r[:1000, :]
logwave_subset = logwave
C = r_subset.T@r_subset # correlation matrix, dimension # wavelengths x # wavelengths


#SVD

U, S, Vh = np.linalg.svd(r_subset, full_matrices=True)
# rows of Vh are eigenvectors
eigvecs_svd = Vh.T
eigvals_svd = S**2
svd_sort = np.argsort(eigvals_svd)[::-1]
eigvecs_svd = eigvecs_svd[:,svd_sort]
eigvals_svd = eigvals_svd[svd_sort]

#Eig
eigvals, eigvecs = sorted_eigs(r_subset, return_eigvalues = True)

plt.figure(2)
[plt.plot(logwave, eigvecs.T[:,i], 'o')for i in range(1000)]
plt.xlabel('Wavelengths (A)', fontsize = 16)
plt.ylabel('Eig eigenvectors', fontsize = 16)

#PCA
def PCA(l, r, project = True):
    """
    Perform PCA dimensionality reduction
    --------------------------------------------------------------------------------------
    """
    eigvector = sorted_eigs(r)
    eigvec=eigvector[:,:l] #sort eigenvectors, only keep l
    reduced_wavelength_data= np.dot(eigvec.T,r.T) #np.dot(eigvec.T, np.dot(eigvec,r.T))
    if project == False:
        return reduced_wavelength_data.T # get the reduced wavelength weights
    else: 
        return np.dot(eigvec, reduced_wavelength_data).T # multiply eigenvectors by 
                                                        # weights to get approximate spectrum

plt.figure(3)
plt.plot(logwave_subset, PCA(5,r_subset)[1,:], label = 'l = 5')
plt.plot(logwave_subset, r_subset[1,:], label = 'original data')

plt.ylabel('normalized 0-mean flux', fontsize = 16)
plt.xlabel('wavelength [$A$]', fontsize = 16)
plt.legend()                                                       


plt.figure(4)
Nc_values = range(20)[1:]
plt.figure(5)
for Nc in Nc_values:
    plt.plot(logwave, r_subset[1,:]- PCA(Nc,r_subset)[1,:] ,label=f'Nc = {Nc}')

plt.xlabel('Log wavelength (Angstrom)')
plt.ylabel('Residuals')
plt.legend()
plt.show()

def get_pc_coefficients(Nc, fl, eigenv):
    # Project centered data onto the first Nc eigenvectors
    coefficients = np.dot(fl, eigenv[:, :Nc])
    return coefficients

coefficients = get_pc_coefficients(Nc, r, eigvecs)

plt.figure(5)
plt.scatter(coefficients[:, 0], coefficients[:, 1], c='blue', label='Coefficient 1')
plt.scatter(coefficients[:, 0], coefficients[:, 2], c='red', label='Coefficient 2')
plt.xlabel('Coefficient 0')
plt.ylabel('Coefficient Value')
plt.title('Coefficient 0 vs Coefficients 1 and 2')
plt.legend()
plt.show()