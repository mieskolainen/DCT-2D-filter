# Simulation tools
#
# Mikael Mieskolainen, 2021

import matplotlib.pyplot as plt
import numpy as np
from skimage import exposure
from PIL import Image


def add_poisson_noise(img, MAXVAL=255, clip=True):

  # Scale intensities to [0, MAXVAL]
  new = img / np.max(img[:]) * MAXVAL

  # Create random Poisson counts with mean given by the pixel values
  for i in range(new.shape[0]):
    for j in range(new.shape[1]):
      new[i,j] = np.random.poisson(lam=new[i,j])
  
  # Clip
  if clip:
    new[new < 0] = 0 # Though should not be negative due to Poisson > 0
    new[new > MAXVAL] = MAXVAL
  
  return new


def add_gaussian_noise(img, sigma=20, MAXVAL=255, clip=True):

  # Scale intensities to [0,1]
  new = img / np.max(img[:])
  
  # Add Gaussian noise
  new += (sigma / MAXVAL)*np.random.normal(loc=0, scale=1, size=img.shape)

  # Scale
  new *= MAXVAL

  # Clip
  if clip:
    new[new < 0] = 0
    new[new > MAXVAL] = MAXVAL

  return new


def readimg(filename, M, N):
    """ Read binary file
    """
    with open(filename, 'rb') as f:
        #M, N = ...  # parse; advance file-pointer to data segment
        data  = np.fromfile(f, dtype='<f8', count=M*N)
        array = np.reshape(data, [M, N], order='C')
    return array

def writeimg(img, filename):

    with open(filename, 'wb') as f:
        img.tofile(f, sep='', format='%s')
    return True


def PSNR(ref, approx):
  """ https://en.wikipedia.org/wiki/Peak_signal-to-noise_ratio
  """
  MSE      = np.mean((ref[:] - approx[:])**2)
  log10MSE = np.log10(MSE) if MSE > 0 else -np.inf

  return 20*np.log10(np.max(ref[:])) - 10*log10MSE
