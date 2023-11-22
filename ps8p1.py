from numpy import loadtxt
from numpy.fft import rfft, irfft
import numpy as np
import pylab as pl

pianodata = np.array(loadtxt('piano.txt'))
trumpetdata =  np.array(loadtxt('trumpet.txt'))
Np = pianodata.shape[0]
Nt = trumpetdata.shape[0]



pl.figure()
pl.plot(range(Np), pianodata, color = 'orange', label =' piano')
pl.plot(range(Nt), trumpetdata, color= 'blue', label = 'trumpet')
pl.title('Signals')
pl.legend()

pft = rfft(pianodata)
ptruncated= pft[0:Np//10]
tft = rfft(trumpetdata)
ttruncated= tft[0:Nt//10]

pl.figure()
pl.plot(range(Np//10),abs(ptruncated), color = 'orange', label =' piano')
pl.plot(range(Np//10),abs(ttruncated), color= 'blue', label = 'trumpet')
pl.title('First 10k components')
pl.legend()



peak_index = np.argmax(pianodata)
dominant_frequency = np.abs(pianodata[peak_index])