#!/bin/env python

import numpy as np
import pylab as py
from waveequation import waveequation as solver

py.ion()

N, M = 64, 64
data = np.zeros((N,M,2))

for i, x in enumerate(np.linspace(-1., 1., N)):
	for j, y in enumerate(np.linspace(-1., 1., M)):
		data[i,j,1] = np.exp(-(x*x+y*y) / (0.01))

fig = py.figure(1)
py.title("t = %0.2f"%0.0)
py.subplot(121)
im = py.imshow(data[:,:,1])
py.colorbar()
py.subplot(122)
pl, = py.plot(np.linspace(-1, 1, N), data[N/2, :, 1])

nstep = 0
t = 0
while True:
  data = solver.step(data,2,1,0.5)
  t += 0.5
  if nstep%2 == 0:
    py.title("t = %0.2f"%t)
    im.set_data(data[:,:,1])
    im.autoscale()
    pl.set_ydata(data[N/2, :, 1].copy())
    py.ylim([pl.get_ydata().min(), pl.get_ydata().max()])
    py.draw()
  nstep+=1
