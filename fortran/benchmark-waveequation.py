#!/bin/env python

import numpy as np
from waveequation import waveequation as solver

N, M = 512, 512
data = np.zeros((N,M,2))

solver.fill(data)

for i, x in enumerate(np.linspace(-1., 1., N)):
	for j, y in enumerate(np.linspace(-1., 1., M)):
		data[i,j,1] = np.exp(-(x*x+y*y) / (0.01))

nstep = 0
t = 0
while 100:
  data = solver.step(data,2,1,0.5)
