#!/bin/env python

import numpy as np
from waveequation import waveequation as solver

N, M = 512, 512
data = np.zeros((N,M,2))

solver.fill(data)

t = 0
while 100:
  data = solver.step(data,2,1,0.5)
