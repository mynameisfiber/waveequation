#!/bin/env pythn

import numpy as np
import pylab as py

class Grid:
  def __init__(self, filename):
    self.filename = filename;
    self.Nx = 0;
    self.Ny = 0;
    self.Nz = 0;
    self.dx = 0.0;
    self.dy = 0.0;
    self.dz = 0.0;
    self.data = np.array(())

    self.loaddata(filename);

  def castline(self, fd, type=float):
    return [type(x) for x in fd.readline().strip().split()]

  def loaddata(self, filename=None):
    if filename is None: 
      filename = self.filename
      
    fd = file(filename, "r+")
    self.Nx, self.Ny, self.Nz = self.castline(fd, int)
    self.dx, self.dy, self.dz = self.castline(fd, float)
    
    self.data = np.array(fd.readlines(), dtype=np.double).reshape((self.Nx, self.Ny, self.Nz, 2))

if __name__ == "__main__":
  data = Grid("output.txt")
  py.imshow(data.data[16,:,:,0])
  py.show()
