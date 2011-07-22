//
// Micha Gorelick, 2011
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "Grid.cpp"
#include <cmath>

#define PI 3.14159

static const double CFL  = 0.8;
static const int MAXSTEP = 20;
static const double Cs   = 1.0;

void gaussian_pulse(Grid *grid, double std, double A, int ic, int jc, int kc);
void rhs(Grid *grid, Grid *dU);
double step(Grid *grid, double CFL);

// Input:
//   grid         -> grid to initialize with gaussian pulse
//   std          -> standard deviation
//   Amp          -> Amplitude
//   ic / jc / kc -> center of pulse
void gaussian_pulse(Grid *grid, double std=.01, double Amp=1.0, int ic=-1, int jc=-1, int kc=-1) 
{
  double x,y,z,r2;
  for (int i=0; i<grid->Nx; i++) {
    x = grid->xdom[0] + i * grid->dx;

    for (int j=0; j<grid->Ny; j++) {
      y = grid->ydom[0] + j * grid->dy;

      for (int k=0; k<grid->Nz; k++) {
        z = grid->zdom[0] + k * grid->dz;
        r2 = x*x + y*y + z*z;

        grid->data[i][j][k].phi = 0.0;
        grid->data[i][j][k].A   = Amp * exp(- r2 / (2.0*std)) / sqrt(2*PI);

      }
    }
  }
}

void rhs(Grid *grid, Grid *dU)
{
  for (int i=0; i<grid->Nx; i++) {
    for (int j=0; j<grid->Ny; j++) {
      for (int k=0; k<grid->Nz; k++) {
        dU->data[i][j][k].phi 
                 = (((*grid)(i+1,j,k).A + (*grid)(i-1,j,k).A - 2 * (*grid)(i,j,k).A) / (grid->dx * grid->dx) +
                    ((*grid)(i,j+1,k).A + (*grid)(i,j-1,k).A - 2 * (*grid)(i,j,k).A) / (grid->dy * grid->dy) +
                    ((*grid)(i,j,k+1).A + (*grid)(i,j,k-1).A - 2 * (*grid)(i,j,k).A) / (grid->dz * grid->dz)
                   ) / Cs;
        dU->data[i][j][k].A
                 = (*grid)(i,j,k).phi;
      }
    }
  }
}

double step(Grid *grid, double CFL=0.8)
{
  double dt = CFL / Cs;

  Grid *dU = new Grid();
  rhs(grid, dU);
  for (int i=0; i<grid->Nx; i++) {
    for (int j=0; j<grid->Ny; j++) {
      for (int k=0; k<grid->Nz; k++) {
        grid->data[i][j][k].A   += dU->data[i][j][k].A   * dt;
        grid->data[i][j][k].phi += dU->data[i][j][k].phi * dt;
      }
    }
  }
  free(dU);

  return dt;
}

int main(int argc, char **argv)
{
  Grid *grid = new Grid(); 
  gaussian_pulse(grid);

  double t = 0;
  char filename[256];
  for (int nstep=0; nstep < MAXSTEP; nstep++)
  {
    double dt = step(grid, CFL);
    t += dt;
    printf("t = %8.4f\r",t);
    sprintf(filename, "output-%08d.txt", nstep);
    grid->dump(filename);
  }
  printf("\n");

  free(grid);
  return 0;
}
