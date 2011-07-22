#include <iostream>
#include <fstream>
#include <cmath>

struct Data
{
	double A;
	double phi;

	double &operator[](int idx)
	{
		switch (idx) {
			case 0: return (this->A);
			case 1: return (this->phi);
			default: return this->A; // TODO: exceptions!
		}
	}

};

struct Grid
{
	int Nx, Ny, Nz;
	double dx, dy, dz;
	double xdom[2], ydom[2], zdom[2];
	Data data[32][32][32];

  inline int _mod(int a, int b)
  {
    while (a < 0 || a > b) {
      a += (a<0?+1:-1)*b;
    }
    return a;
  }

  Data &operator()(int i, int j, int k)
	{
    i = _mod(i, Nx);
    j = _mod(j, Ny);
    k = _mod(k, Nz);

		return this->data[i][j][k];
	}

	Grid()
	{
		this->Nx = 32;
		this->Ny = 32;
		this->Nz = 32;

		this->xdom = {-1., 1.};
		this->ydom = {-1., 1.};
		this->zdom = {-1., 1.};

		this->dx = (this->xdom[1] - this->xdom[0]) / this->Nx ;
		this->dy = (this->ydom[1] - this->ydom[0]) / this->Ny ;
		this->dz = (this->zdom[1] - this->zdom[0]) / this->Nz ;
	}

	double min(int idx=0)
	{
		double min = this->data[0][0][0][idx];
		for(int i=0; i < this->Nx; i++) {
			for(int j=0; j < this->Ny; j++) {
				for(int k=0; k < this->Nz; k++) {
					min = ( min > this->data[i][j][k][idx] ? this->data[i][j][k][idx] : min);
				}
			}
		}
		return min;
	}

	double max(int idx=0)
	{
		double max = this->data[0][0][0][idx];
		for(int i=0; i < this->Nx; i++) {
			for(int j=0; j < this->Ny; j++) {
				for(int k=0; k < this->Nz; k++) {
					max = ( max < this->data[i][j][k][idx] ? this->data[i][j][k][idx] : max);
				}
			}
		}
		return max;
	}


	bool dump(char *filename)
	{
		FILE *out = fopen(filename, "w+"); 
		if (out == NULL) return false;

		if (fprintf(out, "%d\t%d\t%d\n", this->Nx, this->Ny, this->Nz)<0) return false;
		if (fprintf(out, "%f\t%f\t%f\n", this->dx, this->dy, this->dz)<0) return false;
		for(int i=0; i < this->Nx; i++) {
			for(int j=0; j < this->Ny; j++) {
				for(int k=0; k < this->Nz; k++) {
					if(fprintf(out, "%f\n", this->data[i][j][k][0])<0) return false;
					if(fprintf(out, "%f\n", this->data[i][j][k][1])<0) return false;
				}
			}
		}

		fclose(out);
		return true;
	}
};
