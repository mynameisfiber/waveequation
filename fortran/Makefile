waveequation.so: waveequation.f90
        f2py -c -m waveequation --opt="-O3" --f90flags="-Wall -fopenmp" -lgomp -DF2PY_REPORT_ON_ARRAY_COPY=100 waveequation.f90
