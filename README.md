# Introduction

Grid-point Model dynamical CORE (GMCORE) is currently on the latitude-longitude grid, but we also plan to incorporate quasi-uniform grid as backup. The numerics are working on C-grid, with general terrain following vertical coordinate.

# Status

- [ ] Parallelization using MPI:
  - [X] 1D latitudional decomposition (done)
  - [X] 2D decomposition with asynchronized MPI communication (done)
  - [ ] Optimize for X86 (~2025.10)
- [ ] Nesting at middle and low latitudes (~2021.11).
- [ ] Acceleration using GPU (~?).
- [X] Baroclinic version (done).
  - [X] Hydrostatic baroclinic version (done)
    - [X] Rossby-Haurwitz wave test (done)
    - [X] Mountain induced wave test (done)
    - [X] Steady state test (done)
    - [X] Steady state PGF test (done)
    - [X] Baroclinic wave test (done)
    - [X] Held-Suarez test (done)
    - [X] Colliding Modons (done)
  - [X] Nonhydrostatic baroclinic version (done)
    - [X] X-Z version (done)
    - [X] Quasi-2D mountain wave on reduced sphere (done)
    - [X] Circular mountain wave on reduced sphere (done)
    - [X] Internal gravity wave (done)
    - [X] Splitting supercell (done)
- [X] Advection module (done)
  - [X] Update with fully monotonicity in 3D (done)
- [ ] Incorporation with physics parameterisation (2025.05-2025.10).
- [ ] Data assimilation (~?).

# Usage

First make sure you have installed netCDF library, and set `NETCDF_ROOT` environment variable to it. Then clone the repository:

```
$ git clone https://gitee.com/dongli85/GMCORE gmcore
$ cd gmcore
$ ./pull_libs.py
```

You could build the model by following:

```
$ cd build
$ FC=mpiifort cmake ..
$ make -j8
```

There is a Python script `run_tests.py`, which will clone the testbed repository, and run several tests, but it assumes MPI to be installed or SLURM job manager is available:

```
$ ./run_tests.py -w <work_directory> --slurm -q <job_queue> -n <process_number> --ntasks-per-node <n>
```

It will take some time to run the tests. When the tests are finished, cd to `<work_directory>`, and use some visualization tools, such as Panoply, to view the results.

# Authors

- Li Dong <dongli@lasg.iap.ac.cn>
- Jianghao Li

You are welcome to join our team to develop a robust global model!
