## Model Predictive Control with Casadi

```math
J = \phi(X_{N}) + \sum_{k=0}^{N-1}(X_{k}-X_{k, ref})^{T}Q(X_{k}-X_{k, ref})+(U_{k}-U_{k, ref})^{T}R(U_{k}-U_{k, ref})
```


### Set up casadi

1. Git clone the `CasAdi` repository from https://github.com/casadi/casadi
2. Go inside the directory and `mkdir -p build && cd build`
3. Make sure that you have `CMAKE` on your system
4. Run command `cmake .. -DWITH_IPOPT=TRUE` and `make -j4 && make install`

### Set up Eigen

1. Git clone the `Eigen3` repository from https://gitlab.com/libeigen/eigen.git
2. Go inside the directory and `mkdir -p build && cd build`
3. Run command `cmake ..` and `make && make install`


## Build and Compile CMAKE

1. Go to directory either `MECANUM` or `OMNI`
2. Go inside the directory and `mkdir -p build && cd build`
3. Run command `cmake ..` and `make`
4. Go to build direction and run  `./mpc_mencaum` or `./mpc_omni`.