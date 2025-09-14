# Code for General Matrix Separation

## 1D
### H is the average filter
n = 300;
m = n - 1;
p = m;

run `test_1d_avg.m`. 

$$H=\begin{bmatrix}
-1&1&&&\\
&-1&1&&\\
&&-1&\ddots&\\
&&&\ddots&1\\
1&&&&-1
\end{bmatrix}.$$

Matlab output:
```
******* no preconditioning, lasso by ADMM ******
Number of iterations: 500.
Duration: 28.733 seconds.
Relative error of recovering L0: 1.153768e-02
Relative error of recovering S0: 2.470040e-03

******* with preconditioning, lasso by ADMM ******
Number of iterations: 28.
Duration: 2.103 seconds.
Relative error of recovering L0: 2.379997e-05
Relative error of recovering S0: 2.141106e-06

******* no preconditioning, lasso by FISTA ******
Number of iterations: 112.
Duration: 9.504 seconds.
Relative error of recovering L0: 3.850090e-06
Relative error of recovering S0: 2.840306e-07

******* with preconditioning, lasso by FISTA ******
Number of iterations: 28.
Duration: 1.606 seconds.
Relative error of recovering L0: 8.337728e-07
Relative error of recovering S0: 2.190805e-08
```

### H is iid N(0, 1)
n = 300;
m = 270;
p = 266;

run `test_1d_randn.m`. 

Matlab output:
```
******* no preconditioning, lasso by ADMM ******
Number of iterations: 500.
Duration: 34.807 seconds.
Relative error of recovering L0: 8.845180e-01
Relative error of recovering S0: 3.399294e-03

******* with preconditioning, lasso by ADMM ******
Number of iterations: 100.
Duration: 10.627 seconds.
Relative error of recovering L0: 3.270631e-04
Relative error of recovering S0: 8.156555e-07

******* no preconditioning, lasso by FISTA ******
Number of iterations: 500.
Duration: 44.451 seconds.
Relative error of recovering L0: 7.891691e+16
Relative error of recovering S0: 3.929548e+14

******* with preconditioning, lasso by FISTA ******
Number of iterations: 100.
Duration: 2.532 seconds.
Relative error of recovering L0: 1.924595e-04
Relative error of recovering S0: 4.581211e-07
```

## 2D
### Gi are i.i.d N(0, 1)
m1 = 10;
p1 = 9;
m2 = 9;
p2 = 8;
K = 50;

run `test_2d_randn.m`

Matlab output:
```
******* Gi Gaussian, with preconditioning, lasso by ADMM ******
Number of iterations: 54.
Duration: 0.097 seconds.
Relative error of recovering L0: 1.874403e-05
Relative error of recovering S0: 1.323147e-07

******* Gi Gaussian, with preconditioning, lasso by FISTA ******
Number of iterations: 54.
Duration: 0.050 seconds.
Relative error of recovering L0: 2.270315e-06
Relative error of recovering S0: 1.943040e-08
```
### Gi are circulant
run `test_2d_circ.m`

m1 = 20;
p1 = m1;
m2 = 18;
p2 = m2;
K = 50;

Matlab output
```
******* Gi circulant, with preconditioning, lasso by ADMM ******
Number of iterations: 28.
Duration: 0.171 seconds.
Relative error of recovering L0: 1.594163e-05
Relative error of recovering S0: 4.504595e-08

******* Gi circulant, with preconditioning, lasso by FISTA ******
Number of iterations: 28.
Duration: 0.060 seconds.
Relative error of recovering L0: 7.656994e-07
```

### Gi block structured
run `test_2d_Ei.m`

m1 = 20;
m2 = 30;
p1 = m1;
p2 = m2;
K = 100;

Gi = kron(eye(mi/ni), Ei); Each Ei is NOT circulant.

Matlab output
```
******* Gi block structured, with preconditioning, lasso by ADMM ******
Number of iterations: 28.
Duration: 0.608 seconds.
Relative error of recovering L0: 1.440568e-06
Relative error of recovering S0: 3.490566e-08

******* Gi block structured, with preconditioning, lasso by FISTA ******
Number of iterations: 28.
Duration: 0.233 seconds.
Relative error of recovering L0: 5.662892e-08
Relative error of recovering S0: 1.402201e-09
```

### Gi block circulant
run `test_2d_Ei_circ.m`

Gi = kron(eye(mi/ni), Ei); Each Ei is circulant.

Matlab output
```
******* Gi block circ, no preconditioning, lasso by ADMM ******
Number of iterations: 200.
Duration: 4.043 seconds.
Relative error of recovering L0: 3.547646e+00
Relative error of recovering S0: 5.303301e-01

******* Gi block circ, with preconditioning, lasso by ADMM ******
Number of iterations: 25.
Duration: 0.557 seconds.
Relative error of recovering L0: 1.281153e-06
Relative error of recovering S0: 3.440614e-08

******* Gi block circ, no preconditioning, lasso by FISTA ******
Number of iterations: 200.
Duration: 3.902 seconds.
Relative error of recovering L0: 4.189241e-02
Relative error of recovering S0: 1.085966e-02

******* Gi block circ, with preconditioning, lasso by FISTA ******
Number of iterations: 25.
Duration: 0.343 seconds.
Relative error of recovering L0: 6.864156e-08
Relative error of recovering S0: 2.055706e-09
```

### Gi block circulant, 
run `test_2d_Ei_circ2.m`

Gi = kron(eye(mi/ni), Ei); Each Ei is circulant and mi/ni is small.

Matlab output
```
******* Gi block circ, no preconditioning, lasso by ADMM ******
Number of iterations: 56.
Duration: 2.855 seconds.
Relative error of recovering L0: 1.742138e+01
Relative error of recovering S0: 9.963165e-01

******* Gi block circ, with preconditioning, lasso by ADMM ******
Number of iterations: 28.
Duration: 0.695 seconds.
Relative error of recovering L0: 8.780395e-07
Relative error of recovering S0: 3.455678e-08

******* Gi block circ, no preconditioning, lasso by FISTA ******
Number of iterations: 61.
Duration: 1.030 seconds.
Relative error of recovering L0: 1.721960e+01
Relative error of recovering S0: 9.900553e-01

******* Gi block circ, with preconditioning, lasso by FISTA ******
Number of iterations: 28.
Duration: 0.259 seconds.
Relative error of recovering L0: 1.852779e-08
Relative error of recovering S0: 1.846655e-09
```

### video background removal and deconvolution simultaneously
run `test_2d_vid`

Matlab output
```
******* with preconditioning, lasso by ADMM ******
Number of iterations: 100.
Duration: 26.717 seconds.
******* with preconditioning, lasso by FISTA ******
Number of iterations: 100.
Duration: 10.128 seconds.
```

The recovered results are similar, whether by ADMM or FISTA.

![results](figs/frame10.jpg)

run `test_2d_vid_full`

Matlab output
```
****** with preconditioning, lasso by ADMM ******
Number of iterations: 100.
Duration: 704.386 seconds.
******* with preconditioning, lasso by FISTA ******
Number of iterations: 100.
Duration: 265.895 seconds.
```
![full](figs/frame10_full.jpg)