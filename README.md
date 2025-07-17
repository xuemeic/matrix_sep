# Code for General Matrix Separation

## 1D
### H is the average filter
n = 50;
m = 50;
p = 50;

run `test_1d_avg.m`. 

$$H=\begin{bmatrix}
1&1&&&\\
&1&1&&\\
&&1&\ddots&\\
&&&\ddots&1\\
1&&&&1
\end{bmatrix}.$$

Matlab output:
```
Elapsed time is 0.157968 seconds.
Relative error of recovering L0: 1.184148e-03
Relative error of recovering S0: 7.729501e-05
Ran 113 many outer loops.
Input H was circulant? Answer: 1 
****** with preconditioning ******
Elapsed time is 0.066149 seconds.
Relative error of recovering L0: 2.522293e-05
Relative error of recovering S0: 2.402050e-06
Ran 41 many outer loops.
Input H was circulant? Answer: 0 
```

### H is iid N(0, 1)
n = 50;
m = 45;
p = 48;

run `test_1d_randn.m`. 

Matlab output:
```
Elapsed time is 0.647993 seconds.
Relative error of recovering L0: 7.340115e-01
Relative error of recovering S0: 1.049253e-02
Ran 500 many outer loops.
Input H was circulant? Answer: 0 
****** with preconditioning ******
Elapsed time is 0.080081 seconds.
Relative error of recovering L0: 7.035533e-05
Relative error of recovering S0: 2.530329e-06
Ran 95 many outer loops.
Input H was circulant? Answer: 0 
```

## 2D
### Gi are i.i.d N(0, 1)
run `test_2d_randn.m`

Matlab output:
```
Elapsed time is 1.002164 seconds.
Relative error of recovering L0: 3.156142e-01
Relative error of recovering S0: 1.530180e-03
Ran 200 many outer loops.
Input H was circulant? Answer: 0 
****** with preconditioning ******
Elapsed time is 0.128357 seconds.
Relative error of recovering L0: 1.065884e-05
Relative error of recovering S0: 7.374641e-08
Ran 54 many outer loops.
Input H was circulant? Answer: 0 
```
### Gi are circulant
run `test_2d_circ.m`

Matlab output
```
Elapsed time is 1.289915 seconds.
Relative error of recovering L0: 7.721161e-02
Relative error of recovering S0: 5.211517e-04
Ran 200 many outer loops.
Input H was circulant? Answer: 1 
****** with preconditioning ******
Elapsed time is 0.100528 seconds.
Relative error of recovering L0: 3.255672e-06
Relative error of recovering S0: 6.323612e-08
Ran 26 many outer loops.
Input H was circulant? Answer: 0 
```

### Gi block structured
run `test_2d_Ei.m`

Gi = kron(eye(mi/ni), Ei); Each Ei is NOT circulant.

Matlab output
```
Elapsed time is 4.743561 seconds.
Relative error of recovering L0: 1.615380e+01
Relative error of recovering S0: 8.527813e-01
Ran 104 many outer loops.
Input H was circulant? Answer: 0 
****** with preconditioning ******
Elapsed time is 0.832217 seconds.
Relative error of recovering L0: 6.747744e-07
Relative error of recovering S0: 1.704130e-08
Ran 38 many outer loops.
Input H was circulant? Answer: 0 
```

### Gi block circulant
run `test_2d_Ei_circ.m`

Gi = kron(eye(mi/ni), Ei); Each Ei is circulant.

Matlab output
```
Elapsed time is 7.692549 seconds.
Relative error of recovering L0: 3.199401e+01
Relative error of recovering S0: 9.344267e-01
Ran 153 many outer loops.
Input H was circulant? Answer: 0 
****** with preconditioning ******
Elapsed time is 0.903337 seconds.
Relative error of recovering L0: 6.122999e-07
Relative error of recovering S0: 1.717772e-08
Ran 34 many outer loops.
Input H was circulant? Answer: 0 
```

### video background removal and deconvolution simultaneously
run `test_2d_vid`

Matlab output
```
Elapsed time is 38.244819 seconds.
```
![results](figs/frame10.jpg)

run `test_2d_vid_full`

Matlab output
```
Elapsed time is 1188.331929 seconds.
```
![full](figs/frame10_full.jpg)