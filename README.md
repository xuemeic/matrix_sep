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

run `test_1d_randn.m`. No preconditioning.

Matlab output:
```
Elapsed time is 0.575970 seconds.
Relative error of recovering L0: 7.340115e-01
Relative error of recovering S0: 1.049253e-02
Ran 500 many outer loops.
Input H was circulant? Answer: 0 
```

## 2D
### H is i.i.d N(0, 1)
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