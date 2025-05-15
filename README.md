# Code for General Matrix Separation

## 1D
### H is the average filter
n = 50;
m = 50;
p = 50;

run `test_1d_avg.m`. No preconditioning.

$$H=\begin{bmatrix}
1&1&&&\\
&1&1&&\\
&&1&\ddots&\\
&&&\ddots&1\\
1&&&&1
\end{bmatrix}.$$

Matlab output:
```
Elapsed time is 0.148447 seconds.
Relative error of recovering L0: 1.184148e-03
Relative error of recovering S0: 7.729501e-05
Ran 113 many outer loops.
Input H was circulant? Answer: 1 
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