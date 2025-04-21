function B = SVT(A,a)
[U,S,V] = svd(A,'econ');
S2 = SoftThresh(S,a);
B = U*S2*V';
end