function z = SoftThresh(x, kappa)

% soft thresholding operator, works for vectors or matrices
% z(i)=x(i) - kappa if x(i)>kappa
% z(i)=x(i) + kappa if x(i)<-kappa
% z(i)=0 otherwise

    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end