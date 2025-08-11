function [X, num_iter, para] = my_lasso2(G1, G2, B, lam, para)
% returns argmin_x {0.5|(Ax-B|^2 + lam*|x|_1}
% where A = kron(G2, G1)
% input:
%   Gi: mi x pi
%   B: m1 x m2 x K
%   lam: positive scalar
%   para:
%       .method: 'ADMM', 'FISTA'
%       .tol: stopping criteria
%       .max_iter
%       .E
%       .isCirculant
%       .V1: when ADMM and svd
%       .V2: when ADMM and svd
%       .Sig: when ADMM and svd
%       .coef: when ADMM and circulant
%       .rho: step size if method is ADMM
%       .L: lipschitz constant if method is FISTA
%       .isCirculant: whether Gi are both circulant. 
% output:
%   X: (p1 x p2 x K) 
%   num_iter: number of iterations ran

%%%%%%%%%%%%%% notes %%%%%
% if para.E = true, para.E1, para.E2 are provided
% para.method = 'FISTA' does NOT work as of 8/10/2025

% reference for FISTA: Beck, Amir, and Marc Teboulle. 
% "A fast iterative shrinkage-thresholding algorithm for linear inverse problems." 
% SIAM journal on imaging sciences 2.1 (2009): 183-202.

% written by Xuemei Chen 8/9/2025

[m1, p1] = size(G1);
[m2, p2] = size(G2);
[~, ~, K] = size(B);
%% pass parameter
% default
max_iter_default  = 100;
tol_default       = 1e-5;
method_default    = 'ADMM';
% for ADMM method, can prefactor Gi'*Gi using svd
rho_default = 1;

if ~isfield(para, 'max_iter')   
    para.max_iter = max_iter_default; 
end

if ~isfield(para,'tol')   
    para.tol = tol_default ; 
end

if ~isfield(para, 'method')
    para.method = method_default;
end
if ~isfield(para, 'E')
    % if not provided
    para.E = false;
end

max_iter = para.max_iter;
tol = para.tol;
method = para.method;

if para.E
    E1 = para.E1;
    E2 = para.E2;
    [n1, ~] = size(E1);
    [n2, ~] = size(E2);
end

%keyboard
if strcmp(method, 'ADMM')
    if ~isfield(para, 'rho')
        para.rho = rho_default;
    end
    rho = para.rho;
    if ~isfield(para, 'isCirculant')
        % if not provided, evaluate
        if para.E
            para.isCirculant = is_circulant(E1) && is_circulant(E2);
        else
            para.isCirculant = is_circulant(G1) && is_circulant(G2);
        end
    end
    isCirc = para.isCirculant;
    
    if para.E && isCirc
        if ~isfield(para, 'coef')
            % if coef not provided, evaluate
            d1 = abs(fft(E1(:,1))).^2;
            d2 = abs(fft(E2(:,1))).^2;
            Lc = kron(ones(m1/n1, m2/n2), d1 * d2');
            para.coef = Lc; 
        end
        coef = para.coef;
    elseif para.E && (~isCirc)
        if ~isfield(para, 'V1')
            % if ... not provided, evaluate
            [~, eS1, W1] = svd(E1' * E1);
            S1 = kron(eye(m1/n1), eS1);
            S1 = diag(S1); % m1 by 1
            para.V1 = kron(eye(m1/n1), W1);

            [~, eS2, W2] = svd(E2' * E2);
            S2 = kron(eye(m2/n2), eS2);
            S2 = diag(S2); % m2 by 1
            para.V2 = kron(eye(m2/n2), W2);
            para.Sig = S1 * S2'; % m1 by m2
            
        end
        
        V1 = para.V1; 
        V2 = para.V2;
        Sig = para.Sig; %%%%%%%%%%%%%%%%%%%%%
    elseif isCirc % G1, G2 circ
        if (~isfield(para, 'coef'))
            % coef not provided, need to evaluate
            D1 = abs(fft(G1(:,1))).^2;
            D2 = abs(fft(G2(:,1))).^2;
            para.coef = D1 * D2'; 
        end
        coef = para.coef;
    else % G1, G2 not circ
        if  ~isfield(para, 'V1')
            % if decomposed matrices not provided, evaluate
            [para.V1, sig1] = pref(G1);
            [para.V2, sig2] = pref(G2);
            
            para.Sig = sig1 * sig2'; % p1 by p2
            
        end
        
        V1 = para.V1;        
        V2 = para.V2;        
        Sig = para.Sig;

    end

elseif strcmp(method, 'FISTA')
    if ~isfield(para, 'L')
        % if para.L not provided, evaluate
        s1 = svd(G1'*G1);
        s1 = diag(s1);
        s2 = svd(G2'*G2);
        s2 = diag(s2);
        para.L = s1(1)*s2(1);       
    end
    L = para.L; % lip constant
end


%% main part
X = zeros(p1, p2, K);

switch method
    case 'ADMM'
        rhs0 = bilinear_framewise(B, G1', G2);  % precomputed     
        Z = zeros(p1, p2, K);
        U = zeros(p1, p2, K);
        for j = 1:max_iter              
            Xlast = X;

            % update X: solve {A'*A + rho}X = A'*b + rho*(Z - U)
            rhs = rhs0 + rho*(Z - U); 
            
            if isCirc  && (~para.E)                               
                X = real(ifft2(fft2(rhs)./(coef + rho)));
            elseif isCirc && para.E
                X = bfft2(rhs, n1, n2)./(coef + rho);
                X = real(bifft2(X, n1, n2));
            else % svd, common case
                
                temp = bilinear_framewise(rhs, V1', V2);
                temp2 = temp ./ (Sig + rho);        
                X = bilinear_framewise(temp2, V1, V2');            
            end            

            % update Z
            Z = SoftThresh(X + U, lam/rho);

            % update U
            U = U + X - Z;

            % check convergence
            X_change = norm(X - Xlast, 'fro')/(norm(Xlast, 'fro') + 1);
            if X_change < tol
                break
            end
        end
        num_iter = j;
    
    case 'FISTA'                
        t = 1;
        Y = X;
        
        for j = 1:max_iter
            Xlast = X;
            
            Ay = bilinear_framewise(Y, G1, G2');
            
            temp = bilinear_framewise(Ay - B, G1', G2);
            X = SoftThresh(Y - (1/L)*temp, lam/L);
            t = (1+sqrt(1+4*t^2))/2;
            Y = X + (t-1)/t*(X - Xlast);
            X_change = norm(X - Xlast, 'fro')/(norm(Xlast, 'fro') + 1);
            if X_change < tol   
                break
            end
        end
        num_iter = j;
        

end

end

%% helper functions

function v = bilinear_framewise(V, left, right)
    % performs the bilinear frame-wise transformation on a tensor V.
    %
    % Inputs:
    %   V     - an m1 x m2 x k tensor 
    %   left  - a p x m1 matrix to be applied on the left of each frame
    %   right - a m2 x q matrix to be applied on the right of each frame
    %
    % Output:
    %   v     - a p x q x k tensor where each frame is transformed as left * V(:,:,i) * right

    % validate input dimensions
    [m1_V, m2_V, ~] = size(V);
    [~, m1_left] = size(left);
    [m2_right, ~] = size(right);

    if m1_left ~= m1_V
        error('Number of columns of left (%d) must match first dimension of V (%d).', m1_left, m1_V);
    end
    if m2_right ~= m2_V
        error('Number of rows of right (%d) must match second dimension of V (%d).', m2_right, m2_V);
    end

    LV = pagemtimes(left, V);        % Size: [p, m2_V, k]

    v = pagemtimes(LV, right);       % Size: [p, q, k]
end

function Y = bfft2(X, n1, n2)
% X is m1 x m2 x K
% ni divides mi
% X[:,:,k] has (m1/n1) by (m2/n2) blocks and each block is X_ij
% Y has the same shape as X and Y_ij = fft2(X_ij)
Y = zeros(size(X));
[m1, m2, ~] = size(X);
for i = 1:(m1/n1)
    for j = 1:(m2/n2)
        ridx = ((i-1)*n1+1):(i*n1);
        cidx = ((j-1)*n2+1):(j*n2);
        Y(ridx,cidx,:) = fft2(X(ridx,cidx,:));
    end
end

end

function Y = bifft2(X, n1, n2)
% X is m1 x m2 x K
% ni divides mi
% X[:,:,k] has (m1/n1) by (m2/n2) blocks and each block is X_ij
% Y has the same shape as X and Y_ij = ifft2(X_ij)
Y = zeros(size(X));
[m1, m2, ~] = size(X);
for i = 1:(m1/n1)
    for j = 1:(m2/n2)
        ridx = ((i-1)*n1+1):(i*n1);
        cidx = ((j-1)*n2+1):(j*n2);
        Y(ridx,cidx,:) = ifft2(X(ridx,cidx,:));
    end
end

end
