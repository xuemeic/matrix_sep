function [x, num_iter] = my_lasso(A, b, lam, para)
% returns argmin_x {0.5|Ax-b|^2 + lam*|x|_1}
% input:
%   A: m x p
%   b: m x k
%   lam: positive scalar
%   para:
%       .method: 'ADMM', 'ISTA', 'FISTA'
%       .tol: stopping criteria
%       .max_iter
%       .V_A
%       .S2_A
%       .L_A
%       .rho: step size if method is ADMM
%       .L: lipschitz constant if method is FISTA
%       .decomp: 'svd' or 'chol'
%       .isCirculant: whether A is circulant. only used when gen_matrix_sep
%       is called.
% output:
%   .x: (p x k) 
%   .num_iter: number of iterations ran

%%%%%%%% notes %%%%%%%%%%%
% for para.method: choose ADMM or FISTA. ISTA is slow in general.

% reference for FISTA: Beck, Amir, and Marc Teboulle. 
% "A fast iterative shrinkage-thresholding algorithm for linear inverse problems." 
% SIAM journal on imaging sciences 2.1 (2009): 183-202.

% written by Xuemei Chen 8/4/2025

[~, p] = size(A);
[~, k] = size(b);
%% pass parameter
% default
max_iter_default  = 100;
tol_default       = 1e-5;
method_default    = 'FISTA';
decomp_default = 'svd';
% for ADMM, can prefactor A'A using svd or 
% prefactor A'A + rho using cholesky
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


max_iter = para.max_iter;
tol = para.tol;
method = para.method;


if strcmp(method, 'ADMM')
    if ~isfield(para, 'rho')
        para.rho = rho_default;
    end
    
    if ~isfield(para, 'decomp')
        para.decomp = decomp_default;
    end
    if ~isfield(para, 'isCirculant')
        % if not provided, evaluate
        para.isCirculant = is_circulant(A);
    end
    rho = para.rho;
    decomp = para.decomp;
    if para.isCirculant && (~isfield(para, 'coef'))
        % if A is circulant and para.coef not provided, evaluate
        d = abs(fft(A(:,1)));
        para.coef = d.^2;
    end

    if para.isCirculant
        coef = para.coef;
    else % not circulant
        if strcmp(decomp, 'svd') && (~isfield(para, 'V_A'))
        % if use svd and V_A not provided, evaluate
        [para.V_A, para.S2_A] = pref(A);
        end
        if strcmp(decomp, 'svd')
            % if use svd
            V = para.V_A;
            sig = para.S2_A;
        end
        if strcmp(decomp, 'chol') && (~isfield(para, 'L_A'))
            % if use chol and L_A not provided, evaluate
            para.L_A = chol(A'*A + rho*eye(p), 'lower');
        end
        if strcmp(decomp, 'chol')
            % if use chol
            L_A = para.L_A;            
        end

    end
    
elseif strcmp(method, 'FISTA')
    if ~isfield(para, 'L')
        % if para.L not provided, evaluate
        s = svd(A'*A);
        s = diag(s);
        para.L = s(1);        
    end
    L = para.L;
end


%% main part
x = zeros(p, k);

switch method
    case 'ADMM'
        Atb = A'*b; % m by k, Atb precomputed        
        z = zeros(p, k);
        u = zeros(p, k);
        for j = 1:max_iter              
            xlast = x;

            % update x
            rhs = Atb + rho*(z - u); 
            if para.isCirculant                 
              % [m,n]./[m,1] will divide each col
                x = ifft(fft(rhs)./(coef + rho));
                x = real(x);
            elseif strcmp(decomp, 'svd')
                x = (V'*rhs)./(sig + rho);
                x = V*x;
            elseif strcmp(decomp, 'chol')
                y = L_A \ rhs;  
                x = L_A' \ y; 
            end

            % update z
            z = SoftThresh(x + u, lam/rho);

            % update u
            u = u + x - z;

            % check convergence
            x_change = norm(x - xlast, 'fro')/(norm(xlast, 'fro') + 1);
            if x_change < tol
                break
            end
              
        end
        num_iter = j;
    case 'ISTA'
        
        s = svd(A'*A);
        s = diag(s);
        s = s(1);
        t = 1/s;
        for j = 1:max_iter
            xlast = x;
            y = x - t*A'*(A*x - b);
            x = SoftThresh(y, lam*t);
            x_change = norm(x - xlast, 'fro')/(norm(xlast, 'fro') + 1);
            if x_change < tol   
                break
            end
        end
        num_iter = j;
    case 'FISTA'
        
        
        t = 1;
        y = x;
        for j = 1:max_iter
            xlast = x;    
            tlast = t;
            x = SoftThresh(y - (1/L)*A'*(A*y - b), lam/L);
            t = (1+sqrt(1+4*t^2))/2;
            y = x + (tlast - 1)/t*(x - xlast);
            x_change = norm(x - xlast, 'fro')/(norm(xlast, 'fro') + 1);
            if x_change < tol   
                break
            end
        end
        num_iter = j;
        

end

end

%% helper functions
function [V, sig] = pref(A)
% [~, sig, V] = svd(A'*A)

[p, m] = size(A);

if m < p
   [~, sig, V] = svd(A'*A); % V is m by m       
   sig = diag(sig); % m by 1
            
else
            [~, s, V] = svd(A); % V is m by m 
            s = diag(s); % p by 1
            sig = zeros(m, 1);
            sig(1:length(s)) = (abs(s)).^2;
end
end

function flag = is_circulant(A)
    [m, n] = size(A);
    
    % we need A to be square if circulant
    if m ~= n
        flag = false;
        return;
    end

    first_row = A(1, :);
    
    % check if each row is a shift of the first row
    for i = 2:m
        expected_row = circshift(first_row, [0, i-1]);
        if any(A(i, :) ~= expected_row)
            flag = false;
            return;
        end
    end
    flag = true;
end