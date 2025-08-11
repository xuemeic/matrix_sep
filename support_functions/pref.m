function [V, sig] = pref(A)
% prefactorization of A
% A is m by p
% [~, sig, V] = svd(A'*A)
% V is square p by p
% sig is p by 1


[m, p] = size(A);

if p < m % tall
   [~, sig, V] = svd(A'*A); % V is p by p       
   sig = diag(sig); % p by 1
            
else % fat, m smaller
   [~, s, V] = svd(A); % V is p by p 
   s = diag(s); % m by 1
   sig = zeros(p, 1);
   sig(1:length(s)) = (abs(s)).^2;
end
end