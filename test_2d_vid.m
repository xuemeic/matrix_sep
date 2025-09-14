% test on a cropped video where each frame is 48 by 48
% E1 is circulant but 2 by 2
% E2 is circulant but 3 by 3


load('data/video.mat'); % get V: 240 by 320 by 300
V = uint8(V); % values are integers from 0 to 255
V = im2double(V); % ranges from 0 to 1

[h, w, numFrames] = size(V);


% Calculate center crop indices
half_crop = 24;  
center_row = floor(h/2);
center_col = floor(w/2);

row_range = (center_row - half_crop + 1):(center_row + half_crop);
col_range = (center_col - half_crop + 1):(center_col + half_crop);

% Initialize cropped video array
V_cropped = V(row_range, col_range, :);

[m1, m2, K] = size(V_cropped); % 48 x 48 x 300

n1 = 2;
rng(3);

E1r1 = rand(1, n1); 
E1r1 = E1r1/sum(E1r1); % First row of the circulant matrix 
E1 = toeplitz([E1r1(1), fliplr(E1r1(2:end))], E1r1); 
% E1 will have row sum to be 1

n2 = 3;
E2r1 = rand(1, n2); 
E2r1 = E2r1/sum(E2r1); % First row of the circulant matrix 
E2 = toeplitz([E2r1(1), fliplr(E2r1(2:end))], E2r1); 
E2 = E2';
% E2 will have column sum to be 1

G1 = kron(eye(m1/n1), E1);
G2 = kron(eye(m2/n2), E2);

% preconditioning is applied


M0 = pagemtimes(pagemtimes(G1, V_cropped), G2');
% M0 is very blurred

M = -M0 + pagemtimes(pagemtimes(G1, ones(size(M0))), G2');
%CM = pagemtimes(pagemtimes(C1, M), C2');

clear para
para.rho_outer = 1; 
para.rho_inner = 1; 
para.max_iter = 50;
para.lasso_max_iter = 20;
para.tol_outer = 1e-7;
para.lasso_tol = 1e-5;
para.is_E = true; % 
para.E1 = E1;
para.E2 = E2;
para.lasso_method = 'ADMM';
para.preconditioned = false;
para.is_E = false; %%%%%%%%%%% currently false works
lam =  1/sqrt(max(m1*m2, K));

tic
outputa = gen_matrix_sep_2d_con(M, G1, G2, lam, para);
t = toc;

desired_print(outputa, t)


para.lasso_method = 'FISTA';
tic
outputf = gen_matrix_sep_2d_con(M, G1, G2, lam, para);
t2 = toc;

desired_print(outputf, t2)


%%
figure(1)
k = 10;
subplot(1, 4, 1)
imshow(V_cropped(:,:, k))
title('original (unknown)')
subplot(1, 4, 2)
imshow(M0(:,:, k))
title('Blurred: given M0')
subplot(1, 4, 3)
imshow(1-outputa.S(:,:, k), [])
title('deblurred sparse component')
subplot(1, 4, 4)
imshow(-outputa.L(:,:, k), [])
title('background')

%%
%print("-f1", 'figs/frame10', '-djpeg', '-r300')

function desired_print(o, t)
if o.para.preconditioned
    sp = "with preconditioning";
else
    sp = "no preconditioning";
end
fprintf("******* %s, lasso by %s ******\n", sp, o.para.lasso_method)
fprintf("Number of iterations: %g.\n", o.count_outer);
fprintf("Duration: %.3f seconds.\n", t)


end
