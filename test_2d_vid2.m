% test on a cropped video where each frame is 48 by 48
% E1 is circulant 16 by 16
% E2 is circulant 24 by 24
% This experiment runs about 30 seconds
% para.is_E = true or false does not make a difference


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

n1 = 16;
rng(3);

E1r1 = rand(1, n1); 
b = sort(E1r1);
n1_h = floor(n1/2);
newidx = reshape(reshape([linspace(1, n1_h, n1_h), linspace(n1, n1_h+1, n1_h)], [n1_h, 2])', [n1,1]);
E1r1(newidx) = b;
E1r1 = E1r1/sum(E1r1); % First row of the circulant matrix 
E1 = toeplitz([E1r1(1), fliplr(E1r1(2:end))], E1r1); 
% E1 will have row sum to be 1

n2 = 24;
E2r1 = rand(1, n2); 
b = sort(E2r1);
n2_h = floor(n2/2);
newidx = reshape(reshape([linspace(1, n2_h, n2_h), linspace(n2, n2_h+1, n2_h)], [n2_h, 2])', [n2,1]);
E2r1(newidx) = b;

E2r1 = E2r1/sum(E2r1); % First row of the circulant matrix 
E2 = toeplitz([E2r1(1), fliplr(E2r1(2:end))], E2r1); 
E2 = E2';
% E2 will have column sum to be 1

G1 = kron(eye(m1/n1), E1);
G2 = kron(eye(m2/n2), E2);


M0 = pagemtimes(pagemtimes(G1, V_cropped), G2');
% M0 is very blurred

M = - M0 + pagemtimes(pagemtimes(G1, ones(size(M0))), G2');


para.rho_outer = 1;
para.rho_inner = 1;
para.N_outer = 100;
para.N_inner = 30;
para.tol_outer = 1e-8;
para.tol_inner = 1e-6;
para.is_E = true;
para.E1 = E1;
para.E2 = E2;
lam =  1/sqrt(min(m1*m2, K))*0.2;

% preconditioning is applied
tic
outputC = gen_matrix_sep_2d_con(M, G1, G2, lam, para);
toc


%%
figure(1)
k = 10;
subplot(1, 4, 1)
imshow(V_cropped(:,:, k))
title('original')
subplot(1, 4, 2)
imshow(M0(:,:, k))
title('Blurred')
subplot(1, 4, 3)
imshow(1-outputC.S(:,:, k), [])
title('deblurred sparse component')
subplot(1, 4, 4)
imshow(-outputC.L(:,:, k), [])
title('background')