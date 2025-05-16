load('data/video.mat'); % get V: 240 by 320 by 300
V = uint8(V);
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

[m1, m2, K] = size(V_cropped);

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
% E2 will have column sum to be 1

G1 = kron(eye(m1/n1), E1);
G2 = kron(eye(m2/n2), E2);
[CG1, C1] = precondition(G1);
[CG2, C2] = precondition(G2);

M0 = pagemtimes(pagemtimes(G1, V_cropped), G2');
% M0 is very blurred

M = -M0 + pagemtimes(pagemtimes(G1, ones(size(M0))), G2');
CM = pagemtimes(pagemtimes(C1, M), C2');

para.rho_outer = 1;
para.rho_inner = 1;
para.N_outer = 100;
para.N_inner = 30;
para.tol_outer = 1e-8;
para.tol_inner = 1e-6;
para.is_E = false;
lam =  1/sqrt(min(m1*m2, K))*0.2;

tic
outputC = gen_matrix_sep_2d(CM, CG1, CG2, lam, para);
toc
HS_output = pagemtimes(pagemtimes(G1, outputC.S), G2');
L_output = M - HS_output;

%%
figure(1)
k = 10;
subplot(1, 4, 1)
imshow(V_cropped(:,:, k))
title('original')
subplot(1, 4, 2)
imshow(M0(:,:, k))
title('Blurred frame')
subplot(1, 4, 3)
imshow(1-outputC.S(:,:, k), [])
title('deblurred sparse')
subplot(1, 4, 4)
imshow(-L_output(:,:, k), [])
title('background')
