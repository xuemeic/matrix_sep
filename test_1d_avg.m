% Make H

n = 30;
m = 30;

first_row = [-1 1 zeros(1, m-2)]; % First row of the circulant matrix H
H = toeplitz([first_row(1), fliplr(first_row(2:end))], first_row); 