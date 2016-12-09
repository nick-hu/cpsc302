%% *CPSC 302: Assignment 8*
% Nicholas Hu

%% Question 3
% (a)

clear variables; clc; format;

load mandrill;
M = X;

load durer;
D = X;

colormap(gray);
subplot(1, 2, 1);
image(M);
title(sprintf('mandrill: rank %d', rank(M)));
subplot(1, 2, 2);
image(D);
title(sprintf('durer: rank %d', rank(D)));

[UM, SM, VM] = svd(M, 'econ');
[UD, SD, VD] = svd(D, 'econ');

for r = 2.^(1:6)
    figure;
    colormap(gray);
    s = 1:r;
    
    T = UM(:, s) * SM(s, s) * VM(:, s)';
    subplot(1, 2, 1);
    image(T);
    title(sprintf('mandrill: rank %d', r));
    
    T = UD(:, s) * SD(s, s) * VD(:, s)';
    subplot(1, 2, 2);
    image(T);
    title(sprintf('durer: rank %d', r));
end