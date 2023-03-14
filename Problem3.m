% ----------------------------------------------------------------------- %
clear variables;
close all;
clc;
% ----------------------------------------------------------------------- %

% Import the Silicon Wafer Data

% import data
raw = importdata("peas.csv");
data = raw.data;
X = data(:,1:11);
Y = data(:, 12:17);
[t, u, w_star, c, p, R2] = pls_nipals(X,Y,3);