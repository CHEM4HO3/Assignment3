% 4H03 A3 Q1 Starter
% Jake Nease
% McMaster University

% ----------------------------------------------------------------------- %
clear variables;
close all;
clc;
% ----------------------------------------------------------------------- %

% Import the Silicon Wafer Data

% import data
raw = importdata('silicon-wafer-thickness.csv');
data = raw.data;
txt = {'G1','G2','G3','G4','G5','G6','G7','G8','G9'};

% center and scale data
X = (data - mean(data))./std(data);

% Create the subsets for cross-validation

% NOTE that the number of subsets G has been chosen to divide evenly in the
% data set... If it was not, we would round up or down and the last subset
% would absorb the "difference."

G = 4; % <-- number of subsets to be usedsubsets
blocks = length(X)/G; % <-- How many points in each testing set
X = X(randperm(length(X)),:); % <-- Shuffle the data to randomize 

% "train" and "test" are 3D matrices (arrays).
% ROW --> the observation as we know it already
% COLUMN --> the variable as we know it already
% PANE (depth dimension) --> there is a pane for each group in G

% For THIS example: 
% train is 138 x 9 x 4
% test is 46 x 9 x 4

% Example: Call an entire training set for subgroup g as train(:,:,g)

train = zeros((G-1)*blocks,size(X,2),G);
test = zeros(blocks,size(X,2),G);

for i = 1:G
    
    S = X;
    S(blocks*(i-1)+1:blocks*i,:) = [];
    train(:,:,i) = S; % <-- This is a weird one, but MATLAB won't let you delete rows from a 3D array
    test(:,:,i) = X(blocks*(i-1)+1:blocks*i,:);

end

% OK, that should have made you the data sets automatically for training
% and testing. You can use those as you proceed.


% ----------------------------------------------------------------------- %
% YOUR CODE STARTS HERE %
% ----------------------------------------------------------------------- %











