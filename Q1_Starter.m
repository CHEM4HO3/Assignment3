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
X_unshuffled = X;
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
q2 = zeros(5,1);
var_x = sum(sum(X.*X, "omitnan"),"omitnan");

for A = 1:5
    press = 0;
    for j = 1:G
        [t_train, p, r_train] = nipalspca(train(:,:,j),A);
        t_test = test(:,:,j) * p;
        E = test(:,:,j) - t_test * p';
        press = press + sum(sum(E.*E, "omitnan"),"omitnan");
    end
%     press = sum(sse);
    [~, ~, r2] = nipalspca(X,A);
    q2(A) = 1 - (press/var_x);
    fprintf('A = %d\n',A);
    fprintf('Q2 = %f\n', q2(A))
    fprintf('r2 = %f\n', r2);
    fprintf('\n');
    if(A > 1)
      per_change = (q2(A) - q2(A-1)) / q2(A-1);
      if(per_change < 0.01)
          break;
      end
    end
end
A = A-1;

[t, p, ~, res_x] = nipalspca(X_unshuffled,A);

figure
bar(p(:,1:2))
legend(["Component 1", "Component 2"]);
title("Loadings Plot")

score_figure = scoreplot(t(:,1),t(:,2));
title("Score Plot");

t2_base = t./std(t);
T2 = sum(t2_base.*t2_base,2);

N = size(X_unshuffled,1);
t2_95_perc = (((N-1)*(N+1)*A)/(N*(N-A)))*finv(0.95, A, N-A);
t2_99_perc = (((N-1)*(N+1)*A)/(N*(N-A)))*finv(0.99, A, N-A);

figure
plot(T2, '-ok');
yline(t2_95_perc,'--g', '95%');
yline(t2_99_perc,'--r', '99%');
title("T^2");

spe = sum(res_x(:,:,A).*res_x(:,:,A),2);
m_spe = mean(spe);
v_spe = var(spe);
spe_95_perc = (v_spe/(2*m_spe))*chi2inv(0.95, (2*m_spe^2)/v_spe);
spe_99_perc = (v_spe/(2*m_spe))*chi2inv(0.99, (2*m_spe^2)/v_spe);

figure
plot(spe, '-ok');
yline(spe_95_perc,'--g', '95%');
yline(spe_99_perc,'--r', '99%');
title("SPE");

[outlier_rows] = find(T2 > t2_99_perc);
for row = outlier_rows
    X_unshuffled(row,:) = [];
end

[t, p, r2, res_x] = nipalspca(X_unshuffled,A);

t2_base = t./std(t);
T2 = sum(t2_base.*t2_base,2);

N = size(X_unshuffled,1);

t2_95_perc = (((N-1)*(N+1)*A)/(N*(N-A)))*finv(0.95, A, N-A);
t2_99_perc = (((N-1)*(N+1)*A)/(N*(N-A)))*finv(0.99, A, N-A);

figure
plot(T2, '-ok');
yline(t2_95_perc,'--g', '95%');
yline(t2_99_perc,'--r', '99%');
title("T^2 without outliers");

spe = sum(res_x(:,:,A).*res_x(:,:,A),2);
m_spe = mean(spe);
v_spe = var(spe);
spe_95_perc = (v_spe/(2*m_spe))*chi2inv(0.95, (2*m_spe^2)/v_spe);
spe_99_perc = (v_spe/(2*m_spe))*chi2inv(0.99, (2*m_spe^2)/v_spe);

figure
plot(spe, '-ok');
yline(spe_95_perc,'--g', '95%');
yline(spe_99_perc,'--r', '99%');
title("SPE without outliers");










