% 4H03 A3 Q2
% Michael Djurdjevic, djurdjm, 400132129 
% Aron Markandaier, markanda, 400121110
% Harsahib Matharoo, matharoh, 400185871
% McMaster University

clear variables;
close all;
clc;

% import data
raw = importdata('peas.csv');
data = raw.data;
txt = raw.textdata;

% centre + scale
data = (data - mean(data))./std(data);

% grab X and Y portions
X = data(:,1:11);
Y = data(:,12:17);

% X scatterplot matrix

[~,ax] = plotmatrix(X); % ignore all parameters except subplot axes
% ax in this case is a K x K (K columns being plotted) matrix of axes for each subplot
font = 12; % font size
title("Scatterplot Matrix - Pea Measurements");

% Formatting subplots - you can also make your own array of titles like title = {'one', 'two', 'three'...}
for i = 1:length(ax)
ax(i,1).YLabel.String = string(i);
ax(i,1).YLabel.FontSize = font;
ax(end,i).XLabel.String = string(i);
ax(end,i).XLabel.FontSize = font;
end

% 2 component PCA on X + loading plot

figure()

[t,p,R2] = nipalspca(X, 2);
bar(p(:,1:2))
legend(["Component 1", "Component 2"]);
title("Loadings Plot");
xticklabels(txt(1:11));

% PCR model + 3D scatterplot

figure()

% a coefficients
a = (t'*t)^(-1)*t'*Y;

% prediction
yhat = t*a;

% 3D scatter
hold on;
grid on;
scatter3(t(:,1), t(:,2), Y(:,1), 15, "filled", "green");
scatter3(t(:,1), t(:,2), yhat(:,1), 15, "filled", "red");

% Create a grid of points to evaluate the model surface
[t1grid,t2grid] = meshgrid(linspace(min(t(:,1)),max(t(:,1)),25), linspace(min(t(:,2)),max(t(:,2)),25));
Zgrid = zeros(size(t1grid));

% Evaluate the model surface at each point in the grid
for i = 1:size(t1grid,1)
    for j = 1:size(t1grid,2)
        x = [t1grid(i,j), t2grid(i,j)];
        z = x*a;
        Zgrid(i,j) = z(1);
    end
end

% Plot the model surface
surf(t1grid, t2grid, Zgrid, 'FaceAlpha', 0.5);
view(30,30);
title('PCR Model Prediction of Flavour Rating');
xlabel("t1");
ylabel("t2");
zlabel("Flavour Rating Prediction");
legend(["Real", "Predicted"])
hold off;

% Y scatterplot matrix
figure()

[~,ax] = plotmatrix(Y); % ignore all parameters except subplot axes
% ax in this case is a K x K (K columns being plotted) matrix of axes for each subplot
font = 12; % font size
title("Scatterplot Matrix - Pea Ratings");

% Formatting subplots - you can also make your own array of titles like title = {'one', 'two', 'three'...}
for i = 1:length(ax)
ax(i,1).YLabel.String = string(i);
ax(i,1).YLabel.FontSize = font;
ax(end,i).XLabel.String = string(i);
ax(end,i).XLabel.FontSize = font;
end



