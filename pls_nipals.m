function [t, u, w_star, c, p, R2] = pls_nipals(X, Y, num_components)

% center X and Y
X = (X - mean(X))./std(X);
Y = (Y - mean(Y))./std(Y);
X0 = X;
Y0 = Y;

% initialize arrays
t = zeros(size(X, 1), num_components);
u = zeros(size(Y, 1), num_components);
% w_star = zeros(size(X, 2), num_components);
c = zeros(size(Y, 2), num_components);
p = zeros(size(X, 2), num_components);
R2 = zeros(1, num_components);

for i = 1:num_components
    % initialize weight vector
    u = Y(:, 2);
    
    % repeat until convergence
    while true
        w_new = ((1\(u'*u))*(u'*X))';

        w_new = w_new./norm(w_new);

        t_new = (1\(w_new'*w_new))*(X*w_new);

        c_new = ((1\(t_new'*t_new))*(t_new'*Y))';

        u_new = (1\(c_new'*c_new))*(Y*c_new);
        if isnan(norm(u_new - u)/norm(u))
            i
            1 + 1
        end
        % check for convergence
        if norm(u_new - u)/norm(u)  < 1e-6
            break;
        end
        u = u_new;

    end
    
    % calculate loading vector for X
    p_new = (1\(t_new' * t_new))*(X'*t_new);
    p_new =  p_new ./ norm(p_new);
    
    % store results
    t(:, i) = t_new;
    u(:, i) = u_new;
    c(:, i) = c_new;
    w_star = w_new*(1\(p_new'*w_new));
    p(:, i) = p_new ./ norm(p_new);

    % calculate R^2
    Yhat = t * c';
    R2(i) = 1 - sum(sum((Y - Yhat).^2)) / sum(sum(Y0.^2));
    
    % update X and Y for next component
    X = X - t_new * p_new';
    Y = Y - u_new * c_new';
end


end