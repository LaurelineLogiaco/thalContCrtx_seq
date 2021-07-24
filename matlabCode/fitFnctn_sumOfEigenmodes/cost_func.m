function [c, g] = cost_func(x, f, dt)

    %%%% cost function that we want to minimize
    
    % x is the vector of parameters (see EignMd_Approx for definitions)
    % f is the target function
    % dt is the time step for the target function
    
    % c is the cost
    % g is the gradient of the cost relative to the parameters x (same size as x)
    
    %%%%
    
    T = numel(f);
    % Total number of timesteps in the target function

    [f_hat, Phi] = EignMd_Approx(x, T, dt);
    % f_hat is the appromation of the function using the sum of eigenmodes
    % Phi has the unweighted timecourse of the eigenmodes (their real and
    % imaginary parts separately).
    
    err = f_hat - f;
    c = 0.5 * sum(err.^2);
    % c is prop. to the mean squared error.
    
    if nargout > 1
        % then we compute te gradient of c with respect to x: g
        
        n = round(size(x, 1) / 2);
        t = dt * (1:T)';
        g = zeros(size(x));
        g(1:n,1) = sum((err' .* t) .* ...
                        (Phi(:, 1:n) .* x(1:n,2)' + Phi(:, n+1:end) .* x(n+1:end,2)'), ...
                       1)';
        g(n+1:end,1) = sum((err' .* t) .* ...
                            (Phi(:, n+1:end) .* x(1:n,2)' - Phi(:, 1:n) .* x(n+1:end,2)'), ...
                           1)';
        g(:,2) = (err * Phi)';
    end