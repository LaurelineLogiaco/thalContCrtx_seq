function [f, Phi] = EignMd_Approx(x, T, dt)
    
    % Function estimate using a (small) sum of eigenmodes = basis functions
    % which are exponentially-modulated sine waves (or complex exponentials)
    % whose timescales are in the first column of x, and whose amplitudes
    % and phases are encoded in the second column of x.
    
    n = round(size(x, 1) / 2);
    % n is the total number of real basis functions
    
    t = dt * (1:T)';
    % vector of times (unit of membrane timescales for the model).
    % dt is the timestep
    
    a_tmp = 2 * exp(t .* x(1:n,1)');
    % exponential part; x(1:n,1) contains the inverse of the timescales
    % (-> it is the real part of the eigenvalues in a dynamical system perspective).
    % This is a matrix where different rows are different times, and
    % different columns are different exponential timescales.
    
    b_tmp = -t .* x(n+1:end,1)';
    % oscillatory part; x(n+1:end,1) is proportional to the oscillatory
    % frequency (-> it is the imaginary part of the eigenvalues in a dynamical system perspective)
    
    Phi = [a_tmp .* cos(b_tmp), a_tmp .* sin(b_tmp)];
    % the first n column contains the unweighted real part of the real part
    % of the n eigenmodes, and the following n columns contain the
    % unweighted imaginary part of the n eigenmodes.
    
    f = (Phi * x(:,2))';
    % The first n values of vector x(:,2) contains the real part of the n 
    % complex weights, and the following n values contains the imaginary
    % part of these weights (the alpha^{des}_i in the equation for \hat{y}
    % in Methods section 7 (numerics)).