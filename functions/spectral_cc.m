function [xarray, resarray] = spectral_cc(A,T,varargin)

    %%% Computes x = alpha * Ax + (1-alpha) * Tp(x)

    addpath tensor_toolbox/
    n = length(A);

    par = inputParser;
    par.addParameter('x0', ones(n,1));
    par.addParameter('maxiter', 100);
    par.addParameter('gamma', .9);
    par.addParameter('alpha', .5);
    par.addParameter('p', 0);
    par.addParameter('tol', 1e-6);
    par.parse(varargin{:});
    
    x0 = par.Results.x0;
    maxiter = par.Results.maxiter;
    gamma = par.Results.gamma;
    alpha = par.Results.alpha;
    tol = par.Results.tol;
    p = par.Results.p;
    
    x = x0;    xarray = x0;

    for k  = 1 : maxiter
        v = alpha.*A*x; 
        w = Tp(T,x,p);
        w = (1-alpha).*w;
        x = (v+w); x = x./norm(x,1); x = x.^gamma;
        xarray = [xarray x];
        norm(xarray(:,end-1)-x,1)
        if norm(xarray(:,end-1)-x,1) < tol
            break;
        end
    end
    
    if k == maxiter, warning('Maxiter reached!'); end
end

