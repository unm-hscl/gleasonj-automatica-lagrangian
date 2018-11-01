function polyInner = minkSumInner(polyA, polyB, varargin)    
    % Computes the scaling and shifting for the template polytope polyC
    % such that it is completely contained in the polytope polyA + polyB
    % If no template is provided, we use polyA as the template
    
    if nargin > 3
        throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
    elseif nargin == 3
        polyC = varargin{1};
    else
        polyC = polyA;
    end
    
    % To simplify the problem, we will not optimize for the shift but
    % simply use the minkowski sum of Chebyshev centers (deepest points)
    % TODO: Can we guarantee that Minkowski sum's chebyshev center is this
    % point exactly?
    xc_guess = polyA.chebyCenter().x + polyB.chebyCenter().x;
    
    % Bisect
    t_ub = 10;
    t_lb = 0;
    bisect_tol = 1e-3;
    if outer_opt(t_ub, xc_guess, polyA, polyB, polyC) >= 0
        warning('SReachTools:runtime',['10x magnification of the template',...
            ' is a subset of the Minkowski sum. Please use a larger ',...
            'polytope for better results.']);
        t_opt = t_ub;
    else
        while abs(t_ub - t_lb) > bisect_tol
            t_try = (t_ub + t_lb)/2;
            optval = outer_opt(t_try, xc_guess, polyA, polyB, polyC);
            if optval >= 0
                t_lb = t_try;
                t_opt = t_try;
                fprintf('t=%1.3f:   Feasible\n',t_try);
            else
                t_ub = t_try;
                fprintf('t=%1.3f: Infeasible\n',t_try);
            end
        end
    end    
    xc_opt = xc_guess;
    polyInner = t_opt * polyC + xc_opt;
end

function optval = outer_opt(t, xc, polyA, polyB, polyC)
    % Compute the support function vector that can violate containment,
    % i.e. inf_l sup_{lambda, y, z} containment_obj
    % We need this value to be positive for guranteeing containment
    options = optimoptions('fminunc','Display','off','OptimalityTolerance',1e-8);
    [~,optval] = fminunc(@(ell) inner_opt(ell, t, xc, polyA, polyB, polyC), rand(polyA.Dim,1),options);
end

function optval = inner_opt(ell, t, xc, polyA, polyB, polyC)
    % Solve the inner optimization problem that ensures containment for a
    % given support function vector ell, scaling t, shift xc and polytopes
    % polyA, polyB (Minkowski summands) and polyC (template polytope that
    % is to be scaled and shifted)    
    cvx_begin quiet
        variable lambda(size(polyC.H,1),1);
        variable y(polyA.Dim,1);
        variable z(polyB.Dim,1);
        
        maximize (ell' * (y+z-xc) - polyC.b'*lambda)
        subject to
            polyA.A * y <= polyA.b;
            polyB.A * z <= polyB.b;
            lambda >= 0;
            polyC.A' * lambda == t * ell;
    cvx_end
    optval = cvx_optval;
end


%     options = optimoptions('fmincon','Display','iter');
%     yopt = fmincon(@(y) y(1), [0.1;xc_guess],[-1,zeros(1, length(xc_guess))],0,[],[],[],[], @(y) mycon(y, polyA, polyB, polyC),options);
%     t_opt = yopt(1);
%     xc_opt = yopt(2:end);
%
%
% function [c,ceq] = mycon(y, polyA, polyB, polyC)
%      c = -outer_opt(y(1), y(2:end), polyA, polyB, polyC);
%      ceq = [];
% end
