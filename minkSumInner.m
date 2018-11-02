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
    if isEmptySet(polyA) || isEmptySet(polyB) || isEmptySet(polyC)
        throwAsCaller(SrtInvalidArgsError('Empty set provided as arguments'));
    end
    if ~((polyA.Dim == polyB.Dim) && (polyB.Dim == polyC.Dim))
        throwAsCaller(SrtInvalidArgsError('Mismatch in polytope dimensions'));
    end
    % Ensure that polyC is "centered" (This way uniform scaling is not partial)
    polyC = polyC - polyC.chebyCenter.x;    
    
    % To simplify the problem, we will not optimize for the shift but
    % simply use the minkowski sum of Chebyshev centers (deepest points)
    % TODO: Can we guarantee that Minkowski sum's chebyshev center is this
    % point exactly?
    xc_guess = polyA.chebyCenter().x + polyB.chebyCenter().x;
    
    
    % Bisect
    t_ub = 10;
    t_lb = 0;
    bisect_tol = 1e-2;
    myzero_lb = -1e-6;
    if outer_opt(t_ub, xc_guess, polyA, polyB, polyC) >= 0
        warning('SReachTools:runtime',['10x magnification of the template',...
            ' is a subset of the Minkowski sum. Please use a larger ',...
            'polytope for better results.']);
        t_opt = t_ub;
    else
        fprintf('Bisection upper bound of t=%1.3f is valid\n',t_ub);
        t_opt = 0;
        while abs(t_ub - t_lb) > bisect_tol
            t_try = (t_ub + t_lb)/2;
            optval = outer_opt(t_try, xc_guess, polyA, polyB, polyC);
            if optval >= myzero_lb
                t_opt = t_try;
                t_lb = t_try;
                fprintf('t=%1.3f:   Feasible\n',t_try);
            else
                t_ub = t_try;
                fprintf('t=%1.3f: Infeasible\n',t_try);
            end
        end
    end    
    if t_opt < bisect_tol
        warning('SReachTools:runtime', sprintf(['%1.2fx magnification of ',...
            'the template is also not a subset of the Minkowski sum. Please',...
            ' use a smaller polytope for a non-trivial results.'],bisect_tol));
    end
    xc_opt = xc_guess;
    polyInner = t_opt * polyC + xc_opt;
end

function optval = outer_opt(t, xc, polyA, polyB, polyC)
    % Compute the support function vector that can violate containment,
    % i.e. inf_l sup_{lambda, y, z} containment_obj
    % We need this value to be positive for guranteeing containment
    options = optimoptions('fmincon','Display','off');
    [~,optval] = fmincon(@(ell) inner_opt(ell, t, xc, polyA, polyB, polyC),...
        0.5*ones(polyA.Dim,1), [], [], [], [],...
        -ones(polyA.Dim,1), ones(polyA.Dim,1), [], options);
end

function optval = inner_opt(ell, t, xc, polyA, polyB, polyC)
    % Solve the inner optimization problem that ensures containment for a
    % given support function vector ell, scaling t, shift xc and polytopes
    % polyA, polyB (Minkowski summands) and polyC (template polytope that
    % is to be scaled and shifted)    
    cvx_clear
    cvx_begin quiet
        cvx_precision best
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
    if strcmpi(cvx_status,'Solved')        
    elseif strcmpi(cvx_status,'Unbounded')
        % Setting it to a large pos value so that fmincon won't freak out
        optval = 1e10;
    elseif strcmpi(cvx_status,'Infeasible')
        % Setting it to a large neg value so that fmincon won't freak out
        optval = -1e10;
    else
        keyboard
    end
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
