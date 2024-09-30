%% Problem 1 - Hybrid algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A hybrid algorithm that uses bisection and Newton's method
% to locate a root within a given interval [xmin, xmax]. If the 
% number of iterations exceeds for either bisection or Newton's
% method returns 50, the function returns NaN.
%
% Arguments:
%  f:      Function whose root is sought (takes one argument).
%  dfdx:   Derivative function (takes one argument).
%  xmin:   Initial bracket minimum.
%  xmax:   Initial bracket maximum.
%  tol1:   Relative convergence criterion for bisection.
%  tol2:   Relative convergence criterion for Newton iteration.
% Returns:
%  x:      Estimate of root.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = hybrid(f, dfdx, xmin, xmax, tol1, tol2)
    overflow_counter = 0;
    MAX_ITERATIONS = 50;

    % Bisection:
    converged = false;
    fmin = f(xmin);
    while not(converged)
        xmid = (xmin + xmax)/2;
        fmid = f(xmid);
        if fmid == 0
            break
        elseif fmid*fmin < 0
            xmax = xmid;
        else
            xmin = xmid;
            fmin = fmid;
        end
        if (xmax - xmin)/abs(xmid) < tol1
            converged = true;
        end 

        overflow_counter = overflow_counter + 1;
        if overflow_counter == MAX_ITERATIONS
            x = NaN;
            return;
        end 
    end
    bisection_result = xmid;

    overflow_counter = 0;

    % Newton's method:
    converged = false;
    x = bisection_result; 
    xprev = bisection_result;
    while not(converged)
        x = xprev - f(xprev)/dfdx(xprev);
        if abs((x - xprev)/xprev) < tol2
            converged = true;
        end
        xprev = x;

        overflow_counter = overflow_counter + 1;
        if overflow_counter == MAX_ITERATIONS
            x = NaN;
            return;
        end 
    end

end
