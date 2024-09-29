%% Problem 2 - D-dimensional Newton iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton's method for a d-dimensional space.
%
% Arguments:
%  f:   Function which implements the nonlinear system of 
%       equations. Function is of the form f(x) where x is a 
%       length-d vector, and which returns length-d column 
%       vector.
%  jac: Function which is of the form jac(x) where x is a 
%       length-d vector, and which returns the d x d matrix of 
%       Jacobian matrix elements for the nonlinear system defined 
%       by f.
%  x0:  Initial estimate for iteration (length-d column vector).
%  tol: Convergence criterion: routine returns when relative 
%       magnitude of update from iteration to iteration is 
%       <= tol.
% Returns:
%  x:   Estimate of root (length-d column vector).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = newtond(f, jac, x0, tol)
    x = x0;
    res = f(x0);
    dx = jac(x0)\res;
    while rms(dx) > tol
        res = f(x);
        dx = jac(x)\res;
        x = x - dx;
    end 
end
