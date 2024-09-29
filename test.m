%% Test script for Problem 1 and Problem 2

close all; clear; clc;

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script - Problem 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example polynomial function given in problem 1 of Homework 1 
% document.
%
% Arguments:
%  x:  Polynomial independent variable
% Returns:
%  example_f_out:  Function evaluated at x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example_f_out = example_f(x)
    example_f_out = 512*x^10 - 5120*x^9 + 21760*x^8 - 51200*x^7 + ...
    72800*x^6 - 64064*x^5 + 34320*x^4 - 10560*x^3 + 1650*x^2 - 100*x + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of example polynomial function given in problem 1 of
% Homework 1 document.
%
% Arguments:
%  x:  Polynomial independent variable
% Returns:
%  example_dfdx_out:  Derivative evaluated at x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example_dfdx_out = example_dfdx(x)
    example_dfdx_out = 20*(-5 + 165*x - 1584*x^2 + 6864*x^3 - 16016*x^4 ...
    + 21840*x^5 - 17920*x^6 + 8704*x^7 - 2304*x^8 + 256*x^9);
end

% Root finding
roots = ones([1,10]);

BS_tol = 1.0e-2;
NM_tol = 1.0e-12;

roots(1) = hybrid(@example_f, @example_dfdx, 0.0, 0.04, BS_tol, NM_tol);
roots(2) = hybrid(@example_f, @example_dfdx, 0.05, 0.15, BS_tol, NM_tol);
roots(3) = hybrid(@example_f, @example_dfdx, 0.23, 0.35, BS_tol, NM_tol);
roots(4) = hybrid(@example_f, @example_dfdx, 0.47, 0.6, BS_tol, NM_tol);
roots(5) = hybrid(@example_f, @example_dfdx, 0.77, 0.9, BS_tol, NM_tol);
roots(6) = hybrid(@example_f, @example_dfdx, 1.11, 1.22, BS_tol, NM_tol);
roots(7) = hybrid(@example_f, @example_dfdx, 1.65, 1.75, BS_tol, NM_tol);
roots(8) = hybrid(@example_f, @example_dfdx, 1.86, 1.90, BS_tol, NM_tol);
roots(9) = hybrid(@example_f, @example_dfdx, 1.4, 1.5, BS_tol, NM_tol);
roots(10) = hybrid(@example_f, @example_dfdx, 1.98, 2.0, BS_tol, NM_tol);

% Plotting
xvec = linspace(0,2,10000);

fig = figure;
plot(xvec, arrayfun(@example_f, xvec), 'LineWidth', 1);
for i = 1:length(roots)
    xline(roots(i), 'LineWidth', 1);
end
yline(0, 'LineWidth', 1);
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script - Problem 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example nonlinear system given in problem 2 of Homework 1 
% document.
%
% Arguments:
%  x:  Vector of length 3. x, y, z independent variables in the
%      system.
% Returns:
%  example_sys_out:  Column ector of length 3. f1, f2, f3 
%      outputs of each function in the system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example_sys_out = example_sys(x)
    example_sys_out = zeros(3,1);
    example_sys_out(1) = x(1)^2 + x(2)^4 + x(3)^6 - 2;
    example_sys_out(2) = cos(x(1)*x(2)*x(3)^2) - x(1) - x(2) - x(3);
    example_sys_out(3) = x(2)^2 + x(3)^3 - (x(1) + x(2) - x(3))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian matrix of example nonlinear system given in problem 2 
% of Homework 1 document.
%
% Arguments:
%  x:  Vector of length 3. x, y, z independent variables in the
%      System.
% Returns:
%  example_jac_out:  3x3 matrix. Entries of the Jacobian matrix
%      for f1(x,y,z), f2(x,y,z), f3(x,y,z).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function example_jac_out = example_jac(x)
    example_jac_out(1,1) = 2*x(1);
    example_jac_out(1,2) = 4*x(2)^3;
    example_jac_out(1,3) = 6*x(3)^5;
    example_jac_out(2,1) = -x(2)*x(3)^2*sin(x(1)*x(2)*x(3)^2) - 1;
    example_jac_out(2,2) = -x(1)*x(3)^2*sin(x(1)*x(2)*x(3)^2) - 1;
    example_jac_out(2,3) = -2*x(1)*x(2)*x(3)*sin(x(1)*x(2)*x(3)^2) - 1;
    example_jac_out(3,1) = -2*x(1)-2*x(2)+2*x(3);
    example_jac_out(3,2) = -2*x(1)+2*x(3);
    example_jac_out(3,3) = 2*x(1)+2*x(2)+3*x(3)^2-2*x(3);
end

% Root finding
initial_guess = [-1.0; 0.75; 1.50];
NM_3D_tol = 1.0e-6;

roots = newtond(@example_sys, @example_jac, initial_guess, NM_3D_tol)

disp(example_sys(roots));