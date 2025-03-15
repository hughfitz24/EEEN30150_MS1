%% NewtonRaphson.m
% M-file creating the function that implements the 
% Newton-Raphson algorithm for the solution of systems of equations.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of MS1. 

function x = NewtonRaphson(f, J, x0, tol, maxIter)
    %% Inputs:
    %   f       - Function handle for equation system
    %   J       - Function handle for Jacobian matrix
    %   x0      - Initial guess
    %   tol     - Tolerance (1e-9)
    %   maxIter - Maximum iterations (20)
    % Outputs:
    %   x       - Solution vector

    % Initial guess
    x = x0;

    % Newton-Raphson iteration
    for iter = 1:maxIter
        % Evaluate function and Jacobian
        F = f(x);
        Jacobian = J(x);

        % Solve Jacobian * delta_x = -F using Gaussian elimination
        deltaX = GaussianElimination(Jacobian, -F);

        % Update solution
        % disp("--- X: --- ")
        % fprintf("%.10f\n", x)
        % disp("--- Delta X: --- ")
        % fprintf("%.10f\n", deltaX)
        x = x + deltaX;

        % Display iteration progress
        % fprintf('Iteration %d: ', iter);
        % fprintf("%.10f, %.10f\n", x(1), x(2))

        % Check for convergence
        if norm(deltaX) < tol
            fprintf('Converged after %d iterations.\n', iter);
            return;
        end
    end

    error('Newton-Raphson did not converge after %d iterations.', maxIter);
end

