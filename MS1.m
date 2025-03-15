%% MS1.m
% Main file implementing the functions for the solution of the amplifier
% system.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of MS1. 

clear; clc;

%% Constants definitions
Is = 0.96e-12;
Vt = 0.02585202; 
Rl = 1.2e3;
Rb = 22e3;
E = 1.5;
Vcc = 15;

%% Helper anonymous functions
% To assist in keeping the code for the Jacobian matrix DRY, these
% anonymous functions were created. 

Vbe = @(Ib, e) E + e - Rb*Ib;
Vbc = @(Ic, Ib, e) Vbe(Ib,e) - Vcc +Rl*Ic;
g = @(Ib, e) exp((1/Vt)*(Vbe(Ib, e)));
h = @(Ic, Ib, e) exp((1/Vt) * (Vbc(Ic, Ib, e)));


%% Jacobian Matrix
% Partial Derivatives
dF1dIc = @(Ic, Ib, e) 1 - Is*(g(Ib, e)*(-Rl/190) - 1.3333*(h(Ic, Ib, e)*(Rl/Vt)));
dF1dIb = @(Ic, Ib, e) -Is*(g(Ib, e)*(-Rb/Vt)*(1-(Vbc(Ic, Ib, e))/(190)) ...
    + (g(Ib,e)-1)*(Rb/190) - 1.33333*(h(Ic,Ib,e)*(-Rb/Vt)));
dF2dIc = @(Ic, Ib, e) -Is*((2/15)*(h(Ic, Ib, e))*(Rl/Vt));
dF2dIb = @(Ic, Ib, e) 1 - Is*((1/252)*(g(Ib,e)-1)*(-Rb/Vt) + 2/15*(h(Ic, Ib, e))*(-Rb/Vt));

%% Iteration over e(t)
% Amplitude and frequency of e(t)
A = 750e-3;
f = 1.8e3;

% Define small signal variation e(t) and num of samples per cycle
e = @(t) A*sin(2*pi*f*t);
numCycles = 3;
numSamples = 100;
sampleTime = numCycles*(1/f);

% Range of time values to iterate over
t = linspace(0, sampleTime, numCycles*numSamples);

%% Initial guess and algorithm parameters
x0 = [10.7672e-3; 42.7272e-6]; % Initial guess: [Ic; Ib]
tol = 1e-6;        %  Tolerance for convergence
maxIter = 300;      % Maximum iterations

% Store results of iteration
IcVals = zeros(size(t)); 
IbVals = zeros(size(t)); 
VoVals = zeros(size(t)); 

%% Iterate over the time samples, find Ic and Ib at each sample, find output voltage from those values
for i = 1:length(t)
    disp("*****NEW TIME SAMPLE*****")
    currentTime = t(i);
    e_t = e(currentTime);

    f = @(x) bjtSystem(x, E, e_t, Vcc, Rl, Rb);
    J = @(x) [dF1dIc(x(1), x(2), e_t), dF1dIb(x(1), x(2), e_t);
              dF2dIc(x(1), x(2), e_t), dF2dIb(x(1), x(2), e_t)];

    x = NewtonRaphson(f, J, x0, tol, maxIter);

    % Update solution vectors
    IcVals(i) = x(1); % Store collector current
    IbVals(i) = x(2); % Store base current

    VoVals(i) = Rl*x(1);

    % Use the current solution as the initial guess for the next time step
    x0 = x;
    
end
plot(t, VoVals)






