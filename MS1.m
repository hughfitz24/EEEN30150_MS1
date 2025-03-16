%% MS1.m
% Main file implementing the functions for the solution of the amplifier
% system.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of MS1. 

% Clear command window and workspace
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

% Circuit equations, see report for derivation
Vbe = @(Ib, e) E + e - Rb*Ib;
Vbc = @(Ic, Ib, e) Vbe(Ib,e) - Vcc +Rl*Ic;
% Repeated exponentials in system of equations
g = @(Ib, e) exp((1/Vt)*(Vbe(Ib, e)));
h = @(Ic, Ib, e) exp((1/Vt) * (Vbc(Ic, Ib, e)));


%% Jacobian Matrix
% Partial Derivatives
dF1dIc = @(Ic, Ib, e) 1 - Is*(g(Ib, e)*(-Rl/190) ...
    - 1.3333*(h(Ic, Ib, e)*(Rl/Vt)));
dF1dIb = @(Ic, Ib, e) -Is*(g(Ib, e)* ...
    (-Rb/Vt)*(1-(Vbc(Ic, Ib, e))/(190)) ...
    + (g(Ib,e)-1)*(Rb/190) - ...
    1.33333*(h(Ic,Ib,e)*(-Rb/Vt)));
dF2dIc = @(Ic, Ib, e) -Is*((2/15)* ...
    (h(Ic, Ib, e))*(Rl/Vt));
dF2dIb = @(Ic, Ib, e) 1 - Is*((1/252)*(g(Ib,e)-1) ...
    *(-Rb/Vt) + 2/15* ...
    (h(Ic, Ib, e))*(-Rb/Vt));

%% Iteration over e(t)
% Amplitude and frequency of e(t)
A = input("Enter the amplitude of the signal (V): ");
f = 1.8e3;

% Define small signal variation e(t) and num of samples per cycle
e = @(t) A*sin(2*pi*f*t);
numCycles = 3;
numSamples = 1000;
sampleTime = numCycles*(1/f);

% Range of time values to iterate over
t = linspace(0, sampleTime, numCycles*numSamples);

%% Initial guess and algorithm parameters
x0 = [10.7672e-3; 42.7272e-6]; % Initial guess: [Ic; Ib]
tol = 1e-6;        %  Tolerance for convergence
maxIter = 20;      % Maximum iterations

% Store results of iteration
IcVals = zeros(size(t)); 
IbVals = zeros(size(t)); 
VoVals = zeros(size(t)); 
iterationVals = zeros(size(t)); 

%% Iterate over the time samples, find Ic and Ib at each sample, find output voltage from those values
for i = 1:length(t)
    disp("*****NEW TIME SAMPLE*****")
    currentTime = t(i);
    e_t = e(currentTime);

    f = @(x) bjtSystem(x, E, e_t, Vcc, Rl, Rb);
    J = @(x) [dF1dIc(x(1), x(2), e_t), dF1dIb(x(1), x(2), e_t);
              dF2dIc(x(1), x(2), e_t), dF2dIb(x(1), x(2), e_t)];

    [x,iterations] = NewtonRaphson(f, J, x0, tol, maxIter);

    % Update solution vectors
    IcVals(i) = x(1); % Store collector current
    IbVals(i) = x(2); % Store base current
    VoVals(i) = Rl*x(1);
    iterationVals(i) = iterations;

    % Use the current solution as the initial guess for the next time step
    x0 = x;
    
end

%% Improved Plots of Output Voltage and Input Signal
t_ms = t * 1e3; % Convert time to milliseconds

figure;

% Subplot 1: Output Voltage
subplot(3, 1, 1);
plot(t_ms, VoVals, 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Output Voltage (V)', 'FontSize', 12);
title('Output Voltage vs. Time', 'FontSize', 20);
legend('V_o(t)', 'Location', 'best');

% Subplot 2: Input Signal
subplot(3, 1, 2);
plot(t_ms, e(t), 'LineWidth', 1.5, 'Color', 'red');
grid on;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Input Signal e(t) (V)', 'FontSize', 12);
title('Input Signal vs. Time', 'FontSize', 20);
legend('e(t)', 'Location', 'best');

% Subplot 3: Output and Input on Same Axis
subplot(3, 1, 3);
plot(t_ms, VoVals, 'LineWidth', 1.5);
hold on;
plot(t_ms, e(t), 'LineWidth', 1.5, 'Color', 'red');
hold off;
grid on;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Voltage (V)', 'FontSize', 12);
title('Output and Input Signals', 'FontSize', 20);
legend('V_o(t)', 'e(t)', 'Location', 'best');


ylims = input("Enter the Y-axis limits (V) as absolute value (limit will be set at ±input) : ");

figure;
plot(t_ms, VoVals-VoVals(1), 'LineWidth', 1.5);
hold on;
plot(t_ms, e(t), 'LineWidth', 1.5, 'Color', 'red');
hold off;
grid on;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Voltage (V)', 'FontSize', 12);
ylim([-ylims ylims])
title('Output and Input Signals with DC offset removed', 'FontSize', 15);
legend('V_o(t)', 'e(t)', 'Location', 'best');

%% Input and output voltages, with iterations
figure;

% First y-axis, voltage
yyaxis left;
plot(t_ms, VoVals - VoVals(1), 'LineWidth', 1.5, 'Color', 'b'); % Blue for Vo
hold on;
plot(t_ms, e(t), 'LineWidth', 1.5, 'Color', 'r'); % Red for e(t)
ylim([-ylims, ylims]); 
ylabel('Voltage (V)', 'FontSize', 12);
ax = gca;
ax.YColor = 'b'; % match left axis color to voltage

% Second y-axis (iterations)
yyaxis right;
plot(t_ms, iterationVals, 'LineWidth', 1.5, 'Color', 'g'); 
ylabel('Iterations', 'FontSize', 12);
ax = gca;
ax.YColor = 'g'; % match right axis color to iterations

hold off;
grid on;
xlabel('Time (ms)', 'FontSize', 12);
title('Output and Input Signals with DC Offset Removed, with Number of Iterations', 'FontSize', 15);
legend({'V_o(t) (DC Removed)', 'e(t) (Input Signal)', 'Newton-Raphson Iterations'}, ...
    'Location', 'best');











