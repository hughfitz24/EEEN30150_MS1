%% bjtSystem.m
% M-file creating the function that models the 
% system of equations for the amplifier detailed in the assignment brief.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of MS1. 

function F = bjtSystem(x, E, e, Vcc, Rl, Rb)
    % Variables
    Ic = x(1); % Collector current
    Ib = x(2); % Base current

    % Constants
    Is = 0.96e-12; % Saturation current
    a = 36.681686; % 1/Vt

    F = zeros(2, 1);  % Initialize F as a 2x1 column vector

    % system of equations
    F(1) = Ic - Is * ((exp(a * (E + e - Rb * Ib)) - 1) * ...
          (1 - (E + e -Rb*Ib -Vcc + Rl*Ic) / 190) ...
          - 1.33333 * (exp(a * (E + e -Rb*Ib -Vcc + Rl*Ic)) - 1));

    F(2) = Ib - Is * ((1 / 252) * (exp(a * (E + e - Rb * Ib)) - 1) + ...
          (2 / 15) * (exp(a * (E + e -Rb*Ib -Vcc + Rl*Ic)) - 1));
end