%{
    Function to calculate thermodynamic constants from a temperature
%}

function [c, rho, eta, nu, gamma] = calcThermoDynConstants(T)
deltaT = T - 26.85;
c = 3.4723e2 * (1 + 0.00166 * deltaT);      % Speed of sound in air [m/s]
rho = 1.1769 * (1 - 0.00335 * deltaT);      % Density of air [kg·m^{-3}]
eta = 1.846 * (1 + 0.0025 * deltaT);        % Shear viscosity [kg·s^{-1}·m^{-1}]
nu = 0.8410 * (1 - 0.0002 * deltaT);        % Root of Prandtl number [-]
gamma = 1.4017 * (1 - 0.00002 * deltaT);    % Ratio of specific heats [-]
