% This is a silly little sript to determine how the gravitational potential
% between two sensors in our network affects their relative biases.
% SEE: https://en.wikipedia.org/wiki/Error_analysis_for_the_Global_Positioning_System

Re = 6357000;       % radius of sensor (floor)
Rn = 6357001;       % radius of sensor (roof)
Rg = 26541000;      % radius of sensor (GPS)
M  = 5.974e24;      % mass of earth
G  = 6.67408e-11;   % gravitational constant
c  = 2.998e8;       % speed of light
T = 60*60*24*1e12;   % time

% Calculate the time dilation
dG = G * M / (Re * c^2) - G * M / (Rg * c^2);
dG = dG * T;
dN = G * M / (Re * c^2) - G * M / (Rn * c^2);
dN = dN * T;

% Display
disp(['GPS: ps drift per day due to time dilation: ' num2str(dG)]);
disp(['NTB: ps drift per day due to time dilation: ' num2str(dN)]);