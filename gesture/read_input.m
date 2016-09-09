function [acc, grav, gyro, mag, rot] = read_input(path_prefix, varargin)

% read_input() reads sensor input. Right now we have two formats, and this
% function helps and eases the effort of reading files.

try
    [path_prefix '_acc.csv']
    raw_acc = csvread([path_prefix '_acc.csv']);
catch e
    raw_acc = zeros(0, 5);
end

try
    raw_gyro = csvread([path_prefix '_gyro.csv']);
catch e
    raw_gyro = zeros(0, 5);
end

try
    raw_grav = csvread([path_prefix '_grav.csv']);
catch e
    raw_grav = zeros(0, 5);
end

try
    raw_mag = csvread([path_prefix '_mag.csv']);
catch e
    raw_mag = zeros(0, 5);
end

try
    raw_rot = csvread([path_prefix '_rot.csv']);
catch e
    raw_rot = zeros(0, 7);
end

% acc/grav, gyro, mag, and rotation vector share the same sensor time system

% processed sensor data, format: shifted_time(s), x, y, z
acc = raw_acc(:, 2:end);
gyro = raw_gyro(:, 2:end);
grav = raw_grav(:, 2:end);
mag = raw_mag(:, 2:end);
rot = raw_rot(:, 2:end);

offset = min([acc(:, 1) ; gyro(:, 1) ; grav(:, 1) ; mag(:, 1) ; rot(:, 1)]);

acc(:, 1)  = (acc(:, 1)  - offset) * 1e-9;
gyro(:, 1) = (gyro(:, 1) - offset) * 1e-9;
grav(:, 1) = (grav(:, 1) - offset) * 1e-9;
mag(:, 1)  = (mag(:, 1)  - offset) * 1e-9;
rot(:, 1)  = (rot(:, 1)  - offset) * 1e-9;

% filter outliers
acc(:, 2:4) = min(30, acc(:, 2:4));
acc(:, 2:4) = max(-30, acc(:, 2:4));
gyro(:, 2:4) = min(30, gyro(:, 2:4));
gyro(:, 2:4) = max(-30, gyro(:, 2:4));
grav(:, 2:4) = min(30, grav(:, 2:4));
grav(:, 2:4) = max(-30, grav(:, 2:4));

end

