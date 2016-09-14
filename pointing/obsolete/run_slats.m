close all; clear; %pkg load statistics;
addpath('utilities');

load('data.mat');
N = size(positions,1);
D = squareform(pdist(positions));

NSEC_TO_SEC = 1e-9;
SEC_TO_NSEC = 1e9;

% Parameters
A       = [10 10 2.419];  % Area
offset  = 18e9;           % Offset (ns)
bias    = 4000;           % Bias (ns)
noise   = 0.8;            % Gaussian jitter
grw     = 1.6;            % Bias random walk

% Assumption: sensor (TX) clock is perfect, anchor (RX) clocks drift
% State: [sensor_x, sensor_y, offset_1, bias_1, ..., offset_N, bias_N]

% Last update time
T = zeros(1,N);

% State
x = zeros(3+2*N, 1);  
P = eye(length(x));
Q = eye(length(x));
x(1:3) = truth;
P(1,1) = A(1)^2;
P(2,2) = A(2)^2;
P(3,3) = A(3)^2;
Q(1,1) = 0;
Q(2,2) = 0;
Q(3,3) = 0;
for i = 1:N
  x(3+2*(i-1)+1) = 0;                       % Offset
  x(3+2*(i-1)+2) = 0;                       % Bias
  P(3+2*(i-1)+1,3+2*(i-1)+1) = offset^2;    % Offset variance ~17 secs because of overflow
  P(3+2*(i-1)+2,3+2*(i-1)+2) = bias^2;      % Anchor bias
  Q(3+2*(i-1)+1,3+2*(i-1)+1) = 0;           % Anchor offset
  Q(3+2*(i-1)+2,3+2*(i-1)+2) = grw^2;       % Random walk (ADEV for tau=1)
end

% Iterate over dataset
for i = 1:size(prx,1)
  
  % Only consider this data if there is something meaningful there
  idx = find(ptx(i,:) ~= 0);
  if ~isempty(idx)
    
    % Takes care of first time epoch for each measurement axis
    dt = ptx(i,idx) - T(idx);
    
    % Only if we have some meaningful data
    ndx = idx(dt < 2*SEC_TO_NSEC);

    % Check that we have measurements, and that all anchors have previous timesteps
    if (~isempty(ndx) && all(T > 0))
      
      % Takes care of first time epoch for each measurement axis
      dt = ptx(i,ndx) - T(ndx);
      dt = dt(1) / 1e9;

      % Measurement and noise covariance
      z = prx(i,ndx)' - ptx(i,ndx)';
      R = noise^2 * eye(length(z));
      
      % Process model - 
      f = @(x) ekf_proc(x,dt);
      h = @(x) ekf_meas(x,positions,ndx);
      
      % Perform filtering
      [x, P] = ekf(f, x, P, h, z, dt*Q, R);        

      % Debug
      disp(['[' num2str(i) '] Solution: X=' num2str(x(1)) ', Y=' num2str(x(2)) ', Z=' num2str(x(3)) ',  O1=' num2str(x(4)) ', B1=' num2str(x(5)) ', O2=' num2str(x(6)) ', B2=' num2str(x(7)) ...
        ', O3=' num2str(x(8)) ', B3=' num2str(x(9)) ', O4=' num2str(x(10)) ', B4=' num2str(x(11)) ', O5=' num2str(x(12)) ', B5=' num2str(x(13)) ', O6=' num2str(x(14)) ', B6=' num2str(x(15))]);
      
      % Sanity
      %fflush(stdout);
      pause(0.1);

    end

    % Cache update times
    T(idx) = ptx(i,idx);

  end
  
end