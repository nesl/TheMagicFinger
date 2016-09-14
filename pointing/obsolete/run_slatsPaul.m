%% Housekeeping
clc; close all; clear all;
addpath('classes');
addpath('utilities');

%% Raw Data Log Folder
fp = '../../../data/multitwr_30min/multitwr';

%% Node configuration
cf = 'config/nodepositions_nesl';

%% Source Node ID
srcNode = 'ntb-golf';

%% Create data parser
dp = TimingDataParser(cf, fp);
nodeInfo = dp.getNodeInfo();

return;

%% EKF Covariances
% Process Noise
% x = [x y z dx dy dz co1 cs1 co2 cs2 ... coN csN]
N = length(nodeInfo{1});
var_xyz = 0;
var_co = 1e4; % ns
var_cs = 1e3; % ns / s
var_proc = [ var_xyz var_xyz var_xyz ];
for i=1:N
    var_proc = [var_proc var_co var_cs];
end
Q = diag(var_proc);

% Measurement Noise
var_tx = 0.50; % ns
var_rx = 0.50; % ns
var_meas = var_tx;
for i=1:N
    var_meas = [var_meas var_rx];
end
R = diag(var_meas);

% Initial Covariance
vari_xyz = [10 10 2.419];
vari_co = 1e10; % ns
vari_cs = 4e3; % ns
var_init = [vari_xyz(1) vari_xyz(2) vari_xyz(3)];
for i=1:N
    var_init = [var_init vari_co vari_cs];
end
P = diag(var_init);

% Initial State Estimate
x = [0; 0; 0];
for i=1:N
    x = [x; 0; 0];
end

%% Replay data and run EKF
d1 = dp.getSequenceData(1);
t_txlast = d1(1)/1e9; % sec
x_history = [];
z_history = [];

for seq=2:dp.getSequenceCount()
    
    % at sequence 1000, constrain the unknown node's position
    if seq == 1000 && 0
        x(1:3) = [3.440, 4.263, 0.774];
        P(1,1) = 0;
        P(2,2) = 0;
        P(3,3) = 0;
    else
        % get measured data for this sequence
        [meas, ids] = dp.getSequenceData(seq);
        % append measurements to z history
        z_history = [z_history; meas'];
        % get non-zero measurements and corresponding node IDs
        z_ids = ids(find(meas));
        z_pos = [];
        for i=1:length(z_ids)
            id = z_ids(i);
            z_pos = [z_pos; dp.getNodePos(id)];
        end
        z = meas(find(meas));
        
        % Update R covariance matrix based on observed measurements
        R = diag([var_tx var_rx*ones(1,length(z_ids)-1)]);
        
        % get delta time
        t_tx = z(1);
        dt = (t_tx/1e9) - t_txlast;
        t_txlast = (t_tx/1e9);
        
        % strap process and measurement functions with time + node info
        f = @(x) ekf_procPaul(x,dt);
        h = @(x) ekf_measPaul(x,t_tx, z_ids,z_pos);
        
        % update state estimate
        [x, P] = ekf(f, x, P, h, z, dt*Q, R);
        
        fprintf('x: %.1f, y: %.1f, z: %.1f', x(1), x(2), x(3) );
        for i=0:(N-1)
            fprintf(', o[%d]:%d / %d', i, x(4 + 2*i), x(5 + 2*i) );
        end
        fprintf('\n');
        
        % append state estimate to history
        x_history = [x_history; x'];
    end
    
end

%% Plot Markers & Colors
pltopts = {'-o', '-s', '-^', '-x', '-*', '--o', '--s', '--^', '--x', '--*'};
pltcolors = [
    1 0 0
    0 0 1
    0 0 0
    0.8 0.8 0
    1 0 1
    0 0.7 0.7
    0 1 0
    0.5 0.5 0.5
    ];

%% Plot Bias Convergence
cfigure(20,12);
handles = [];
nodenames = nodeInfo{1};

for i=1:length(nodenames)
    h = plot( x_history(1:1:end,5 + 2*(i-1)), pltopts{i}, 'Color', pltcolors(i,:), 'MarkerSize', 4);
    hold on;
    handles = [handles; h];
end

grid on;
xlabel('Sequence Number', 'FontSize', 14);
ylabel('Bias (ppb)', 'FontSize', 14);
legend(handles, nodenames, 'Location', 'SE');
saveplot('output/bias_convergence');

%% Plot Offset Divergence
cfigure(20,12);
handles = [];
nodenames = nodeInfo{1};

for i=1:length(nodenames)
    ofst = x_history(:, 4 + 2*(i-1));
    h = plot( ofst - ofst(1), pltopts{i}, 'Color', pltcolors(i,:) );
    hold on;
    handles = [handles; h];
end

grid on;
xlabel('Sequence Number', 'FontSize', 14);
ylabel('Relative Offset (ns)', 'FontSize', 14);
legend(handles, nodenames, 'Location', 'SW');
saveplot('output/offset_divergence');













