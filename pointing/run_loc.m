%% Clean up console and variables
clc; close all; clear all;
addpath('classes');
addpath('utilities');

%% Raw Data Log Folder
logfolder = 'C:\Users\Alanwar\Dropbox\Pointing\test03\';

%% Node/Network configuration
configfile = 'config/nodepositions_nesl_mobile';

%% Create Network Manager
% NetworkManager(nodeconfig, datafolder, <owr_corrections>, <twr_corrections>) 
iot = IotManager(configfile, logfolder, 'config/antennacorrections_mocap', 'config/twrcorrections_mocap_new', 'config/twrcorrections_mocap_new' );
node_ids = iot.getNodeIds();
node_names = iot.getNodeNames();

%% Rigid body IDs
iot.setRigidBodyId('ntb-mobile', 1);

%% Select reference node
iot.setReferenceNode('ntb-alpha');

%% Measurement set up
iot.enableMessageType(iot.MSGTYPE3, true);
iot.enableMeasurementType(iot.MEASTYPE_R, true);

%% Get the true pointing events
[point_idxs, point_times, point_start, point_stop, anchor_ids, anchor_angles] = getPointEvents(iot);

%% For each pointing segment, gather the range data from each anchor node
SRCID = 8;
WLEN = 0.6;

true_features = [];
true_traces = {};
false_features = [];
false_traces = {};

for i=1:length(point_idxs)
    % pointing info
    xyz_start = point_start(i,:)';
    xyz_stop = point_stop(i,:)';
    aid = anchor_ids(i);
    angle = anchor_angles(i);
    ptime = point_times(i);
    
    anc_xyz = [];
    anc_rng_s = [];
    anc_rng_f = [];
    
    % get ranges for each anchor
    for j=1:length(iot.nodes)
        n = iot.nodes{j};
        [~,rs] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+0.3 );
        [~,rf] = iot.getRangesInWindow( SRCID, n.id, ptime+0.4, ptime+0.7 );
        
        if ~isempty(rs) && ~n.isMobile()
            anc_xyz = [anc_xyz; n.true_p'];
            anc_rng_s = [anc_rng_s; median(rs)];
            anc_rng_f = [anc_rng_f; median(rf)];
        end
    end
    
    % least squares for start and stop location
    fcns = @(xyz) sqrt(sum((repmat(xyz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_s;
    fcnf = @(xyz) sqrt(sum((repmat(xyz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_f;

    LB = [-10 -10 -10];
    UB = [10 10 10];
    xyz_s = lsqnonlin(fcns, [0 0 0], LB, UB, optimset('display','off'));
    xyz_f = lsqnonlin(fcnf, [0 0 0], LB, UB, optimset('display','off'));
    
    % find anchor whose angle deviates the least from the vector defined by
    % start -> finish (ignoring vertical Y axis--too noisy)
    ps = [xyz_s(1) 0 xyz_s(3)];
    pf = [xyz_f(1) 0 xyz_f(3)];
    v1 = (pf-ps)/norm(pf-ps);
    minAngle = inf;
    minId = -1;
    
    for j=1:length(iot.nodes)
        n = iot.nodes{j};
        panc = [n.true_p(1) 0 n.true_p(3)];
        v2 = (panc-ps)./norm(panc-ps);
        dAngle = atan2(norm(cross(v1,v2)), dot(v1,v2));
        if dAngle < minAngle
            minAngle = dAngle;
            minId = n.id;
        end
    end
    
    fprintf('minId: %d (%d)\n', minId, aid);

end


