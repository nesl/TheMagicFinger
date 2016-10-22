%% Clean up console and variables
clc; close all; clear all;
addpath('classes');
addpath('utilities');

%% Raw Data Log Folder
%logfolder = 'C:\Users\Alanwar\Dropbox\Pointing\test03\';
logfolder = 'C:\Users\Alanwar\Dropbox\gesture\';

%% Node/Network configuration
configfile = 'config/nodepositions_nesl_mobile_gesture';

%% Create Network Manager
% NetworkManager(nodeconfig, datafolder, <owr_corrections>, <twr_corrections>) 
iot = IotManager(configfile, logfolder, 'config/antennacorrections_mocap', 'config/twrcorrections_mocap_new', 'config/twrcorrections_mocap_new' );
node_ids = iot.getNodeIds();
node_names = iot.getNodeNames();

%% Rigid body IDs
iot.setRigidBodyId('ntb-mobile', 2);

%% Select reference node
iot.setReferenceNode('ntb-alpha');

%% Measurement set up
iot.enableMessageType(iot.MSGTYPE3, true);
iot.enableMeasurementType(iot.MEASTYPE_R, true);

%% Get the true pointing events 
% anchor ids are the correct chosen IDS with these anchor angles
% The point idexs in the mocap.cvs
% point_start : xyz of the start point.
% point stop : xyz of the end point.
[peaktimes_z_p, peaktimes_z_n, peaktimes_x_p, peaktimes_x_n, locs_z_p, locs_z_n,locs_x_p,locs_x_n] = getGestureEvents(iot);

numofpoints = length(peaktimes_z_p) + length(peaktimes_z_n) + length(peaktimes_x_p) + length(peaktimes_x_n);
SRCID = 9;
WLEN = 0.6;

true_features = [];
true_label =[];
%% positive z -- zp
R1=0.5;
R2=0.9;
for i=1:length(peaktimes_z_p)
    % pointing info
%    xyz_start = point_start(i,:)';
%    xyz_stop = point_stop(i,:)';
%    aid = anchor_ids(i);
%    angle = anchor_angles(i);
    ptime = peaktimes_z_p(i);
    
    anc_xyz = [];
    anc_rng_s = [];
    anc_rng_f = [];
    
    % get ranges for each anchor 
    for j=1:length(iot.nodes)
        n = iot.nodes{j};%SRCID mobile node id
        %looping over the actual ranging from mobile node id 8 to all others
        %ranging start rs, ranging end rf
        [~,rs] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+R1 );
        [~,rf] = iot.getRangesInWindow( SRCID, n.id, ptime+R1, ptime+R2 );
        
        if ~isempty(rs) && ~n.isMobile() && ~isempty(rf)
            anc_xyz = [anc_xyz; n.true_p'];
            anc_rng_s = [anc_rng_s; median(rs)];
            anc_rng_f = [anc_rng_f; median(rf)];
        end
    end
    
    % least squares to get the start and stop location of pointing action
    fcns = @(xyz) sqrt(sum((repmat(xyz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_s;
    fcnf = @(xyz) sqrt(sum((repmat(xyz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_f;

    if(i==1)
        i
    end
    LB = [-10 -10 -10]; %lower bound
    UB = [10 10 10]; %upper bound
    xyz_s = lsqnonlin(fcns, [0,0,0], LB, UB, optimset('display','off'));
    xyz_f = lsqnonlin(fcnf, [0,0,0], LB, UB, optimset('display','off'));
    
    % find anchor whose angle deviates the least from the vector defined by
    % start -> finish (ignoring vertical Y axis--too noisy)
    ps_z_p = [xyz_s(1) 0 xyz_s(3)];
    pf_z_p = [xyz_f(1) 0 xyz_f(3)];
    
    feature = [ xyz_f(1) - xyz_s(1) , xyz_f(3) -  xyz_s(3)];
        % -- store features --
        % zp 0 , zn 1 , xp 2, xn 3, none 4
    true_features = [true_features; feature];
    true_label = [true_label; 0];
      
end

%% negative z -- zn
for i=1:length(peaktimes_z_n)
    % pointing info
%    xyz_start = point_start(i,:)';
%    xyz_stop = point_stop(i,:)';
%    aid = anchor_ids(i);
%    angle = anchor_angles(i);
    ptime = peaktimes_z_n(i);
    
    anc_xyz = [];
    anc_rng_s = [];
    anc_rng_f = [];
    
    % get ranges for each anchor 
    for j=1:length(iot.nodes)
        n = iot.nodes{j};%SRCID mobile node id
        %looping over the actual ranging from mobile node id 8 to all others
        %ranging start rs, ranging end rf
        [~,rs] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+R1 );
        [~,rf] = iot.getRangesInWindow( SRCID, n.id, ptime+R1, ptime+R2 );
        
        if ~isempty(rs) && ~n.isMobile() && ~isempty(rf)
            anc_xyz = [anc_xyz; n.true_p'];
            anc_rng_s = [anc_rng_s; median(rs)];
            anc_rng_f = [anc_rng_f; median(rf)];
        end
    end
    
    % least squares to get the start and stop location of pointing action
    fcns = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_s;
    fcnf = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_f;

    LB = [-10 -10 -10]; %lower bound
    UB = [10 10 10]; %upper bound

    xyz_s = lsqnonlin(fcns, [-2.2036, 1.2664, -2.5278], LB, UB, optimset('display','off'));
    xyz_f = lsqnonlin(fcnf, [-2.2036, 1.2664, -2.5278], LB, UB, optimset('display','off'));
    
    % find anchor whose angle deviates the least from the vector defined by
    % start -> finish (ignoring vertical Y axis--too noisy)
    ps_z_p = [xyz_s(1) 0 xyz_s(3)];
    pf_z_p = [xyz_f(1) 0 xyz_f(3)];
    
    feature = [ xyz_f(1) - xyz_s(1) , xyz_f(3) -  xyz_s(3)];
        % -- store features --
        % zp 0 , zn 1 , xp 2, xn 3, none 4
    true_features = [true_features; feature];
    true_label = [true_label; 1];

end

%% positive  x -- xp
for i=1:length(peaktimes_x_p)
    % pointing info
%    xyz_start = point_start(i,:)';
%    xyz_stop = point_stop(i,:)';
%    aid = anchor_ids(i);
%    angle = anchor_angles(i);
    ptime = peaktimes_x_p(i);
    
    anc_xyz = [];
    anc_rng_s = [];
    anc_rng_f = [];
    
    % get ranges for each anchor 
    for j=1:length(iot.nodes)
        n = iot.nodes{j};%SRCID mobile node id
        %looping over the actual ranging from mobile node id 8 to all others
        %ranging start rs, ranging end rf
        [~,rs] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+R1 );
        [~,rf] = iot.getRangesInWindow( SRCID, n.id, ptime+R1, ptime+R2 );
        
        if ~isempty(rs) && ~n.isMobile() && ~isempty(rf)
            anc_xyz = [anc_xyz; n.true_p'];
            anc_rng_s = [anc_rng_s; median(rs)];
            anc_rng_f = [anc_rng_f; median(rf)];
        end
    end
    
    % least squares to get the start and stop location of pointing action
    fcns = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_s;
    fcnf = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_f;

    LB = [-10 -10 -10]; %lower bound
    UB = [10 10 10]; %upper bound
    xyz_s = lsqnonlin(fcns, [-2.2036, 1.2664, -2.5278], LB, UB, optimset('display','off'));
    xyz_f = lsqnonlin(fcnf, [-2.2036, 1.2664, -2.5278], LB, UB, optimset('display','off'));
    
    % find anchor whose angle deviates the least from the vector defined by
    % start -> finish (ignoring vertical Y axis--too noisy)
    ps_z_p = [xyz_s(1) 0 xyz_s(3)];
    pf_z_p = [xyz_f(1) 0 xyz_f(3)];
    
    feature = [ xyz_f(1) - xyz_s(1) , xyz_f(3) -  xyz_s(3)];
        % -- store features --
        % zp 0 , zn 1 , xp 2, xn 3, none 4
    true_features = [true_features; feature];
    true_label = [true_label; 2];

end

%% negative  x -- xn
for i=1:length(peaktimes_x_n)
    % pointing info
%    xyz_start = point_start(i,:)';
%    xyz_stop = point_stop(i,:)';
%    aid = anchor_ids(i);
%    angle = anchor_angles(i);
    ptime = peaktimes_x_n(i);
    
    anc_xyz = [];
    anc_rng_s = [];
    anc_rng_f = [];
    
    % get ranges for each anchor 
    for j=1:length(iot.nodes)
        n = iot.nodes{j};%SRCID mobile node id
        %looping over the actual ranging from mobile node id 8 to all others
        %ranging start rs, ranging end rf
        [~,rs] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+R1 );
        [~,rf] = iot.getRangesInWindow( SRCID, n.id, ptime+R1, ptime+R2 );
        
        if ~isempty(rs) && ~n.isMobile() && ~isempty(rf)
            anc_xyz = [anc_xyz; n.true_p'];
            anc_rng_s = [anc_rng_s; median(rs)];
            anc_rng_f = [anc_rng_f; median(rf)];
        end
    end
    
    % least squares to get the start and stop location of pointing action
    fcns = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_s;
    fcnf = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_f;

    LB = [-10 -10 -10]; %lower bound
    UB = [10 10 10]; %upper bound
    xyz_s = lsqnonlin(fcns, [-2.2036, 1.2664, -2.5278], LB, UB, optimset('display','off'));
    xyz_f = lsqnonlin(fcnf, [-2.2036, 1.2664, -2.5278], LB, UB, optimset('display','off'));
    
    % find anchor whose angle deviates the least from the vector defined by
    % start -> finish (ignoring vertical Y axis--too noisy)
    ps_z_p = [xyz_s(1) 0 xyz_s(3)];
    pf_z_p = [xyz_f(1) 0 xyz_f(3)];
    
    feature = [ xyz_f(1) - xyz_s(1) , xyz_f(3) -  xyz_s(3)];
        % -- store features --
        % zp 0 , zn 1 , xp 2, xn 3, none 4
    true_features = [true_features; feature];
    true_label = [true_label; 3];

end


true=0;
false=0;
for i=1:length(true_label)
    dx=true_features(1);
    dz=true_features(2);
    %zp
   if (true_label(i)== 0)
       if ( dz > 0 && abs(dz) > abs(dx))
         true  = true+1;
       else
         false = false+1;
       end
   end
    %zn
   if (true_label(i)== 1)
       if ( dz < 0 && abs(dz) > abs(dx))
         true  = true+1;
       else
         false = false+1;
       end
   end   
   %xp
   if (true_label(i)== 2)
       if ( dx > 0 && abs(dz) < abs(dx))
         true  = true+1;
       else
         false = false+1;
       end
   end
           
  %xp
   if (true_label(i)== 3)
       if ( dx < 0 && abs(dz) < abs(dx))
         true  = true+1;
       else
         false = false+1;
       end
   end         
end
fprintf('True  =%d\n',true);
fprintf('False =%d\n',false);

true2=0;
false2=0;
for i=1:length(true_label)
    dx=true_features(1);
    dz=true_features(2);
    %zp
   if (true_label(i)== 0 || true_label(i)==1)
       if (  abs(dz) > abs(dx))
         true2  = true2+1;
       else
         false2 = false2+1;
       end
   end
   
   if (true_label(i)== 2 || true_label(i)==3)
       if (  abs(dz) < abs(dx))
         true2  = true2+1;
       else
         false2 = false2+1;
       end
   end
end

fprintf('True2  =%d\n',true2);
fprintf('False2 =%d\n',false2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train SVM
all_features = true_features;
all_labels = true_label;
N = size(all_features,1) ;
P = cvpartition(N, 'Holdout', 0.30);
%svmStruct = svmtrain( all_features(P.training,:), all_labels(P.training,:), 'kernel_function', 'poly');
RF = TreeBagger(100,all_features(P.training,:),all_labels(P.training,:), 'Method', 'regression', 'oobpred','on');

disp('===== Non-cooperative detection =====');
fprintf('Training on %d, Testing on %d\n', length(find(P.training)), length(find(P.test)));

%% Test SVM
% [C,F] = svmclassifydist(svmStruct, all_features(P.test,:));
% 
% n_truepos = sum( all_labels(P.test) == C & all_labels(P.test) == 1 );
% n_trueneg = sum( all_labels(P.test) == C & all_labels(P.test) == 0 );
% n_falsepos = sum( all_labels(P.test) ~= C & C == 1 );
% n_falseneg = sum( all_labels(P.test) ~= C & C == 0 );
% 
% precision: how many false positives
% precision = n_truepos / (n_truepos + n_falsepos);
% recall: how many positives missed
% recall = n_truepos / (n_truepos + n_falseneg);
% 
% fprintf('===== SVM =====\n');
% fprintf('Precision: %.2f\n', precision);
% fprintf('Recall: %.2f\n', recall);
% fprintf('F score: %.2f\n', 2*precision*recall/(precision + recall));

% %%{
Cvals = RF.predict(all_features(P.test,:));
C = Cvals > 0.20;
n_truepos = sum( all_labels(P.test) == C & all_labels(P.test) == 1 );
n_trueneg = sum( all_labels(P.test) == C & all_labels(P.test) == 0 );
n_falsepos = sum( all_labels(P.test) ~= C & C == 1 );
n_falseneg = sum( all_labels(P.test) ~= C & C == 0 );

% precision: how many false positives
precision = n_truepos / (n_truepos + n_falsepos);
% recall: how many positives missed
recall = n_truepos / (n_truepos + n_falseneg);

fprintf('===== RF =====\n');
fprintf('Precision: %.2f\n', precision);
fprintf('Recall: %.2f\n', recall);
fprintf('F score: %.2f\n', 2*precision*recall/(precision + recall));
%%}
