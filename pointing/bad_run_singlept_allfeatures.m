%% Clean up console and variables
% long right 3254032540
clc; close all; clear all;
addpath('classes');
addpath('utilities');

%% Raw Data Log Folder
%logfolder = 'C:\Users\Alanwar\Dropbox\Pointing\test01\';
logfolder = 'C:\Users\Alanwar\Dropbox\gesture\newpointing\';

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
[point_idxs, point_times, point_start, point_stop, anchor_ids, anchor_angles] = getSinglePoints(iot);

%point_start start in x,y,z
%%{
figure();
hold on;
scatter(point_start(:,1), point_start(:,3), 'bo');
scatter(point_stop(:,1), point_stop(:,3), 'rs');
scatter(point_start(1,1), point_start(1,3), 'g^');

for i=1:length(iot.nodes)
    n = iot.nodes{i};
    xyz = n.true_p;
    text(xyz(1), xyz(3)+0.2, num2str(n.id));
    plot(xyz(1), xyz(3), 'ms','MarkerSize',10,'MarkerFaceColor', 'm');
end
xlabel('X');
ylabel('Z');
grid on;
axis equal;
%%}

%% For each pointing segment, gather the range data from each anchor node
SRCID = 9;
WLEN = 0.6;

all_features = zeros(length(point_idxs),45);
true_traces = {};
false_features = [];
false_traces = {};
true_label =[];
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
        n = iot.nodes{j};%SRCID mobile node id
        %looping over the actual ranging from mobile node id 8 to all others
        %ranging start rs, ranging end rf
        [~,rs] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+0.3 );
        [~,rf] = iot.getRangesInWindow( SRCID, n.id, ptime+0.4, ptime+0.7 );
        
        if ~isempty(rs) && ~n.isMobile() && ~isempty(rf)
            if ~isnan(median(rs)) && ~isnan(median(rf))
            anc_xyz = [anc_xyz; n.true_p'];
            anc_rng_s = [anc_rng_s; median(rs)];
            anc_rng_f = [anc_rng_f; median(rf)];
            end
        end
    end
    
    % least squares to get the start and stop location of pointing action
    fcns = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_s;
    fcnf = @(xyzz) sqrt(sum((repmat(xyzz,size(anc_xyz,1),1) - anc_xyz).^2, 2)) - anc_rng_f;

    LB = [-10 -10 -10]; %lower bound
    UB = [10 10 10]; %upper bound
    xyz_s = lsqnonlin(fcns, [0 0 0], LB, UB, optimset('display','off'));
    xyz_f = lsqnonlin(fcnf, [0 0 0], LB, UB, optimset('display','off'));
    
    % find anchor whose angle deviates the least from the vector defined by
    % start -> finish (ignoring vertical Y axis--too noisy)
    ps = [xyz_s(1) 0 xyz_s(3)];
    pf = [xyz_f(1) 0 xyz_f(3)];
    v1 = (pf-ps)/norm(pf-ps);
    minAngle = inf;
    minId = -1;
    features = [];
    for j=1:length(iot.nodes)
        n = iot.nodes{j};
        [t,r] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+WLEN );
        [fploss] = iot.getFplossInWindow( SRCID, n.id, ptime, ptime+WLEN );
        
         if length(t) < 2
             fvec= [0 0 0 0 0];
             features = [features , fvec];
             continue;
        end
         
        % sanity check for range errors
%         if max(r) > 20
%              fvec= [0.1 0.1 0.1 0.1 0.1];
%              features = [features , fvec];
%             continue;
%         end
        
        
        % -- calculate additional features --
        % delta_range
        dr = r(end)-r(1);
        % linear fit
        p = polyfit(t-t(1), r, 1);
        residuals = polyval(p,t-t(1)) - r;
        noise = mean(residuals.^2);
        % angle change
        panc = [n.true_p(1) 0 n.true_p(3)];
        v2 = (panc-ps)./norm(panc-ps);
        dAngle = atan2(norm(cross(v1,v2)), dot(v1,v2));
        % delta_range, start_range, lin fit, path loss
        fvec = [dr, r(1), noise, mean(fploss), dAngle];
        
        % -- store features --
        %if aid == n.id
        features = [features , fvec];
    end
    
    all_features(i,:) = features;
    true_label = [ true_label; aid];
    
end

N = size(all_features,1);
results = multisvm( all_features(1:70,:), true_label(1:70), all_features(71:end,:) );

right=0;
wrong=0;
for i=1:length(results)
   if(results(i)==true_label(70+i) )
    right=right+1;
   else
    wrong=wrong+1;
   end
end


% nCorrect = sum(aid_guesses == anchor_ids);
fprintf('svm Accuracy: %f\n', right/length(results));

RF = TreeBagger(100,all_features(1:50,:), true_label(1:50), 'Method', 'regression', 'oobpred','on');
C = RF.predict( all_features(51:end,:) );

results = ceil(C);
right=0;
wrong=0;
for i=1:length(results)
   if(results(i)==true_label(50+i) )
    right=right+1;
   else
    wrong=wrong+1;
   end
end

% nCorrect = sum(aid_guesses == anchor_ids);
fprintf('RF Accuracy: %f\n', right/length(results));


%%
% cfigure(20,10);
% i=1;
% for i=1:size(false_traces,1)
%    hf = plot(false_traces{i} - false_traces{i}(1),'ro-');
%    hold on;
% end
% 
% 
% for i=1:size(true_traces,1)
%     ht = plot(true_traces{i} - true_traces{i}(1), 'bs-');
%     hold on;
% end
% 
% ylim([-1 1]);
% ylabel('Change in Range Estimate (m)');
% xlabel('Index Number');
% grid on;
% legend([ht,hf], {'True Traces', 'False Traces'});

%% Compare true and false delta ranges
% cfigure(20,12);
% [c,f] = calculatePdf(true_features(:,1), 15);
% ht = plot(c,f,'bs-');
% hold on;
% [c,f] = calculatePdf(false_features(:,1), 15);
% hf = plot(c,f,'ro-');
% 
% xlabel('Range Estimate Change (m)');
% ylabel('Frequency');
% grid on;
% legend([ht,hf], {'True \DeltaR', 'False \DeltaR'});

