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
    
    for j=1:length(iot.nodes)
        n = iot.nodes{j};
        [t,r] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+WLEN );
        [fploss] = iot.getFplossInWindow( SRCID, n.id, ptime, ptime+WLEN );
        
        if length(t) < 2
            continue;
        end
        
        % sanity check for range errors
        if max(r) > 20
            continue;
        end
        
        
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
        %fvec = [dr, r(1), noise, mean(fploss), dAngle];
        fvec = [dr, r(1), noise, mean(fploss)];
        % -- store features --
        if aid == n.id
            true_features = [true_features; fvec];
            true_traces = [true_traces; r];
        else
            false_features = [false_features; fvec];
            false_traces = [false_traces; r];
        end
        
    end
end

%% Compare true and false angle vs dR
%{
cfigure(20,12);
ht = plot(true_features(:,1), true_features(:,5),'bs');
hold on;
hf = plot(false_features(:,1), false_features(:,5),'ro');

xlabel('Estimated Range Change (m)');
ylabel('Estimated Angle Change (rad)');
grid on;
legend([ht,hf], {'True \Delta\theta', 'False \Delta\theta'});
%}

%% Compare true and false angles
%
% cfigure(20,12);
% [c,f] = calculatePdf(true_features(:,5), 15);
% ht = plot(c,f,'bs-');
% hold on;
% [c,f] = calculatePdf(false_features(:,5), 15);
% hf = plot(c,f,'ro-');
% 
% xlabel('Estimated Angular Change (rad.)');
% ylabel('Frequency');
% grid on;
% legend([ht,hf], {'True \Delta\theta', 'False \Delta\theta'});



% cfigure(20,10);
% max_t=0;
% max_t_index=1;
% min_t=0;
% min_t_index=1;
% mean_t_index=1;
% mean_t=1;
% for i=1:size(false_traces,1)
%     if(max_t < max(false_traces{i} ))
%         max_t =max(false_traces{i});
%         max_t_index = i;
%     end
%     if(min_t > min(false_traces{i} ))
%         min_t =min(false_traces{i});
%         min_t_index = i;
%     end
% end
%    
% hf = plot(false_traces{max_t_index} - false_traces{max_t_index}(1),'ro-');
% hold on
% ht = plot(false_traces{min_t_index} - false_traces{min_t_index}(1),'bs-');
%% Compare true and false traces
% for i=1:size(false_traces,1)
%     hf = plot(false_traces{i} - false_traces{i}(1),'ro-');
%     hold on;
% end
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
% cfigure(20,10);
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

%% Compare true and false fploss
% cfigure(20,12);
% [c,f] = calculatePdf(true_features(:,4), 20);
% ht = plot(c,f,'bs-');
% hold on;
% [c,f] = calculatePdf(false_features(:,4), 20);
% hf = plot(c,f,'ro-');
% 
% xlabel('Range Estimate Change (m)');
% ylabel('Frequency');
% grid on;
% legend([ht,hf], {'True \DeltaR', 'False \DeltaR'});


%% Train SVM
all_features = [true_features; false_features];
all_labels = [ones(size(true_features,1),1); zeros(size(false_features,1),1)];
N = size(true_features,1) + size(false_features,1);
P = cvpartition(N, 'Holdout', 0.30);
type = 'rbf';

svmStruct = svmtrain( all_features(P.training,:), all_labels(P.training,:), 'kernel_function', type);
RF = TreeBagger(100,all_features(P.training,:),all_labels(P.training,:), 'Method', 'regression', 'oobpred','on');
disp('===== Non-cooperative detection =====');
fprintf('Training on %d, Testing on %d\n', length(find(P.training)), length(find(P.test)));
%fprintf('Type =%s\n',type);
%% Test SVM
[C,F ] = svmclassifydist(svmStruct, all_features(P.test,:));

n_truepos = sum( all_labels(P.test) == C & all_labels(P.test) == 1 );
n_trueneg = sum( all_labels(P.test) == C & all_labels(P.test) == 0 );
n_falsepos = sum( all_labels(P.test) ~= C & C == 1 );
n_falseneg = sum( all_labels(P.test) ~= C & C == 0 );

% precision: how many false positives
precision = n_truepos / (n_truepos + n_falsepos);
% recall: how many positives missed
recall = n_truepos / (n_truepos + n_falseneg);

fprintf('===== SVM %s=====\n',type);
fprintf('Precision: %f\n', precision);
fprintf('Recall: %f\n', recall);
fprintf('F score: %f\n', 2*precision*recall/(precision + recall));

%%{
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
fprintf('Precision: %f\n', precision);
fprintf('Recall: %f\n', recall);
fprintf('F score: %f\n', 2*precision*recall/(precision + recall));
%%}

%% RF Error vs. Trees
%{
oobErrorBaggedEnsemble = oobError(RF);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';
%}

%% Cooperative detection, minimum overhead
disp('===== Min-cooperative detection =====');
aid_guesses = [];

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
            if ~isnan(median(rs)) && ~isnan(median(rf))
                anc_xyz = [anc_xyz; n.true_p'];
            anc_rng_s = [anc_rng_s; median(rs)];
            anc_rng_f = [anc_rng_f; median(rf)];
            end
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
    
    Cvals = [];
    Cs = [];
    Aids = [];
    
    for j=1:length(iot.nodes)
        n = iot.nodes{j};
        [t,r] = iot.getRangesInWindow( SRCID, n.id, ptime, ptime+WLEN );
        [fploss] = iot.getFplossInWindow( SRCID, n.id, ptime, ptime+WLEN );
        
        if length(t) < 2
            continue;
        end
        
        % sanity check for range errors
        if max(r) > 20
            continue;
        end
        
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
        fvec = [dr, r(1), noise, mean(fploss)];
        %fvec = [dr, r(1), noise, mean(fploss), dAngle];
        
        % classify with RF
        %cval = RF.predict(fvec);
        [c,cval] = svmclassifydist(svmStruct, fvec);

        
        %Cvals = [Cvals; cval];
        Cs = [Cs; c];
        Aids = [Aids; n.id];
        
    end
    
    % collaborative guess
    [~,idx] = max(Cs); %Cvals
    aid_guess = Aids(idx);
    aid_guesses = [aid_guesses; aid_guess];
    
end

nCorrect = sum(aid_guesses == anchor_ids);
fprintf('Accuracy: %f\n', nCorrect/length(anchor_ids));


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

