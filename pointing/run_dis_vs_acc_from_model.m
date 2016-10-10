%% Clean up console and variables
% long right 3254032540
clc; close all; clear all;
addpath('classes');
addpath('utilities');
load('cache\model');
%% Raw Data Log Folder
%logfolder = 'C:\Users\Alanwar\Dropbox\Pointing\test01\';
%logfolder = 'C:\Users\Alanwar\Dropbox\gesture\newpointing\';
%logfolder = 'C:\Users\Alanwar\Dropbox\gesture\Bojang\';
%logfolder = 'C:\Users\Alanwar\Dropbox\gesture\differentpoints\';
logfolder = 'C:\Users\Alanwar\Dropbox\gesture\distance\80\';

%% Node/Network configuration
configfile = 'config/nodepositions_nesl_distance';

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
[point_idxs_old, point_times_old, point_start_old, point_stop_old, anchor_ids_old, anchor_angles_old,all_angles_old] = getSinglePoints_manual(iot);

point_idxs=[];
point_times=[];
point_start=[];
point_stop =[];
anchor_ids=[]; 
anchor_angles=[];
all_angles =[];
index =1;
for i=1:length(point_idxs_old)
    
 %   if(anchor_angles_old(i)*180/pi < 70 && anchor_ids_old(i)~=1)
        point_idxs=[point_idxs;point_idxs_old(i)];
        point_times=[point_times;point_times_old(i)];
        point_start(index,:)=point_start_old(i,:);
        point_stop(index,:)=point_stop_old(i,:);
        anchor_ids=[anchor_ids;anchor_ids_old(i)];
        anchor_angles=[anchor_angles;anchor_angles_old(i)];
        all_angles(index,:) =all_angles_old(i,:);
        index = index +1;
  %  end
    
end
%15cm Accuracy: 0.733333
%anchor_ids = [ 6*ones(1,7) 7*ones(1,8) 6*ones(1,9) 7*ones(1,9)];

%50cm Accuracy: 0.833
%anchor_ids = [ 6*ones(1,7) 7*ones(1,7) 6*ones(1,7) 7*ones(1,9)];

%1m Accuracy: 0.806452
%anchor_ids = [ 6*ones(1,8) 7*ones(1,8) 6*ones(1,7) 7*ones(1,8)];

%1.5m Accuracy: 0.827586  Accuracy: 0.862069
%anchor_ids = [ 6*ones(1,7) 7*ones(1,8) 6*ones(1,7) 7*ones(1,8)];

%2m Accuracy: 0.838710  0.903
%anchor_ids = [ 6*ones(1,10) 7*ones(1,8) 6*ones(1,8) 7*ones(1,7)];

%%
%70 cm
%2m 1.000000 0.857143
%anchor_ids = [ 6*ones(1,7) 7*ones(1,3) 6*ones(1,7) 7*ones(1,4)];

%1.5m 0.846154 0.846154
%anchor_ids = [ 6*ones(1,6) 7*ones(1,2) 6*ones(1,3) 7*ones(1,3)];

%1m 0.909091
%anchor_ids = [ 6*ones(1,6) 7*ones(1,5) 6*ones(1,6) 7*ones(1,5)];

%1m_v2 0.576923
%anchor_ids = [ 6*ones(1,6) 7*ones(1,6) 6*ones(1,7) 7*ones(1,7)];

%0.3 0.807692
%anchor_ids = [ 6*ones(1,7) 7*ones(1,6) 6*ones(1,7) 7*ones(1,6)];

cfigure(20,12);
[ff,xxi] = ksdensity( anchor_angles * 180/pi);
plot(xxi,ff)
grid on
xlabel('Angular Divergence')
ylabel('Frequency')

%point_start start in x,y,z
%%{
% figure();
% hold on;
% scatter(point_start(:,1), point_start(:,3), 'bo');
% scatter(point_stop(:,1), point_stop(:,3), 'rs');
% scatter(point_start(1,1), point_start(1,3), 'g^');
% 
% for i=1:length(iot.nodes)
%     n = iot.nodes{i};
%     xyz = n.true_p;
%     text(xyz(1), xyz(3)+0.2, num2str(n.id));
%     plot(xyz(1), xyz(3), 'ms','MarkerSize',10,'MarkerFaceColor', 'm');
% end
% xlabel('X');
% ylabel('Z');
% grid on;
% axis equal;
%%}

%% For each pointing segment, gather the range data from each anchor node
SRCID = 9;
WLEN = 0.6;


%% Cooperative detection, minimum overhead
disp('===== Min-cooperative detection =====');
aid_guesses = [];
anchor_ids_final =[];
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
    
    if  size(anc_xyz,1) <=2
        continue;
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
    
    for j=1:length(iot.nodes) %
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
        
        fvec = [dr, r(1), noise, mean(fploss), dAngle];
        
        % classify with RF
        cval = RF.predict(fvec);
        %[c,cval] = svmclassifydist(svmStruct, fvec);

        
        Cvals = [Cvals; cval];
        %Cs = [Cs; c];
        Aids = [Aids; n.id];
        
    end
    
    % collaborative guess
    [~,idx] = max(Cvals); %Cs
    aid_guess = Aids(idx);
    aid_guesses = [aid_guesses; aid_guess];
    anchor_ids_final = [anchor_ids_final;aid];
    
end

nCorrect = sum(aid_guesses == anchor_ids_final);
fprintf('Accuracy: %f\n', nCorrect/length(anchor_ids_final));


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

