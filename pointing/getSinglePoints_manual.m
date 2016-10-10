function [point_idxs, point_times, point_start, point_stop, anchor_ids, anchor_angles,all_angles] = getSinglePoints_angles(iot)

%% Load
dp = iot.dataparser;
tstamps = dp.mocaplogs(:,1);
xyz = dp.mocaplogs(:,3:5);

%% Find pointing sessions
% a point is defined as a short motion in a line followed by a retrograde
MINPOINTDIST = 0.20; % m
MINPOINTDWELL = 0.1;  % s
POINTBACKOFF = 1.3;   % s
WINSIZE = 0.5; % s
WINOVERLAP = 0.4; % s

tstart = tstamps(1);
tstop = tstamps(end);

displacements_old = [];
disptimes = [];

for t=tstart:(WINSIZE-WINOVERLAP):tstop
    % get the window indices
    idxs = find(tstamps >= t & tstamps < (t+WINSIZE));
    if isempty(idxs)
        continue
    end
    
    % get total displacement in this window
    dx = abs( xyz(idxs(end),1) - xyz(idxs(1),1) );
    dy = abs( xyz(idxs(end),2) - xyz(idxs(1),2) );
    dz = abs( xyz(idxs(end),3) - xyz(idxs(1),3) );
    dxyz = dx + dy + dz;
    
    displacements_old = [displacements_old; dxyz];
    disptimes = [disptimes; t];
end

displacements = (displacements_old);

plot(disptimes-tstart,displacements)
xlabel('Time (sec)');
ylabel('Displacement (m)');

% get positive displacement peaks
[peaks, locs] = findpeaks(displacements, 'MinPeakHeight', MINPOINTDIST);
peaktimes = disptimes(locs);

plot(disptimes , displacements, 'b-');

hold on;
plot(peaktimes, peaks, 'rs');

% detect when the points first began
point_times = [];
point_idxs = [];
GROUPWIN = 3;
GROUPOVR = 2;
t = tstart;
while t < tstop
    t = t + (GROUPWIN-GROUPOVR);
    % how many peaks in window
    pt = find(peaktimes >= t & peaktimes <= t + GROUPWIN);
    if length(pt) >= 2
        time_pstart = peaktimes(pt(1));
        % find closest points in mocap data
        point_idx = dsearchn(tstamps, time_pstart);
        point_times = [point_times; time_pstart];
        point_idxs = [point_idxs; point_idx];
        % advance in time to ignore other points in this window
        t = t + GROUPWIN;
    end
end

plot(point_times, ones(length(point_times)), 'm^');

% pointing start and stop points
point_start = [];
point_stop = [];
for i=1:length(point_idxs)
    idx_start = point_idxs(i);
    point_start = [point_start; xyz(idx_start,:)];
    
    idx_stop = dsearchn(tstamps, point_times(i) + WINSIZE);
    point_stop = [point_stop; xyz(idx_stop,:)];
end

% find closest anchor for each point. Must ensure anchor is in
% forward-pointing direction as well
anchor_ids = [];
anchor_angles = [];
all_angles =[];
for i=1:length(point_idxs)
    xyz_start = point_start(i,:)';
    xyz_stop = point_stop(i,:)';
    closest_angle = inf;
    closest_id = -1;
    angles =[];
    for j=1:length(iot.nodes)
        node = iot.nodes{j};
        xyz = node.true_p;
        % if the node is mobile (not an anchor) ignore it
        if node.isMobile()
            continue;
        end
        % get angle between this node and the pointing vector
        v1 = (xyz_stop-xyz_start)./norm(xyz_stop-xyz_start);
        v2 = (xyz-xyz_stop)./norm(xyz-xyz_stop);
        angle = atan2(norm(cross(v1,v2)), dot(v1,v2));
        if (angle < closest_angle )%&& node.id ~= 6 && node.id ~= 7)
            closest_angle = angle;
            closest_id = node.id;
        end
        angles=[angles,angle];
    end
    anchor_ids = [anchor_ids; closest_id];
    anchor_angles = [anchor_angles; closest_angle];
    all_angles =[ all_angles; angles];
end
