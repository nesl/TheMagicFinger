%function [point_idxs, point_times, point_start, point_stop, anchor_ids, anchor_angles] = getPointEvents(iot)

%% Load
dp = iot.dataparser;
tstamps = dp.mocaplogs(:,1);
xyz = dp.mocaplogs(:,3:5);

%% Find pointing sessions
% a point is defined as a short motion in a line followed by a retrograde
MINPOINTDIST = 0.60; % m
MINPOINTDWELL = 0.1;  % s
POINTBACKOFF = 1.3;   % s
WINSIZE = 0.5; % s
WINOVERLAP = 0.40; % s

tstart = tstamps(1);
tstop = tstamps(end);

displacements = [];
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
    
    displacements = [displacements; dxyz];
    disptimes = [disptimes; t];
end

plot(disptimes-tstart,displacements)
xlabel('Time (sec)');
ylabel('Displacement (m)');

% get positive displacement peaks
[peaks, locs] = findpeaks(displacements, 'MinPeakHeight', MINPOINTDIST);
peaktimes = disptimes(locs);

% remove retrogrades using backoff timer
retro_idxs = [];
for i=1:length(peaks)
    % if retro point, skip
    if ~isempty(find(retro_idxs == i))
        continue;
    end
    
    t = disptimes(locs(i));
    tretro = t + POINTBACKOFF;
    % find retro points
    idxs = find(peaktimes > t & peaktimes <= tretro);
    % if there IS no retro point, schedule this point for removal as well
    if isempty(idxs)
        idxs = i;
    end
    retro_idxs = [retro_idxs; idxs];
end
forward_peaks = peaks;
forward_locs = locs;
forward_peaks(retro_idxs) = [];
forward_locs(retro_idxs) = [];
forward_times = disptimes(forward_locs);

hold on;
plot(forward_times-tstart, forward_peaks, 'gs');

% find groups of 5 as he points five times
GROUPTIMEOUT = 2.5; % s
point_times = [];

i = 1;
while i < length(forward_times)
    % create candidate group
    group = [forward_times(i)];
    
    % time to next point
    dt = forward_times(i+1) - forward_times(i);
    while dt < GROUPTIMEOUT
        group = [group; forward_times(i+1)];
        i = i+1;
        if i == length(forward_times)
            break;
        end
        dt = forward_times(i+1) - forward_times(i);
    end
    
    % append if group of 5
    if length(group) == 5
        point_times = [point_times; group];
    end
    
    i = i+1;
end

plot(point_times-tstart, ones(size(point_times)), 'ro');

% find closest points in mocap data
point_idxs = dsearchn(tstamps, point_times);

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

for i=1:length(point_idxs)
    xyz_start = point_start(i,:)';
    xyz_stop = point_stop(i,:)';
    closest_angle = inf;
    closest_id = -1;
    for j=1:length(iot.nodes)
        node = iot.nodes{j};
        xyz = node.true_p;
        % if the node is mobile (not an anchor) ignore it
        if node.isMobile()
            continue;
        end
        % get angle between this node and the pointing vector
        v1 = (xyz_stop-xyz_start)./norm(xyz_stop-xyz_start);
        v2 = (xyz-xyz_start)./norm(xyz-xyz_start);
        angle = atan2(norm(cross(v1,v2)), dot(v1,v2));
        if angle < closest_angle
            closest_angle = angle;
            closest_id = node.id;
        end
    end
    anchor_ids = [anchor_ids; closest_id];
    anchor_angles = [anchor_angles; closest_angle];
end

%%{


%%}

