function [peaktimes_z_p, peaktimes_x_p, locs_z_p,locs_x_p] = getGestureEvents(iot)

%% Load
dp = iot.dataparser;
tstamps = dp.mocaplogs(:,1);
xyz = dp.mocaplogs(:,3:5);

%% Find pointing sessions
% a point is defined as a short motion in a line followed by a retrograde
MINPOINTDIST = 0.5; % m
MINPOINTDWELL = 0.1;  % s
POINTBACKOFF = 1.3;   % s
WINSIZE = 2; % s
WINOVERLAP = 0.40; % s

tstart = tstamps(1);
tstop = tstamps(end);

disp_z = [];
disp_z_times = [];
disp_x = [];
disp_x_times =[];
for t=tstart:(WINSIZE-WINOVERLAP):tstop
    % get the window indices
    idxs = find(tstamps >= t & tstamps < (t+WINSIZE));
    if isempty(idxs)
        continue
    end
    
    % get total displacement in this window
    dx = mean( xyz(idxs(1:end),1)  );;
    dy = mean( xyz(idxs(1:end),2)  );
    dz = mean( xyz(idxs(1:end),3)  );
    
 %   if(dz > dy && dz> dx)
     disp_z = [disp_z; dz];
     disp_z_times = [disp_z_times; t];
 %   else if(dy > dx && dy> dz)
     disp_x = [disp_x; dx];
     disp_x_times = [disp_x_times; t];           
  %      end
  %  end
end

% get positive displacement peaks
[peaks_z_p, locs_z_p] = findpeaks( -1.*disp_z -mean(-1.*disp_z) , 'MinPeakHeight', MINPOINTDIST);
peaktimes_z_p = disp_z_times(locs_z_p);

[peaks_x_p, locs_x_p] = findpeaks(-1.*disp_x -mean(-1.*disp_x), 'MinPeakHeight', MINPOINTDIST);
peaktimes_x_p = disp_x_times(locs_x_p);

