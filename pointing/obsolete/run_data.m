close all; clear; %pkg load statistics;
addpath('utilities');

% Sensor positions
truth =     [ 3.440, 4.263, 0.774];     % golf

% Anchor position
positions = [ 7.111, 2.884, 2.419;      % alpha
              7.094, 7.708, 2.419;      % bravo
              1.606, 7.674, 2.419;      % charlie
              1.614, 2.812, 2.419;      % delta
              5.284, 5.303, 2.419;      % echo
              3.442, 5.312, 2.419 ];    % foxtrot

% Raw data
prefix = '../../../data/golftx_60min/singletx-ntb-golf/ntb-';
files = {
    'alpha', ...       % alpha
    'bravo', ...       % bravo
    'charlie', ...     % charlie
    'delta', ...       % delta
    'echo', ...        % echo
    'foxtrot', ...     % foxtrot
};
limit = 3600;
draw = true;

if draw
  hA = figure; grid; hold on;
  xlabel('Time (seconds)');
  ylabel('RX clock bias (nanoseconds) with respect to TX clock');
  title('Cumulative error walk of Decawave clocks :: 10 Hz, 3600s');
  hR = figure; grid; hold on;
  xlabel('Time (seconds)');
  ylabel('RX clock bias (nanoseconds) with respect to TX clock');
  title('Error perturbations of Decawave clocks :: 10 Hz, 3600s');
  hD = figure; grid; hold on;
  xlabel('Tau (seconds)');
  ylabel('Standard deviation (nanoseconds)');
  title('Allan deviation of Decawave clocks :: 10 Hz, 3600s');
  hP = figure; grid; hold on;
  xlabel('Time (seconds)');
  ylabel('dBm');
  title('First path power of Decawave clocks :: 10 Hz, 3600s');
  hL = figure; grid; hold on;
  xlabel('Time (seconds)');
  ylabel('dBm');
  title('First path loss of Decawave clocks :: 10 Hz, 3600s');
end

% Work out the maximum packet number
mval = 0;
for i = 1:length(files)
    [prefix files{i}]
    a = importdata([prefix files{i}]);
    ts = a(:, 3);
    ts = ts + cumsum([0; diff(ts)<0]) .* 2^8;
    if (max(ts) > mval)
        mval = max(ts);
   end
end
ptx = zeros(mval, length(files));
prx = zeros(mval, length(files));
ppp = zeros(mval, length(files));
ppl = zeros(mval, length(files));

% Build the TX and RX matrices
for i = 1:length(files)

    % Column conversion
    a = importdata([prefix files{i}]);
    ts = a(:, 3);
    tx = a(:, 4);
    rx = a(:, 5);
    pp = a(:, 6);
    pl = a(:, 8);

    % Simple preprocessing
    ts = ts + cumsum([0; diff(ts)<0]) .* 2^8;
    tx = tx + cumsum([0; diff(tx)<0]) .* 2^40;
    rx = rx + cumsum([0; diff(rx)<0]) .* 2^40;
    rx = rx / 63.8976; % to ns
    tx = tx / 63.8976; % to ns  
    ts = ts(tx - tx(1) < limit*1e9);
    rx = rx(tx - tx(1) < limit*1e9); 
    tx = tx(tx - tx(1) < limit*1e9);
    pp = pp(tx - tx(1) < limit*1e9);
    pl = pl(tx - tx(1) < limit*1e9);

    % Save in nanoseconds
    idx = sub2ind(size(ptx), ts, i*ones(length(ts),1));
    ptx(idx) = tx;
    prx(idx) = rx;
    ppp(idx) = pp;
    ppl(idx) = pl;

    if draw    
      f = files{i};
      switch(f)
          case 'alpha'
              c = [1 0 0];
          case 'bravo'
              c = [0 1 0];
          case 'charlie'
              c = [0 0 1];
          case 'delta'
              c = [0 1 1];
          case 'echo'
              c = [1 1 0];
          case 'foxtrot'
              c = [1 0 1];
          case 'golf'
              c = [0 0 0];
      end
      x = tx(2:end) / 1e9;
      y = (diff(rx) ./ diff(tx) - 1) * 1e9;     

      [ad,~,~,tau] = allan(struct('freq',y,'time',x),logspace(-1,3,50),files{i},0);
      rw = ad(tau==1);    % Random walk
      nf = min(ad);       % Noise floor
      disp(['[' files{i} '] TAUM=' num2str(nf) ', TAU1=' num2str(rw)]);

      figure(hD);
      plot(tau,ad,'Color',c,'LineWidth',2);
      figure(hA);
      plot(x,y,'Color',c,'LineWidth',2);
      figure(hR);
      plot(x(2:end),diff(y),'Color',c,'LineWidth',2);   
      figure(hP);
      plot(tx/1e9,pp,'Color',c,'LineWidth',2);   
      figure(hL);
      plot(tx/1e9,pl,'Color',c,'LineWidth',2);   
    end

end

if draw
  figure(hA);
  legend(files);
  figure(hR);
  legend(files);
  figure(hD);
  set(gca,'XScale','log');
  set(gca,'YScale','log');
  legend(files,'Location','northwest');
end

save('data.mat','truth','positions','ptx','prx','ppp','ppl');
