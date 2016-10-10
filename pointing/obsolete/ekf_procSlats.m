function [ snew ] = ekf_procSlats( s, dt )

snew = s;

% perform process update over each node separately
stateSize = 5;
for i=1:stateSize:length(s)
    % position does not move
    % s(i ... i+2)
    
    % offset is incremented by bias
    snew(i+3) = s(i+3) + s(i+4)*dt;
    
    % bias doesn't change
end


end

