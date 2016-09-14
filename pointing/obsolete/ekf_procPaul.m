function [ xnew ] = ekf_procPaul( x, dt )

% start from old state
xnew = x;

% new xyz is old xyz, so ignore
% ...

% update clock offsets
for i=4:2:length(x)
    % clock offset and deriv.
    ofst =  x(i);
    d_ofst = x(i+1);
    xnew(i) = ofst + dt*d_ofst;
    % new deriv doesn't change
    % ...
end

end

