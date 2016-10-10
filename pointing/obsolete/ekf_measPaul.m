function [ z ] = ekf_measPaul( x,t_tx, z_ids,z_pos )

% measurement is of the form [tx, rx1, rx2, ..., rxN]
z = zeros(length(z_ids), 1);

% tx should just be the tx_time
z(1) = t_tx;

% each rx time should be the tx plus the offset and plus the propagation
% delay (ignore first, since that's tx time)
for i=2:length(z_ids)
    id = z_ids(i);
    % get clock offset
    ofst = x( 4 + 2*id );
    % get distance between node and anchor + prop. delay
    xyz_anc = z_pos(i,:)';
    xyz_state = x(1:3);
    dist = norm(xyz_anc - xyz_state);
    prop = dist/0.2998; % m to ns
    
    z(i) = t_tx + ofst + prop;
end


end

