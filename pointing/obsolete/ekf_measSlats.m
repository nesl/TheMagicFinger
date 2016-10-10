function [ y ] = ekf_measSlats( x, idxi, idxj, Trsp0, Trsp1 )

c = 299792458; % m/s

% find state indices for nodes with index i and j
stateSize = 5;
start_i = (idxi-1)*stateSize + 1;
start_j = (idxj-1)*stateSize + 1;

% extract the needed state info
pi = x(start_i:(start_i+2));
pj = x(start_j:(start_j+2));
oi = x(start_i+3);
oj = x(start_j+3);
bi = x(start_i+4);
bj = x(start_j+4);

% predicted measurement y
% format: y = [d_ij; r_ij; B_ij; R_ij];
y = zeros(4,1);

% In an ideal world we'd have a measurement equation of the following form:
%   z = f(norm(pi-pj)) - k
% But, because its difficult to lienarize around a square root, we'll need
% to minimize the following function instead 
%   z = (norm(pi-pj)/c + k).^2
%     = [sum((pi-pj).^2)] / [c.^2] + k.^2 + 2*sum((pi-pj)*k/c

% d_ij = (range + offset)^2
a = sum((pi-pj).^2)/(c.^2); 
b = (oj - oi);
y(1) = a + b.^2 + 2 * sqrt(a) * b;

% r_ij = (range + t_error)^2
a = sum((pj-pi).^2); 
b = ((c/2)*(bj - bi)*Trsp1);
y(2) = a + b.^2 + 2 * sqrt(a) * b;

% B_ij = bias_j - bias_i
y(3) = bj - bi;

% R_ij = (range + t_error)^2
a = sum((pj-pi).^2); 
b = ((c/4)*(bj - bi)*(Trsp1-Trsp0));
y(4) = a + b.^2 + 2 * sqrt(a) * b;


%{
y(1) = sqrt(sum((pj-pi).^2))/c + (oj - oi);
y(2) = sqrt(sum((pj-pi).^2)) + (c/2)*(bj-bi)*Trsp1;
y(3) = sqrt(sum((pj-pi).^2)) + (c/4)*(bj-bi)*(Trsp1-Trsp0);
y(4) = bj - bi;
%}

end

