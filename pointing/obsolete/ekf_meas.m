function [y] = ekf_meas(x,pos,ndx)
	y = x(3+2*(ndx-1)+1) + sqrt((pos(ndx,1) - x(1)).^2 + (pos(ndx,2) - x(2)).^2 + (pos(ndx,3) - x(3)).^2) / 0.299792458;;
end