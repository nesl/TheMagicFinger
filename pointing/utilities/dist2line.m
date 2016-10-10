function d = dist2line(pt, v1, v2)
% d = dist2line(pt, v1, v2)
%
% pt = [x;y;z]
% v1 = [x1;y1;z1] on line
% v2 = [x2;y2;z2] on line
%
% d is closest distance between pt to line

      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end