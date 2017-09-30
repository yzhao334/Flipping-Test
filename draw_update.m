function [ drone ] = draw_update( drone,trans,quad,z )
X = z(1,:);
Y = z(2,:);
Z = z(3,:);
ay = z(4,:);
ax = z(5,:);
az = z(6,:);
% phi = z(4,:);
% the = z(5,:);
% psi = z(6,:);
% R=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1].'...
%     *[cos(the) 0 -sin(the);0 1 0;sin(the) 0 cos(the)].'...
%     *[1 0 0;0 cos(psi) sin(psi);0 -sin(psi) cos(psi)].';
R = [cos(ay) 0 -sin(ay);0 1 0;sin(ay) 0 cos(ay)].'...
    *[1 0 0;0 cos(ax) sin(ax);0 -sin(ax) cos(ax)].'...
    *[cos(az) sin(az) 0;-sin(az) cos(az) 0;0 0 1].';
T=[X;Y;Z];
trans=trans*[R,T;0 0 0 1];


drone.pts=drone.pts0*trans.';
drone.rpts=drone.rpts0*trans.';
drone.x=[trans(1:3,4),trans(1:3,4)+quad.d*trans(1:3,1)];
drone.y=[trans(1:3,4),trans(1:3,4)+quad.d*trans(1:3,2)];
drone.z=[trans(1:3,4),trans(1:3,4)+quad.d*trans(1:3,3)];

end

