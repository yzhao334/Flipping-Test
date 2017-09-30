function dz = droneDynamics(z,u,quadpar)
% drone dynamics -- customized version
% dz = droneDynamics(z,u,p)
%
% This function computes the first-order dynamics of the quadrator/drone.
%
% INPUTS:
%   x = [12, n] = [P;O;dp;do] = state of the system
%           P = [x;y;z], position in inertia frame
%           O = [ay;ax;az], euler angles
%           dp = [dx;dy;dz], vel in inertia frame
%           do = [wx;wy;wz], anguler velocity
%   u = [4, n] =  control to the drone
%   quad = parameter struct
% orientation: y-x-z angle
% angles: ay, ax, az
%       
% OUTPUTS:
%   dz = dz/dt = time derivative of state
%
% 
% R = 
% [ cos(ay)*cos(az) + sin(ax)*sin(ay)*sin(az), cos(az)*sin(ax)*sin(ay) - cos(ay)*sin(az), cos(ax)*sin(ay)]
% [                           cos(ax)*sin(az),                           cos(ax)*cos(az),        -sin(ax)]
% [ cos(ay)*sin(ax)*sin(az) - cos(az)*sin(ay), sin(ay)*sin(az) + cos(ay)*cos(az)*sin(ax), cos(ax)*cos(ay)]
% 
% X = z(1,:);% x pos
% Y = z(2,:);% y pos
% Z = z(3,:);% z pos
ay = z(4,:);% x vel
ax = z(5,:);% y vel
az = z(6,:);% z vel
dx = z(7,:);
dy = z(8,:);
dz = z(9,:);
wx = z(10,:);
wy = z(11,:);
wz = z(12,:);

w1=u(1,:);
w2=u(2,:);
w3=u(3,:);
w4=u(4,:);

ct=quadpar.Ct*quadpar.rho*quadpar.A*quadpar.r^2;
cQ=quadpar.Cq*quadpar.rho*quadpar.A*quadpar.r^3;
d=sqrt(2)/2*quadpar.d;
M =...
[  -ct,   -ct,   -ct,   -ct;...
  ct*d, -ct*d, -ct*d,  ct*d;...
  ct*d,  ct*d, -ct*d, -ct*d;...
   -cQ,    cQ,   -cQ,    cQ];

Trqs = M*[w1.^2;w2.^2;w3.^2;w4.^2];

Tz=Trqs(1,:); % force, body frame
Tay=Trqs(2,:); % torque 1, body frame
Tax=Trqs(3,:); % torque 2, body frame
Taz=Trqs(4,:); % torque 3, body frame

% in inertia frame
% tau = 
%  Tay*(cos(ay)*cos(az) + sin(ax)*sin(ay)*sin(az)) - Tax*(cos(ay)*sin(az) - cos(az)*sin(ax)*sin(ay)) + Taz*cos(ax)*sin(ay)
%                                                                  Tax*cos(ax)*cos(az) - Taz*sin(ax) + Tay*cos(ax)*sin(az)
%  Tax*(sin(ay)*sin(az) + cos(ay)*cos(az)*sin(ax)) - Tay*(cos(az)*sin(ay) - cos(ay)*sin(ax)*sin(az)) + Taz*cos(ax)*cos(ay)
% tz = 
%  cos(ax)*sin(ay)*Tz
%         -sin(ax)*Tz
%  cos(ax)*cos(ay)*Tz
 

ddx = 1/quadpar.M*(cos(ax).*sin(ay).*Tz);
ddy = 1/quadpar.M*(-sin(ax).*Tz);
ddz = 1/quadpar.M*(cos(ax).*cos(ay).*Tz)+quadpar.g;

day = (wy.*cos(ax) + wz.*cos(ay).*sin(ax) + wx.*sin(ax).*sin(ay))./cos(ax);
dax = wx.*cos(ay) - wz.*sin(ay);
daz = (wz.*cos(ay) + wx.*sin(ay))./cos(ax);

J1 = quadpar.J(1,1);
J2 = quadpar.J(2,2);
J3 = quadpar.J(3,3);

% -Omega x I Omega
%                                                   -(dax*sin(ay) - daz*cos(ax)*cos(ay))*(J2 - J3)*(day - daz*sin(ax))
%  (J1 - J3)*(dax^2*cos(ay)*sin(ay) + dax*daz*cos(ax) - 2*dax*daz*cos(ax)*cos(ay)^2 - daz^2*cos(ax)^2*cos(ay)*sin(ay))
%                                                    (dax*cos(ay) + daz*cos(ax)*sin(ay))*(J1 - J2)*(day - daz*sin(ax))

dw =   [-(dax.*sin(ay) - daz.*cos(ax).*cos(ay)).*(day - daz.*sin(ax))*(J2 - J3) + ...
        Tay.*(cos(ay).*cos(az) + sin(ax).*sin(ay).*sin(az)) - Tax.*(cos(ay).*sin(az) - cos(az).*sin(ax).*sin(ay)) + Taz.*cos(ax).*sin(ay);...
        (J1 - J3)*(dax.^2.*cos(ay).*sin(ay) + dax.*daz.*cos(ax) - 2*dax.*daz.*cos(ax).*cos(ay).^2 - daz.^2.*cos(ax).^2.*cos(ay).*sin(ay)) + ...
        Tax.*cos(ax).*cos(az) - Taz.*sin(ax) + Tay.*cos(ax).*sin(az);...
        (J1 - J2)*(dax.*cos(ay) + daz.*cos(ax).*sin(ay)).*(day - daz.*sin(ax)) + ...
        Tax.*(sin(ay).*sin(az) + cos(ay).*cos(az).*sin(ax)) - Tay.*(cos(az).*sin(ay) - cos(ay).*sin(ax).*sin(az)) + Taz.*cos(ax).*cos(ay)];
% dw = quadpar.J\dw; % dw = [dwx;dwy;dwz]
dwx = dw(1,:)/J1;
dwy = dw(2,:)/J2;
dwz = dw(3,:)/J3;




dz = [dx;dy;dz;...
      day;dax;daz;...
      ddx;ddy;ddz;...
      dwx;dwy;dwz];

% [dX,dY,dZ,dphi,dthe,dpsi,ddx,ddy,ddz,dp,dq,dr] = autoGen_droneDynamics(phi,the,psi,dx,dy,dz,p,q,r, w1,w2,w3,w4, quadpar);

% dz = [dX;dY;dZ;dphi;dthe;dpsi;ddx;ddy;ddz;dp;dq;dr];

end