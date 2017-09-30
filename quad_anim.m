function quad_anim( z,quad,freq )
% animation
clf;pause(0.1);
trans=[[1 0 0;0 -1 0; 0 0 -1].',[0;0;0]; 0 0 0 1];
T=size(z,2)-1;% total horizon
d=quad.d/sqrt(2);
drone.N=20;% blade discritization
% faces
drone.top=[1 2 3 4];
drone.body=[5 6 7 8;2 6 7 3;3 7 8 4;4 8 5 1];
drone.face=[1 5 6 2];
for i=1:4
    drone.blade(i,:)=(i-1)*drone.N+1:i*drone.N;
end
% vertices
drone.pts0=[d -d quad.h;d d quad.h; -d d quad.h;-d -d quad.h;...
           d -d 0;     d d 0;      -d d 0;     -d -d 0];
drone.pts0=[drone.pts0,ones(8,1)];
drone.rpts0=[ones(drone.N,1)*(drone.pts0(1,1:3)+[0,0,-0.005])+quad.r*[cos(linspace(0,-2*pi,drone.N)).',sin(linspace(0,-2*pi,drone.N)).',zeros(drone.N,1)];...
            ones(drone.N,1)*(drone.pts0(2,1:3)+[0,0,-0.005])+quad.r*[cos(linspace(0,-2*pi,drone.N)).',sin(linspace(0,-2*pi,drone.N)).',zeros(drone.N,1)];...
            ones(drone.N,1)*(drone.pts0(3,1:3)+[0,0,-0.005])+quad.r*[cos(linspace(0,-2*pi,drone.N)).',sin(linspace(0,-2*pi,drone.N)).',zeros(drone.N,1)];...
            ones(drone.N,1)*(drone.pts0(4,1:3)+[0,0,-0.005])+quad.r*[cos(linspace(0,-2*pi,drone.N)).',sin(linspace(0,-2*pi,drone.N)).',zeros(drone.N,1)]];
drone.rpts0=[drone.rpts0,ones(4*drone.N,1)];        
% update data
drone=draw_update(drone,trans,quad,z(:,1));
% initialize drawing
hold on;
view([37.5,30]);
axis equal vis3d;box on;grid on;
min_p=min(z(1:3,:),[],2);
max_p=max(z(1:3,:),[],2);
xlim([min_p(1)-quad.d-quad.r,max_p(1)+quad.d+quad.r]);
ylim([-(max_p(2)+quad.d+quad.r),-(min_p(2)-quad.d-quad.r)]);
zlim([-(max_p(3)+quad.d+quad.r),-(min_p(3)-quad.d-quad.r)]);
xlabel('X');ylabel('Y');zlabel('Z');


draw=draw_drone(drone);

for t=0:T        
    drone=draw_update(drone,trans,quad,z(:,t+1));
    set(draw.body,'vertices',drone.pts(:,1:end-1));
    set(draw.top,'vertices',drone.pts(:,1:end-1));
    set(draw.face,'vertices',drone.pts(:,1:end-1));
    set(draw.blade,'vertices',drone.rpts(:,1:end-1));
    set(draw.x,'xdata',drone.x(1,:),'ydata',drone.x(2,:),'zdata',drone.x(3,:));
    set(draw.y,'xdata',drone.y(1,:),'ydata',drone.y(2,:),'zdata',drone.y(3,:));
    set(draw.z,'xdata',drone.z(1,:),'ydata',drone.z(2,:),'zdata',drone.z(3,:));
    pause(1/freq); % show slow animation
end


end

