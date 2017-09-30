function draw = draw_drone( drone )
draw.body=patch('faces',drone.body,'vertices',drone.pts(:,1:end-1));
draw.top=patch('faces',drone.top,'vertices',drone.pts(:,1:end-1));
draw.face=patch('faces',drone.face,'vertices',drone.pts(:,1:end-1));
draw.blade=patch('faces',drone.blade,'vertices',drone.rpts(:,1:end-1));

draw.x=plot3(drone.x(1,:),drone.x(2,:),drone.x(3,:),'linewidth',2,'color','r');
draw.y=plot3(drone.y(1,:),drone.y(2,:),drone.y(3,:),'linewidth',2,'color','g');
draw.z=plot3(drone.z(1,:),drone.z(2,:),drone.z(3,:),'linewidth',2,'color','b');

set(draw.body,'facecolor',[0.1,0.1,0.1],'facealpha',0.8);
set(draw.top,'facecolor','r','facealpha',0.8);
set(draw.face,'facecolor','b','facealpha',0.8);
set(draw.blade,'facecolor',0.7*ones(1,3),'facealpha',0.9,'edgecolor','none');



end

