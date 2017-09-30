function [cost] = costTest_long(Xc,Uc,dist_target,dist_time,tc,ww,scale,u0)
% Xc = reshape(Xc,12,[]);
Xc = reshape(Xc,12,numel(Xc)/12);
Uc = reshape(Uc,4,numel(Uc)/4);
% Uc = reshape(Uc,4,[]);
target = interp1(dist_time(:),dist_target,tc(:));
target = target.';
temp = Xc - target;
% temp = diag([10000*ones(1,3),ones(1,3),...
%     10000*ones(1,3),ones(1,3)])*temp;

% cost=sum(ww(:).'.*(sum(temp.^2,1)+sum(Uc.^2,1)/100))*scale;
cost=sum(ww(:).'.*(sum(temp.^2,1)+sum(Uc.^2,1)/2000^2/100))*scale;
% cost=sum(ww(:).'.*(sum(temp.^2,1)))*scale;
% cost=sum(ww(:).'.*(sum(temp.^2,1)+sum(Uc.^2,1)/2000^2/1000))*scale;
% temp = Uc - u0(:)*ones(1,numel(Uc)/4);
% tercost = sum(sum(Uc(:,[1:3,numel(Uc)/4-10:numel(Uc)/4]).^2))*1000;
tercost2 = sum(sum(Xc(7:12,[floor(numel(Xc)/12/2):numel(Xc)/12]).^2))*10000*(2000^2);
% tercost2 = sum(sum(Xc(7:12,[numel(Xc)/12-4:numel(Xc)/12]).^2))*100*(2000^2);
% tercost2 = sum(sum(Xc(7:12,[numel(Xc)/12-4:numel(Xc)/12]).^2))*100;
% cost=sum( ww(:).'.*(sum(Uc.^2,1)+sum(Xc(7:12,:).^2*2000^2,1)) )*scale;
% tercost = sum(sum(temp(:,[1:3,numel(Uc)/4-5:numel(Uc)/4]).^2))*1000;
% cost=sum(ww(:).'.*(sum(temp.^2,1)))*scale;
% cost = cost+tercost;

cost = cost+tercost2;

% cost = cost+tercost+tercost2;

end

