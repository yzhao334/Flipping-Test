function defects = dyntest(Xc,Uc,dynfcn,D,scale)
% Xc = reshape(Xc,12,[]);
% Uc = reshape(Uc,4,[]);
Xc = reshape(Xc,12,numel(Xc)/12);
Uc = reshape(Uc,4,numel(Uc)/4);

defects = (D*(Xc.')/scale).' - dynfcn(0,Xc,Uc);
defects = defects(:);

end

