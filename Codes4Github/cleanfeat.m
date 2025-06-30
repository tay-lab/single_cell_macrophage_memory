function [fE2,fL2]=cleanfeat(MIfeatE,MIfeatL,fE,fL,datvioE,datvioL,display,bb)

[miL,idxpL]=sort(MIfeatL(:,bb),"descend");
[miE,idxpE]=sort(MIfeatE(:,bb),"descend");
if bb==1
zzE=miE>.65;
else
zzE=miE>.45;
end
idxpE=idxpE(zzE);
miE=miE(zzE);
fE2=fE(idxpE);

if bb==1
zzL=miL>.16;
else
zzL=miL>.02;
end
idxpL=idxpL(zzL);
miL=miL(zzL);
fL2=fL(idxpL);

aa=[];
dat2c=[];
for i=1:length(fE2)
    if anynan(datvioE{idxpE(i),bb,2})
        aa=[aa; i];
    else
    dat2c=[dat2c datvioE{idxpE(i),bb,2}];
    end
end
crE=abs(corr(dat2c));
crE=triu(crE,1);

cc=[];
dat2c=[];
for i=1:length(fL2)
    if anynan(datvioL{idxpL(i),bb,2})
        cc=[cc; i];
    else
    dat2c=[dat2c datvioL{idxpL(i),bb,2}];
    end
end
crL=abs(corr(dat2c));
crL=triu(crL,1);

fE2(aa)=[];
fL2(cc)=[];

if display
    subplot(1,2,1)
h=heatmap(crE, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.XDisplayLabels=fE2;
h.YDisplayLabels=fE2;

subplot(1,2,2)
h=heatmap(crL, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.XDisplayLabels=fL2;
h.YDisplayLabels=fL2;
end

[ii,jj]=find(crE>.95);
if ~isempty(jj)
    fE2(jj)=[];
end

[ii,jj]=find(crL>.95);
if ~isempty(jj)
    fL2(jj)=[];
end

end
