%% Fig 7
addpath functions\
% Create synthetic cells
dt=60;
tspan=dt:dt:240*60;
tt=0:480;    % Time course in simulation
nj=1000;     % # synthetic cells
[Ms1,Ms2]=genReceptors(nj);
pp=0;       %%p==0 No plots, p==1 plot traces, p==2 plot heat map
[Model,tit]=getScore(Ms1,Ms2,pp,tt,nj,dt,tspan);
%% (12) PIC->TNF, (7) PIC(sTNFR)->TNF"
flag=1;
clim=[0 1.5];
a2p=[12 7];
nl=length(a2p);
col=hsv(10);
mod=Model{flag};
Nfkb=mod.NfkB';
norm=max(Nfkb(1,:));

for i=1:nl
mod=Model{a2p(i)};
Nfkb=mod.NfkB';
Nfkb=Nfkb/norm;
tnfr=mod.TNFR'/2000;
ifn=mod.IFNb';
auc=trapz(Nfkb(:,1:30),2);
[~,idx]=sort(auc,"descend");
Nfkb=Nfkb(idx,:);
tnfr=tnfr(idx,:);
ifn=ifn(idx,:);

subplot(1,nl,i)
imagesc(tt,1:nj,Nfkb,clim)
xticks([0:120:480])
xlabel("Time (min)")
yticks([])
xlim([0 480])
% xlim([0 240])
set(gca,"FontSize",8)
title(tit(a2p(i)))

end

%% Figure 7-II
%%(13) PIC->Ctrl", (3) PIC(TNF-IFN-KO)->Ctrl", (6) PIC(TNF-KO)->TNF", (10) PIC(IFNa-KO)->TNF"
% Create synthetic cells
dt=60;
tspan=dt:dt:360*60;
tt=0:360;    % Time course in simulation
nj=1000;     % # synthetic cells
[Ms1,Ms2]=genReceptors(nj);
pp=0;       %%p==0 No plots, p==1 plot traces, p==2 plot heat map
[Model,tit]=getScore(Ms1,Ms2,pp,tt,nj,dt,tspan);
%%
flag=3;
clim=[0 2];
a2p=[13 3 6 10];
nl=length(a2p);
col=hsv(10);
mod=Model{flag};
Nfkb=mod.NfkB';
norm=mean(max(Nfkb(:,1:240),[],2));

for i=1:nl
mod=Model{a2p(i)};
Nfkb=mod.NfkB';
Nfkb=Nfkb/norm;
tnfr=mod.TNFR'/2000;
ifn=mod.IFNb';
tnf=mod.TNF';
auc=trapz(Nfkb(:,1:30),2);
[~,idx]=sort(auc,"descend");
Nfkb=Nfkb(idx,:);
tnf=tnf(idx,:);
ifn=ifn(idx,:);

subplot(2,nl,i)
imagesc(tt,1:nj,Nfkb,clim)
xticks([0:120:480])
xlabel("Time (min)")
yticks([])
xlim([0 360])
% xlim([0 240])
set(gca,"FontSize",12)
title(tit(a2p(i)))

subplot(2,nl,i+nl)
yyaxis left
plot(tt,mean(Nfkb),"LineWidth",2)
xticks([0:120:480])
xlabel("Time (min)")
ylabel("Norm. Nuclear NfkB")
ylim([0 1.5])

yyaxis right
plot(tt,mean(ifn),"LineWidth",2,"Color",col(1,:))
hold on
plot(tt,mean(tnf),"LineWidth",2,"Color",col(2,:),"LineStyle","-")
hold off

ylabel("Secreted Cytokine (uM)")
ylim([0 1.2e-4])
xlim([0 360])
set(gca,"FontSize",12)
legend(["NfkB","IFNb","TNF"],"Location","best")

end


