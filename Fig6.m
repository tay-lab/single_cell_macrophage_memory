%% Fig 6
addpath functions\
% Create synthetic cells
dt=60;
tspan=dt:dt:240*60;
tt=0:480;    % Time course in simulation
nj=1000;     % # synthetic cells
[Ms1,Ms2]=genReceptors(nj);
pp=0;       %%p==0 No plots, p==1 plot traces, p==2 plot heat map
[Model,tit]=getScore(Ms1,Ms2,pp,tt,nj,dt,tspan);
%%
%%(1) Ctrl->TNF, (2) PIC(TNF-IFN-KO)->Ctrl, (4) PIC(TNF-KO+sIFNaR)->TNF,
%%(5) PIC(TNF-KO+ParaIFN)->TNF, (6) PIC(TNF-KO)->TNF"

flag=1;
clim=[0 1.5];
a2p=[1 2 4 5 6];
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

subplot(2,nl,i)
imagesc(tt,1:nj,Nfkb,clim)
xticks([0:120:480])
xlabel("Time (min)")
yticks([])
xlim([0 480])
% xlim([0 240])
set(gca,"FontSize",12)
title(tit(a2p(i)))

subplot(2,nl,i+nl)
yyaxis left
plot(tt,mean(Nfkb),"LineWidth",2)
xticks([0:120:480])
xlabel("Time (min)")
ylabel("Norm. Nuclear NfkB")
ylim(clim)

yyaxis right
plot(tt,mean(tnfr),"LineWidth",2,"Color",col(1,:))
ylabel("Norm. Total Receptor")
ylim([1 4])
xlim([0 480])
set(gca,"FontSize",12)
legend(["NfkB","TNFR"],"Location","best")

end

%%
[Ms1,Ms2]=genReceptors(1000);
Ms1=prctile(Ms1,[10,50,90]);
Ms2=prctile(Ms2,[10,50,90]);
nj=3;

[Model2]=getScore(Ms1,Ms2,0,tt,nj,dt,tspan);
tit2=["10th","50th","90th"];

clim=[0 2];
col=hsv(10);
mod=Model2{6};
Nfkb=mod.NfkB';
Nfkb=Nfkb/norm;
tnfr=mod.TNFR'/2000;
Tlr3=mod.TLR3';

fig3=figure;
col2=hsv(3);
for i=1:3
subplot(2,2,i)
plot(tt,Nfkb(i,:),"LineWidth",2)
xticks([0:120:480])
xlabel("Time (min)")
ylabel("Norm. Nuclear NfkB")
set(gca,"FontSize",12)
ylim(clim)
xlim([0 480])

yyaxis right
plot(tt,tnfr(i,:),"LineWidth",2,"Color",col(1,:))
ylabel("Norm. Total Receptor")
ylim([1 8])
set(gca,"FontSize",12)
xlim([0 480])
title(tit2(i))

subplot(2,2,4)
hold on
yyaxis left
plot(tt,Tlr3(i,:),"LineWidth",2,"LineStyle","-","Color",col2(i,:))
hold off
xticks([0:120:480])
xlabel("Time (min)")
ylabel("Active Receptor")
set(gca,"FontSize",12)
ylim([0 30])

yyaxis right
hold on
zz=tt>240;
plot(tt(zz),Nfkb(i,zz),"LineWidth",2,"LineStyle","--","Color",col2(i,:))
hold off
ylim(clim)
ylabel("Norm. Nuclear NFkB")

title("TLR 3")

end
legend(tit2)

