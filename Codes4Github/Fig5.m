%% Fig 5
% Create synthetic cells
dt=60;
tspan=dt:dt:240*60;
tt=0:480;    % Time course in simulation
nj=1000;     % # synthetic cells
[Ms1,Ms2]=genReceptors(nj);
pp=0;       %%p==0 No plots, p==1 plot traces, p==2 plot heat map
[Model,tit]=getScore(Ms1,Ms2,pp,tt,nj,dt,tspan);
%% Figure 5B-F
%%(14) Ctrl->LPS, (15) LPS(TNF-IFN-KO)->Ctrl, (16) CpG(TNF-IFN-KO)->LPS
flag=14;
clim=[0 1];
a2p=[14 15 16];
nl=length(a2p);
col=hsv(10);
mod=Model{flag};
Nfkb=mod.NfkB';
norm=max(Nfkb(1,:));

for i=1:nl
mod=Model{a2p(i)};
Nfkb=mod.NfkB';
Nfkb=Nfkb/norm;
a20=mod.A20';
Ikba=mod.Ikba';
auc=trapz(Nfkb(:,1:30),2);
[~,idx]=sort(auc,"descend");
Nfkb=Nfkb(idx,:);
a20=a20(idx,:);
Ikba=Ikba(idx,:);

subplot(2,nl,i)
imagesc(tt,1:nj,Nfkb,clim)
xticks([0:120:480])
xlabel("Time (min)")
yticks([])
xlim([0 480])
% xlim([0 240])
set(gca,"FontSize",8)
title(tit(a2p(i)))

subplot(2,nl,i+nl)
yyaxis left
plot(tt,mean(Nfkb),"LineWidth",2)
xticks([0:120:480])
xlabel("Time (min)")
ylabel("Norm. Nuclear NfkB")
ylim(clim)

yyaxis right
plot(tt,mean(a20),"LineWidth",2,"Color",col(1,:))
hold on
plot(tt,mean(Ikba),"LineWidth",2,"Color",col(2,:),"LineStyle","-")
hold off
ylabel("Neg. FB  (uM)")
ylim([0 .04])
xlim([0 480])
set(gca,"FontSize",8)
legend(["NfkB","A20","Ikba"],"Location","best")

end

%%
[Ms1,Ms2]=genReceptors(1000);
Ms1=prctile(Ms1,[10,50,90]);
Ms2=prctile(Ms2,[10,50,90]);
nj=3;
[Model2]=getScore(Ms1,Ms2,0,tt,nj,dt,tspan);
tit2=["10th","50th","90th"];

clim=[0 1.2];
col=hsv(10);
mod=Model2{16};
Nfkb=mod.NfkB';
Nfkb=Nfkb/norm;
a20=mod.A20';
Ikba=mod.Ikba';
Tlr9=mod.TLR9';

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
plot(tt,a20(i,:),"LineWidth",2,"Color",col(1,:))
hold on
plot(tt,Ikba(i,:),"LineWidth",2,"Color",col(2,:),"LineStyle","-")
hold off
ylabel("Neg. FB (uM)")
set(gca,"FontSize",12)
ylim([0 .06])
xlim([0 480])
title(tit2(i))

subplot(2,2,4)
hold on
yyaxis left
plot(tt,Tlr9(i,:),"LineWidth",2,"LineStyle","-","Color",col2(i,:))
hold off
xticks([0:120:480])
xlabel("Time (min)")
ylabel("Active Receptor")
set(gca,"FontSize",12)
ylim([0 150])

yyaxis right
hold on
zz=tt>240;
plot(tt(zz),Nfkb(i,zz),"LineWidth",2,"LineStyle","--","Color",col2(i,:))
hold off
ylim(clim)
ylabel("Norm. Nuclear NFkB")

title("TLR 9")

end
legend(tit2)

