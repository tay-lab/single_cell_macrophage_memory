%% Fig 7
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


%%
%%% Ms1=PolyIC, Ms2=CpG
[Ms1,Ms2]=genReceptors(10000);

fig2=figure;
% histogram(Ms1)
rx=[1:100:20000];
histogram(Ms2,30,"Normalization","pdf","FaceColor","r")
hold on
plot(rx,gampdf(rx,2.2,1420),"LineWidth",3,"Color","b")
% plot(rx,gampdf(rx,2.2,1420),"LineWidth",3,"Color","b")
hold off
set(gca,"FontSize",12)
title("TLR 9")
ylabel("Probability")
xlabel("Receptor Level")

%% Mean CpG , sc and Global TNF
nj=500;     % # synthetic cells
Par(9)=1*10^-6;       %sc IFN translation
Par(10)=.005;     %Global IFN production
Par(11)=10^-5;      %sc TNF translation
Par(12)=.5;      %Global TNF translation

%%If flag=1 is Pic, flag=4 is CpG
%%%  clim=[0 2],   clim=[0 1]
flag=1;
clim=[0 2];
[Ms1,Ms2]=genReceptors(1000);
[s,Model]=getScore(Par,Ms1,Ms2,0,mnExp,tt,nj);

fig1=figure;
col=hsv(10);
mod=Model{flag};
Nfkb=mod.NfkB';
norm=max(Nfkb(1,:));

%%
for i=1:3
mod=Model{flag+i};

Nfkb=mod.NfkB';
Nfkb=Nfkb/norm;
a20=mod.A20';
Ikba=mod.Ikba';
IKK=mod.IKK';
auc=trapz(Nfkb(:,1:30),2);
[~,idx]=sort(auc,"descend");
Nfkb=Nfkb(idx,:);
a20=a20(idx,:);
Ikba=Ikba(idx,:);
IKK=IKK(idx,:);

subplot(2,2,1)
hold on
plot(tt,mean(Nfkb),"LineWidth",2)
xticks([0:120:480])
xlabel("Time (min)")
yticks([])
set(gca,"FontSize",12)


end