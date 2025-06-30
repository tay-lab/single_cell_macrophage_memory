%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt=0:480;
% Par=ParFt;
% Par(1)=.1;          %kb R activation
% Par(2)=.001;       %kc R inactivation
Par(3)=5*10^-5;       %ka Adaptor activation 
Par(4)=.01;       %ki Adaptor inactivation 
Par(5)=.0002;       %kap IKK neutral->active 
Par(6)=.0005;       %kip IKK active->inactive 
Par(7)=.002;       %kp IKK inactive->neutral

Par(8)=1*10^-5;       %A20 degradation
Par(9)=5*10^-5;       %sc IFN translation
Par(10)=.5;     %Global IFN production
Par(11)=10^-3;      %sc TNF translation
Par(12)=.5;      %Global TNF translation
Par(13)=.004;      %gamma      

Par(14)= 1.4*10^-7; %kt
Par(15) = .01; %1/s knin
Par(16) = .05; %1/s klin 

Par(17) = 2; %1/s n1
Par(18) = 2; %1/s n2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ms1,Ms2]=genReceptors(1000);
Ms1=prctile(Ms1,[15,25,50,75,95]);
Ms2=prctile(Ms2,[15,25,50,75,95]);
nj=5;
[s,Model]=getScore(Par,Ms1,Ms2,1,[],tt,nj); 
% save ParFit Par tt 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create synthetic cells
tt=0:480;
nj=1000;     % # synthetic cells
[Ms1,Ms2]=genReceptors(nj);
[s,Model,tit]=getScore(Par,Ms1,Ms2,2,[],tt,nj);
%% Figure 5
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


[Ms1,Ms2]=genReceptors(1000);
Ms1=prctile(Ms1,[10,50,90]);
Ms2=prctile(Ms2,[10,50,90]);
nj=3;
[s,Model2]=getScore(Par,Ms1,Ms2,0,[],tt,nj);
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
legend(tit2)

title("TLR 9")

end

%% Figure 6
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

[Ms1,Ms2]=genReceptors(1000);
Ms1=prctile(Ms1,[10,50,90]);
Ms2=prctile(Ms2,[10,50,90]);
nj=3;
[s,Model2]=getScore(Par,Ms1,Ms2,0,[],tt,nj);
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
legend(tit2)

title("TLR 3")

end



%% Figure 7
%%(12) PIC->TNF, (7) PIC(sTNFR)->TNF"
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
tt=0:680;
nj=1000;     % # synthetic cells
[Ms1,Ms2]=genReceptors(nj);
[s,Model,tit]=getScore(Par,Ms1,Ms2,2,[],tt,nj);
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