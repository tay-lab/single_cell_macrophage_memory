%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7 related to 
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
%
% Wang AG, Mercado-Vásquez G, Son M, et al
%
% Corresponding author: tays@uchicago.edu
%
% Quorum amplification by TNF-alpha and IFN-beta secretion after polyI:C (Figure 7B-G).
% Section 1: polyI:C -> TNF with and without TNF secretion / sTNFR (7B-C). Section "Figure 7-II":
% polyI:C with both, one or no secreted cytokine (7D-G). nj = 1000 synthetic cells.
% Uses: functions/getScore.m, TNF_Model_2.m, Params.m, genReceptors.m. Runtime ~1.2 min (14 workers).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 7
addpath functions simulation_results classifier_results data/flow   % run from the main code folder: all paths are relative to the current folder
runSimulation = false;   % false: load the saved 1000-cell results from simulation_results/ (seconds);
                         % true: re-run the simulation (~1-2 min with a parallel pool) and overwrite them
% Create synthetic cells
pp=2;       %%p==0 No plots, p==1 plot traces, p==2 plot heat map
dt=60;
tspan=dt:dt:240*60;
tt=0:480;    % Time course in simulation
nj=1000;     % # synthetic cells   % 1000 synthetic cells, as reported in the figure legends
Ms=genReceptors(nj);
M=2000;
Par(1)=2.5*10^-5;       %sc IFN translation
Par(2)=1*10^-5;          %sc TNF translation
Par(3)=1;               %Total NF-kB
Par(4)=2;               %Total IKK
Par(5)=1;               %Total Ikkk

yin=zeros(23,1);
yin(1) = M;   %Receptors TNF
yin(2) = M;   %Receptors LPS

yin(13) = 0;   %Nuclear NFkB
yin(14) = 0;  %mRNA of target gene (uM)
yin(15) = 0.001;   %Free Ikba (uM)
yin(16) = 0.0001;   %A20 (uM)
yin(17) = 0;    %Ligand TNF
yin(18) = 0;    %Ligand LPS
yin(19) = 0;    %Ligand PolyIC
yin(20) = 0;    %Ligand CpG
yin(21) = 0;    %Ligand IFNg

[Model,tit]=simOrLoad('Fig7_BC', runSimulation, @() getScore(Ms,pp,tt,nj,dt,tspan,yin,Par,M), [1 7 12], {'NfkB','TNFR','IFNb','TNF'});   % Fig 7B-C, 1000 synthetic cells
figure("Name","Figure 7B-C","Position",[20 20 1000 600])   % explicit figure (needed when results are loaded instead of simulated)
%% (12) PIC->TNF, (7) PIC(sTNFR)->TNF"
flag=1;
clim=[0 3];
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
pp=0;
dt=60;
tspan=dt:dt:360*60;
tt=0:360;    % Time course in simulation
nj=1000;     % # synthetic cells   % 1000 synthetic cells, as reported in the figure legends
Ms=genReceptors(nj);
M=2000;

[Model2,tit]=simOrLoad('Fig7_DG', runSimulation, @() getScore(Ms,pp,tt,nj,dt,tspan,yin,Par,M), [3 6 10 13], {'NfkB','TNFR','IFNb','TNF'});   % Fig 7D-G, 1000 synthetic cells
figure("Name","Figure 7D-G","Position",[20 20 1000 600])   % explicit figure (needed when results are loaded instead of simulated)

flag=3;
clim=[0 2];
a2p=[13 3 6 10];
nl=length(a2p);
col=hsv(10);
mod=Model2{flag};
Nfkb=mod.NfkB';
norm=mean(max(Nfkb(:,1:240),[],2));

for i=1:nl
mod=Model2{a2p(i)};
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
ylim([0 2])   % set to match the published panel (Fig 7D-G left axis, 0-2)

yyaxis right
plot(tt,1e3*mean(ifn),"LineWidth",2,"Color",col(1,:))   % uM -> nM (x1000), as in the published panel
hold on
plot(tt,1e3*mean(tnf),"LineWidth",2,"Color",col(2,:),"LineStyle","-")   % uM -> nM (x1000)
hold off

ylabel("Secreted Cytokine (nM)")
ylim([0 .02])   % set to match the published panel (Fig 7D-G right axis, 0-0.02 nM)
xlim([0 360])
set(gca,"FontSize",12)
legend(["NfkB","IFNb","TNF"],"Location","best")

end


