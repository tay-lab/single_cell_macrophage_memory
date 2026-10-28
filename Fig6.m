%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6 related to 
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
%
% Wang AG, Mercado-Vásquez G, Son M, et al
%
% Corresponding author: tays@uchicago.edu
%
% Single-cell model of polyI:C-induced priming with IFN-beta feedback (Figure 6C-H): autocrine,
% paracrine and mixed IFN-beta signaling; nj = 1000 synthetic cells (receptor levels from genReceptors).
% Uses: functions/getScore.m, TNF_Model_2.m, Params.m, genReceptors.m. Runtime ~40 s (14 workers).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 6
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

[Model,tit]=simOrLoad('Fig6_cells', runSimulation, @() getScore(Ms,pp,tt,nj,dt,tspan,yin,Par,M), [1 2 4 5 6], {'NfkB','TNFR','IFNb','TLR3'});   % 1000 synthetic cells
%%
%%(1) Ctrl->TNF, (2) PIC(TNF-IFN-KO)->Ctrl, (4) PIC(TNF-KO+sIFNaR)->TNF,
%%(5) PIC(TNF-KO+ParaIFN)->TNF, (6) PIC(TNF-KO)->TNF"
close all
figure("Name","Figure 6C-F","Position",[20 20 1600 700])   % wider so the titles do not overlap (was [2000 20 1000 700], off-screen)
flag=1;
clim=[0 2.5];   % set to match the published panel (Fig 6D-F color scale and Fig 6C left axis, 0-2.5)
a2p=[1 2 4 5 6];
% Panel labels as in the published figure (protocol numbers from getScore):
% 1 Ctrl -> TNF, 2 PIC -> Ctrl, 4 autocrine IFN-beta (PIC Auto IFN -> TNF, Fig 6D),
% 5 paracrine IFN-beta (PIC Para IFN -> TNF, Fig 6E), 6 autocrine + paracrine (PIC Para+Auto -> TNF, Fig 6F)
labC=["Ctrl \rightarrow TNF","PIC \rightarrow Ctrl","PIC Auto IFN \rightarrow TNF","PIC Para IFN \rightarrow TNF","PIC Para+Auto \rightarrow TNF"];   % Fig 6C
labTop=["Ctrl \rightarrow TNF","PIC \rightarrow Ctrl","D) Autocrine IFN-\beta","E) Paracrine IFN-\beta","F) Autocrine + Paracrine IFN-\beta"];   % tracemaps
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
title(labTop(i),"FontSize",11)   % CHANGED 2026-10-18: published labels

subplot(2,nl,i+nl)
yyaxis left
plot(tt,mean(Nfkb),"LineWidth",2)
xticks([0:120:480])
xlabel("Time (min)")
if i==1, ylabel("Norm. Nuclear NF-\kappaB"), end   % CHANGED 2026-10-18: label only the outer panels (avoids overlapping labels)
ylim(clim)

yyaxis right
plot(tt,mean(tnfr),"LineWidth",2,"Color",col(1,:))
if i==nl, ylabel("Norm. Total Receptor"), end   % CHANGED 2026-10-18: label only the outer panels (avoids overlapping labels)
ylim([1 3])   % set to match the published panel (Fig 6C right axis, 1-3)
xlim([0 480])
set(gca,"FontSize",12)
if i==1, legend(["NF-\kappaB","TNFR"],"Location","best"); end; title(labC(i),"FontSize",11)   % CHANGED 2026-10-18: published labels

end

%%
Ms=genReceptors(1000);
Ms=prctile(Ms,[10,50,90]);
nj=3;
[Model2]=simOrLoad('Fig6_percentiles', runSimulation, @() getScore(Ms,pp,tt,nj,dt,tspan,yin,Par,M), 6, {'NfkB','TNFR','IFNb','TLR3'});   % 10th/50th/90th percentile cells
tit2=["10th","50th","90th"];

clim=[0 3];
col=hsv(10);
mod=Model2{6};
Nfkb=mod.NfkB';
Nfkb=Nfkb/norm;
tnfr=mod.TNFR'/2000;
Tlr3=mod.TLR3';

fig3=figure("Name","Figure 6G-H");
col2=hsv(3);
for i=1:3
subplot(2,2,i)
plot(tt,Nfkb(i,:),"LineWidth",2)
xticks([0:120:480])
xlabel("Time (min)")
ylabel("Norm. Nuclear NF-\kappaB")
set(gca,"FontSize",12)
ylim(clim)
xlim([0 480])

yyaxis right
plot(tt,tnfr(i,:),"LineWidth",2,"Color",col(1,:))
ylabel("Norm. Total Receptor")
ylim([1 3])   % set to match the published panel (Fig 6G right axis, 1-3)
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
ylim([0 80])

yyaxis right
hold on
zz=tt>240;
plot(tt(zz),Nfkb(i,zz),"LineWidth",2,"LineStyle","--","Color",col2(i,:))
hold off
ylim(clim)
ylabel("Norm. Nuclear NF-\kappaB")

title("TLR3")

end
legend(tit2)

