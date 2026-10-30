%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5 related to 
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
%
% Wang AG, Mercado-Vásquez G, Son M, et al
%
% Corresponding author: tays@uchicago.edu
%
% Single-cell model of CpG-induced tolerance (Figure 5B-F).
% Simulates nj = 1000 synthetic cells with initial TLR9 levels drawn from Gamma(2.2, 1420), the fit
% to the TLR9 flow cytometry (Figure 5C), and plots the population means (5B), the single-cell
% tracemap (5D), the 10th/50th/90th-percentile cells (5E) and their active receptor (5F).
% Uses: functions/getScore.m, TNF_Model_2.m, Params.m, genReceptors.m. Runtime ~1.5 min (14 workers).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Load functions
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
% Par=ParS;
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

[Model,tit]=simOrLoad('Fig5_cells', runSimulation, @() getScore(Ms,pp,tt,nj,dt,tspan,yin,Par,M), [14 15 16], {'NfkB','A20','Ikba','TLR9'});   % 1000 synthetic cells
%% Figure 5B-F
%%(14) Ctrl->LPS, (15) LPS(TNF-IFN-KO)->Ctrl, (16) CpG(TNF-IFN-KO)->LPS
figure
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
ylim([0 .02])   % set to match the published panel (Fig 5B right axis, 0-0.02 uM)
xlim([0 480])
set(gca,"FontSize",8)
legend(["NfkB","A20","Ikba"],"Location","best")

end

%%
Ms=genReceptors(1000);
Ms=prctile(Ms,[10,50,90]);
nj=3;
[Model2]=simOrLoad('Fig5_percentiles', runSimulation, @() getScore(Ms,pp,tt,nj,dt,tspan,yin,Par,M), 16, {'NfkB','A20','Ikba','TLR9'});   % 10th/50th/90th percentile cells

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
ylim([0 .02])   % set to match the published panel (Fig 5E right axis, 0-0.02 uM)
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

%% Figure 5C: TLR9 levels in naive BMDMs (flow cytometry) and Gamma fit
% Sample 6 (TLR9-APC, fixed + permeabilized), "Single Cells" gate exported from FlowJo/FACSDiva:
% data/flow/export_Sample_6_TLR9_Fix+Perm_Single Cells.csv (raw .fcs files in data/flow/fcs).
% Gamma fit: fitdist(...,"Gamma") on the positive APC intensities.
flowFile = 'data/flow/export_Sample_6_TLR9_Fix+Perm_Single Cells.csv';   % relative to the main code folder
receptordist = readtable(flowFile);
receptordist = table2array(receptordist(:,7));        % column 7 = 640nm_670/30_APC-A (TLR9)
receptordist = receptordist(receptordist > 0);
pd = fitdist(receptordist,"Gamma");
fprintf("Figure 5C: n = %d cells, Gamma shape a = %.2f, scale b = %.0f\n", numel(receptordist), pd.a, pd.b);
figure("Name","Figure 5C")
histogram(receptordist,'Normalization','pdf','FaceColor',[0.98 0.6 0.6]); hold on
x=0:10:2e4; plot(x,pdf(pd,x),"LineWidth",2,"Color","b"); hold off
xlim([0 2e4]); xlabel("Receptor Level"); ylabel("Probability"); title("TLR9 distribution")
text(0.55*2e4,0.9*max(pdf(pd,x)),sprintf("\\Gamma(%.1f,%.0f)",pd.a,pd.b))
