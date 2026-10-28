%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S6 related to
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
% Wang AG, Mercado-Vásquez G, Son M, et al. Corresponding author: tays@uchicago.edu
%
% Initial protein expression variability and Jensen-Shannon distance (Figure S6A-C).
% Experimental data: phase II / phase I NF-kB peak ratios of single cells (genExtFeat.mat; polyI:C -> TNF
% n = 1795 cells, CpG -> LPS n = 973 cells). For each protein (TLR9/TLR3, IkBa, A20/negative feedback,
% total NF-kB, total IKK, total IKKK/IRF/STAT, TLR4/TNFR) the initial levels of 1000 synthetic cells are
% drawn from a Gamma distribution with the mean fixed to the model value (shape = mean/theta,
% scale = theta); theta is scanned and the simulated peak-ratio distribution is compared with the
% experimental one using the Jensen-Shannon distance (log2; 0 = identical, 1 = maximally different) on
% histograms with Scott's-rule bins (1st-99th percentiles of the pooled data).
% Uses: functions/getScoreCondTot.m, TNF_Model_2.m, Params.m, genExtFeat.mat. Long run
% (7 scans x ~20 theta values x 1000 cells, parfor).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions simulation_results classifier_results data/flow   % run from the main code folder: all paths are relative to the current folder
runSimulation = false;   % false: load the saved scans from simulation_results/FigS6_JSD.mat (seconds);
                         % true: re-run the 7 scans (~7 min with a parallel pool, 1000 cells per theta) and overwrite
resFile = 'simulation_results/FigS6_JSD.mat';   % relative to the main code folder
if ~runSimulation && isfile(resFile), SR = load(resFile); else, runSimulation = true; end
load genExtFeat.mat data0
mx=prctile(max([data0.("Phase1 NfkB") data0.("Phase2 NfkB")]),50);
myExp=[];
myExp2=[];
nfkbExp={};
rpkE=[];

cats=["Control -> TNFnaive","Control -> LPSnaive","PolyIC -> TNFcond","CpG -> LPScond"];
rpkED={};
for i=1:2
zz=data0.Category==cats(i);
myExp(i,1)=median(data0.("Phase2 Fourier")(zz));
myExp(i,2)=median(data0.("Phase2 Peak")(zz)/mx);
myExp(i,3)=median(data0.("Phase2 Time to Peak")(zz));
myExp(i,4)=median(data0.("Phase2 Speed")(zz));
rpkED{i}=data0.("Phase2 Peak")(zz)/mx;


nfkbExp{i}=data0.("Phase2 NfkB")(zz,:)/mx;
myExp2(i,:)=median(nfkbExp{i});%-min(mean(nfkbExp));
end

for i=3:4
zz=data0.Category==cats(i);
myExp(i,1)=median(data0.("Phase1 Fourier")(zz));
myExp(i,2)=median(data0.("Phase1 Peak")(zz)/mx);
myExp(i,3)=median(data0.("Phase1 Time to Peak")(zz));
myExp(i,4)=median(data0.("Phase1 Speed")(zz));
rpkED{i}=data0.("Phase1 Peak")(zz)/mx;

nfkbExp{i}=data0.("Phase1 NfkB")(zz,:)/mx;
myExp2(i,:)=median(nfkbExp{i});%-min(mean(nfkbExp));
end

zz=data0.Category=="PolyIC -> TNFcond";
nfkbExp{5}=[data0.("Phase1 NfkB")(zz,:) data0.("Phase2 NfkB")(zz,:)]/mx;
myExp3(5,:)=median(nfkbExp{5});%-min(mean(nfkbExp));
pk1=data0.("Phase1 Peak")(zz)/mx;
pk2=data0.("Phase2 Peak")(zz)/mx;
tk1=data0.("Phase1 Time to Peak")(zz);
tk2=data0.("Phase2 Time to Peak")(zz);
% yy=pk1<=0 | pk2<0 ;
rpkED{5}=[pk1, pk2];

pk=rpkED{5};
rpkEDr{5}=pk(:,2)./pk(:,1);
rpkE(1,:)=[median(pk(:,1)) median(rpkEDr{5}) median(rpkED{1}) median(tk1) median(tk2)];


zz=data0.Category=="CpG -> LPScond";
nfkbExp{6}=[data0.("Phase1 NfkB")(zz,:) data0.("Phase2 NfkB")(zz,:)]/mx;
myExp3(6,:)=median(nfkbExp{6});%-min(mean(nfkbExp));
pk1=data0.("Phase1 Peak")(zz)/mx;
pk2=data0.("Phase2 Peak")(zz)/mx;
tk1=data0.("Phase1 Time to Peak")(zz);
tk2=data0.("Phase2 Time to Peak")(zz);
% yy=pk1<=0 | pk2<0;
rpkED{6}=[pk1, pk2];

pk=rpkED{6};
rpkEDr{6}=pk(:,2)./pk(:,1);
rpkE(2,:)=[median(pk(:,1)) median(rpkEDr{6}) median(rpkED{2}) median(tk1) median(tk2)];

Par(1)=2.5*10^-5;       %sc IFN translation
Par(2)=1*10^-5;          %sc TNF translation
%% Shannon Distance for TLR 9 variability  -> Figure S6A (Gamma distributions) and S6B (JS distance) panel
rng(0)
tt=0:6:474;
nfkb=[];
mm=1000;
Ms=2000;
[~, edges] = histcounts(rmoutliers([rpkEDr{5}' rpkEDr{6}'],"percentiles",[1 99]), 'BinMethod', 'scott');

P=[];
P(:,1)= get_distribution(rpkEDr{5}, edges);
P(:,2)= get_distribution(rpkEDr{6}, edges);

Q=[];
d1=[];
ksc1=10.^[0:5/20:5];

Par=[];
% Par(1)=2.5*10^-5;       %sc IFN translation
% Par(2)=1*10^-5;          %sc TNF translation
Par(3)=1;               %Total NF-kB
Par(4)=2;               %Total IKK
Par(5)=1;               %Total Ikkk

if runSimulation   % simulate 1000 cells for each theta (skipped when loading saved results)
for k=1:length(ksc1)
    k
ssM=ksc1(k);

mnM=3000;
Ms2=gamrnd(mnM/ssM,ssM,[mm 1]);
parfor i=1:mm
    yin=zeros(23,1);
    yin(1) = Ms;   %Receptors TNF
    yin(2) = Ms;   %Receptors LPS
    yin(13) = 0;   %Nuclear NFkB
    yin(15) = 0.001;   %Free Ikba (uM)
    yin(16) = 0.0001;   %A20 (uM)

[~,Nf]=getScoreCondTot(Par,0,rpkE,myExp3,Ms2(i),yin);
for j=1:2
nfkb(:,i,j)=Nf(j,:);
end
end

for j=1:2
for i=1:mm
    fpeak1=max(nfkb(tt<200,i,j));
    fpeak2=max(nfkb(tt>240,i,j));
    rpkNr(i,j)=fpeak2./fpeak1;
end
end

rpkNDr{1}=rpkNr(:,1);
rpkNDr{2}=rpkNr(:,2);

for j=1:2
Q(:,j,k)= get_distribution(rpkNDr{j}, edges);
d1(j,k) = jsd(P(:,j), Q(:,j,k));

end
end

    Q1 = Q;   % keep the TLR9/TLR3 distributions for Figure S6C
else   % load the saved scan
    d1=SR.d1; Q=SR.Q1;
end
%%
fig=figure("Position",[20 100 800 800]);
subplot(1,2,1)
lnd=[];
for k=1:4:length(ksc1)-3
plot(edges(2:end),Q(:,1,k),"LineWidth",2)
xlim([0 5])
hold on
lnd=[lnd ["\theta="+num2str(ksc1(k))]];
end
bb=bar(edges(2:end),P(:,1),"LineWidth",.5)
uistack(bb, 'bottom');
lnd=["Experimental" lnd];
ylabel("Probability density")
xlabel("Peak in phase 2/Peak in phase 1")
axis square
legend(lnd)
set(gca,"FontSize",10)
title("PIC->TNF")

subplot(1,2,2)
for k=1:4:length(ksc1)-3
plot(edges(2:end),Q(:,2,k),"LineWidth",2)
xlim([0 5])
hold on
end
bb=bar(edges(2:end),P(:,2),"LineWidth",.5)
uistack(bb, 'bottom');
axis square
legend(lnd)
ylabel("Probability density")
xlabel("Peak in phase 2/Peak in phase 1")
set(gca,"FontSize",10)
title("CpG->LPS")   % was "CpG->TNF" (typo; this panel is CpG -> LPS)

%%
fig=figure("Position",[20 100 1400 800]);

subplot(2,3,1)
lnd=[];
for k=1:4:length(ksc1)-3
ssM=ksc1(k);
mnM=3000;
xx=[0:10:8000];
    plot(xx,gampdf(xx,mnM/ssM,ssM),"LineWidth",2)
    xlim([0 6000])
    ylim([0 .005])
ylabel("Probability density")
xlabel("Receptor expression")
hold on
lnd=[lnd ["\theta="+num2str(ksc1(k))]];
end
axis square
legend(lnd)
set(gca,"FontSize",10)

subplot(2,3,4)
plot(ksc1, d1',"LineWidth",2)
xlabel("Scale parameter")
ylabel("Jensen-Shannon distance")
ylim([0 1])
xlim([1 100000])
xticks([10,1000,100000])
legend({"PIC->TNF","CpG->LPS"})
axis square
set(gca,"XScale","log","FontSize",10)
title("Initial TLR9/TLR3 variability")

%% Shannon Distance for Ikba variability  -> Figure S6A (Gamma distributions) and S6B (JS distance) panel
rng(0)
tt=0:6:474;
nfkb=[];
mm=1000;
mnM=3000;
Ms=2000;
[~, edges] = histcounts(rmoutliers([rpkEDr{5}' rpkEDr{6}'],"percentiles",[1 99]), 'BinMethod', 'scott');

P=[];
P(:,1)= get_distribution(rpkEDr{5}, edges);
P(:,2)= get_distribution(rpkEDr{6}, edges);

Q=[];
d2=[];
ikbDis=[];

ksc2=10.^[-5:5/30:-2];
if runSimulation   % simulate 1000 cells for each theta (skipped when loading saved results)
for k=1:length(ksc2)
    k
% Ms=gamrnd(mnM/ssM,ssM,[mm 1]);   %TNF Receptor (Measured)

ssM=ksc2(k);
yinikb=0.001;
Ms3=gamrnd(yinikb/ssM,ssM,[mm 1]);
parfor i=1:mm
    yin=zeros(23,1);
    yin(1) = Ms;   %Receptors TNF
    yin(2) = Ms;   %Receptors LPS
    yin(13) = 0;   %Nuclear NFkB
    yin(15) = Ms3(i);   %Free Ikba (uM)
    yin(16) = 0.0001;   %A20 (uM)

[~,Nf]=getScoreCondTot(Par,0,rpkE,myExp3,mnM,yin);
for j=1:2
nfkb(:,i,j)=Nf(j,:);
end
end

for j=1:2
for i=1:mm
    fpeak1=max(nfkb(tt<200,i,j));
    fpeak2=max(nfkb(tt>240,i,j));
    rpkNr(i,j)=fpeak2./fpeak1;
end
end

rpkNDr{1}=rpkNr(:,1);
rpkNDr{2}=rpkNr(:,2);

for j=1:2
Q(:,j,k)= get_distribution(rpkNDr{j}, edges);
d2(j,k) = jsd(P(:,j), Q(:,j,k));

end
end

else   % load the saved scan
    d2=SR.d2;
end
%%
subplot(2,3,2)
lnd=[];
ksc2=10.^[-5:5/30:-2];
for k=1:4:length(ksc2)
ssM=ksc2(k);
yinikb=0.001;
xx=[.0001:.00001:.004];
    plot(xx,gampdf(xx,yinikb/ssM,ssM),"LineWidth",2)
    % xlim([0.2 1])
ylabel("Probability density")
xlabel("Initial Ikba expression (\muM)")
hold on
lnd=[lnd ["\theta="+num2str(ksc2(k),'%.2g')]];
end
axis square
hold off
legend(lnd)
set(gca,"FontSize",10)

subplot(2,3,5)
plot(ksc2, d2',"LineWidth",2)
xlabel("Scale parameter")
ylabel("Jensen-Shannon distance")
ylim([0 1])
% xlim([10 100000])
% xticks([10,100,1000,10000,100000])
legend({"PIC->TNF","CpG->LPS"})
axis square
set(gca,"XScale","log","FontSize",10)
title("Initial IkBa variability")

%% Shannon Distance for A20 variability  -> Figure S6A (Gamma distributions) and S6B (JS distance) panel
rng(0)
tt=0:6:474;
nfkb=[];
mm=1000;
mnM=3000;
[~, edges] = histcounts(rmoutliers([rpkEDr{5}' rpkEDr{6}'],"percentiles",[1 99]), 'BinMethod', 'scott');

P=[];
P(:,1)= get_distribution(rpkEDr{5}, edges);
P(:,2)= get_distribution(rpkEDr{6}, edges);

Q=[];
d3=[];
a20Dis=[];

Ms=2000;
ksc3=10.^[-6:5/20:-2];   % CHANGED 2026-10-18: theta range 1e-6 to 1e-2, as in the published Figure S6
if runSimulation   % simulate 1000 cells for each theta (skipped when loading saved results)
for k=1:length(ksc3)
    k
% Ms=gamrnd(mnM/ssM,ssM,[mm 1]);   %TNF Receptor (Measured)

ssM=ksc3(k);
yina20=0.0001;   % CHANGED 2026-10-18: mean initial negative feedback = 1e-4 uM (Table 1); 1e-3 did not reproduce Figure S6B
Ms3=gamrnd(yina20/ssM,ssM,[mm 1]);
parfor i=1:mm
    yin=zeros(23,1);
    yin(1) = Ms;   %Receptors TNF
    yin(2) = Ms;   %Receptors LPS
    yin(13) = 0;   %Nuclear NFkB
    yin(15) = 0.001;   %Free Ikba (uM)
    yin(16) = Ms3(i);   %A20 (uM)

[~,Nf]=getScoreCondTot(Par,0,rpkE,myExp3,mnM,yin);
for j=1:2
nfkb(:,i,j)=Nf(j,:);
end
end

for j=1:2
for i=1:mm
    fpeak1=max(nfkb(tt<200,i,j));
    fpeak2=max(nfkb(tt>240,i,j));
    rpkNr(i,j)=fpeak2./fpeak1;
end
end

rpkNDr{1}=rpkNr(:,1);
rpkNDr{2}=rpkNr(:,2);

for j=1:2
Q(:,j,k)= get_distribution(rpkNDr{j}, edges);
d3(j,k) = jsd(P(:,j), Q(:,j,k));

end
end

else   % load the saved scan
    d3=SR.d3;
end
%%
subplot(2,3,3)
lnd=[];
for k=1:4:length(ksc3)
ssM=ksc3(k);
yina20=0.0001;   % CHANGED 2026-10-18: mean initial negative feedback = 1e-4 uM (Table 1); 1e-3 did not reproduce Figure S6B
xx=[1e-7:1e-7:3e-4];   % CHANGED 2026-10-18: plot range of the published panel
    plot(xx,gampdf(xx,yina20/ssM,ssM),"LineWidth",2)
    xlim([0 3e-4])   % CHANGED 2026-10-18: as in the published panel
    ylim([0 5e4])   % CHANGED 2026-10-18: as in the published panel
ylabel("Probability density") 
xlabel("Initial Neg. FB expression (\muM)")   % CHANGED 2026-10-18: label as in the published panel
hold on
lnd=[lnd ["\theta="+num2str(ksc3(k))]];
end
axis square
legend(lnd,"Location","Best")
set(gca,"FontSize",10)

subplot(2,3,6)
plot(ksc3, d3',"LineWidth",2)
xlabel("Scale parameter")
ylabel("Jensen-Shannon distance")
ylim([0 1])
% xlim([10 100000])
% xticks([10,100,1000,10000,100000])
legend({"PIC->TNF","CpG->LPS"},"Location","Best")
axis square
set(gca,"XScale","log","FontSize",10)
title("Initial Neg. FB variability")   % CHANGED 2026-10-18: title as in the published panel

%% Shannon Distance for Total NFkB variability  -> Figure S6A (Gamma distributions) and S6B (JS distance) panel

rng(0)
tt=0:6:474;
nfkb=[];
mm=1000;
mnM=3000;
[~, edges] = histcounts(rmoutliers([rpkEDr{5}' rpkEDr{6}'],"percentiles",[1 99]), 'BinMethod', 'scott');

P=[];
P(:,1)= get_distribution(rpkEDr{5}, edges);
P(:,2)= get_distribution(rpkEDr{6}, edges);

Q=[];
d4=[];
nfkbDis=[];

Ms=2000;
ksc4=10.^[-3:5/20:1];

    yin=zeros(23,1);
    yin(1) = Ms;   %Receptors TNF
    yin(2) = Ms;   %Receptors LPS
    yin(13) = 0;   %Nuclear NFkB
    yin(15) = 0.001;   %Free Ikba (uM)
    yin(16) = 0.0001;   %A20 (uM)

if runSimulation   % simulate 1000 cells for each theta (skipped when loading saved results)
for k=1:length(ksc4)
    k

ssM=ksc4(k);
yinnfkb=1;
Ms3=gamrnd(yinnfkb/ssM,ssM,[mm 1]);
parfor i=1:mm
    Par=[];
    % Par(1)=2.5*10^-5;       %sc IFN translation
    % Par(2)=1*10^-5;          %sc TNF translation
    Par(3)=Ms3(i);               %Total NF-kB
    Par(4)=2;               %Total IKK
    Par(5)=1;               %Total Ikkk

[~,Nf]=getScoreCondTot(Par,0,rpkE,myExp3,mnM,yin);
for j=1:2
nfkb(:,i,j)=Nf(j,:);
end
end

for j=1:2
for i=1:mm
    fpeak1=max(nfkb(tt<200,i,j));
    fpeak2=max(nfkb(tt>240,i,j));
    rpkNr(i,j)=fpeak2./fpeak1;
end
end

rpkNDr{1}=rpkNr(:,1);
rpkNDr{2}=rpkNr(:,2);

for j=1:2
Q(:,j,k)= get_distribution(rpkNDr{j}, edges);
d4(j,k) = jsd(P(:,j), Q(:,j,k));

end
end

else   % load the saved scan
    d4=SR.d4;
end
%%
fig=figure("Position",[20 100 1400 800]);
subplot(2,3,1)
lnd=[];

for k=1:4:length(ksc4)
ssM=ksc4(k);
yinnfkb=1;
xx=[.001:.001:2];
    plot(xx,gampdf(xx,yinnfkb/ssM,ssM),"LineWidth",2)
    xlim([0 2])
    ylim([0 15])
ylabel("Probability density") 
xlabel("Total NF-kB expression (\muM)")
hold on
lnd=[lnd ["\theta="+num2str(ksc4(k))]];
end
axis square
hold off
legend(lnd,"Location","Best")
set(gca,"FontSize",10)

subplot(2,3,4)
plot(ksc4, d4',"LineWidth",2)
xlabel("Scale parameter")
ylabel("Jensen-Shannon distance")
ylim([0 1])
xlim([10^-3 10])
% xticks([10,100,1000,10000,100000])
legend({"PIC->TNF","CpG->LPS"},"Location","Best")
axis square
set(gca,"XScale","log","FontSize",10)
title("Total NF-kB variability")

%% Shannon Distance for Total IKK variability  -> Figure S6A (Gamma distributions) and S6B (JS distance) panel

rng(0)
tt=0:6:474;
nfkb=[];
mm=1000;
mnM=3000;
[~, edges] = histcounts(rmoutliers([rpkEDr{5}' rpkEDr{6}'],"percentiles",[1 99]), 'BinMethod', 'scott');

P=[];
P(:,1)= get_distribution(rpkEDr{5}, edges);
P(:,2)= get_distribution(rpkEDr{6}, edges);

Q=[];
d5=[];
ikkDis=[];

Ms=2000;
ksc5=10.^[-3:5/20:2];

    yin=zeros(23,1);
    yin(1) = Ms;   %Receptors TNF
    yin(2) = Ms;   %Receptors LPS
    yin(13) = 0;   %Nuclear NFkB
    yin(15) = 0.001;   %Free Ikba (uM)
    yin(16) = 0.0001;   %A20 (uM)

if runSimulation   % simulate 1000 cells for each theta (skipped when loading saved results)
for k=1:length(ksc5)
    k
ssM=ksc5(k);
yinIkk=2;
Ms3=gamrnd(yinIkk/ssM,ssM,[mm 1]);
parfor i=1:mm
    Par=[];
    % Par(1)=2.5*10^-5;       %sc IFN translation
    % Par(2)=1*10^-5;          %sc TNF translation
    Par(3)=1;               %Total NF-kB
    Par(4)=Ms3(i);               %Total IKK
    Par(5)=1;               %Total Ikkk

[~,Nf]=getScoreCondTot(Par,0,rpkE,myExp3,mnM,yin);
for j=1:2
nfkb(:,i,j)=Nf(j,:);
end
end

for j=1:2
for i=1:mm
    fpeak1=max(nfkb(tt<200,i,j));
    fpeak2=max(nfkb(tt>240,i,j));
    rpkNr(i,j)=fpeak2./fpeak1;
end
end

rpkNDr{1}=rpkNr(:,1);
rpkNDr{2}=rpkNr(:,2);

for j=1:2
Q(:,j,k)= get_distribution(rpkNDr{j}, edges);
d5(j,k) = jsd(P(:,j), Q(:,j,k));

end
end

else   % load the saved scan
    d5=SR.d5;
end
%%
% fig=figure("Position",[2000 100 1400 800]);
subplot(2,3,2)
lnd=[];

for k=1:4:length(ksc5)
ssM=ksc5(k);
yinIkk=1;
xx=[.001:.001:2];
    plot(xx,gampdf(xx,yinIkk/ssM,ssM),"LineWidth",2)
    xlim([0 2])
    ylim([0 15])
ylabel("Probability density") 
xlabel("Total IKK expression (\muM)")
hold on
lnd=[lnd ["\theta="+num2str(ksc5(k))]];
end
axis square
hold off
legend(lnd,"Location","Best")
set(gca,"FontSize",10)

subplot(2,3,5)
plot(ksc5, d5',"LineWidth",2)
xlabel("Scale parameter")
ylabel("Jensen-Shannon distance")
ylim([0 1])
xlim([10^-3 100])
% xticks([10,100,1000,10000,100000])
legend({"PIC->TNF","CpG->LPS"},"Location","Best")
axis square
set(gca,"XScale","log","FontSize",10)
title("Total IKK variability")

%% Shannon Distance for Total IKKK, IRF and STAT variability  -> Figure S6A (Gamma distributions) and S6B (JS distance) panel

rng(0)
tt=0:6:474;
nfkb=[];
mm=1000;
mnM=3000;
[~, edges] = histcounts(rmoutliers([rpkEDr{5}' rpkEDr{6}'],"percentiles",[1 99]), 'BinMethod', 'scott');

P=[];
P(:,1)= get_distribution(rpkEDr{5}, edges);
P(:,2)= get_distribution(rpkEDr{6}, edges);

Q=[];
d6=[];
ikkDis=[];

Ms=2000;
ksc6=10.^[-3:5/20:2];

    yin=zeros(23,1);
    yin(1) = Ms;   %Receptors TNF
    yin(2) = Ms;   %Receptors LPS
    yin(13) = 0;   %Nuclear NFkB
    yin(15) = 0.001;   %Free Ikba (uM)
    yin(16) = 0.0001;   %A20 (uM)

if runSimulation   % simulate 1000 cells for each theta (skipped when loading saved results)
for k=1:length(ksc6)
    k

ssM=ksc6(k);
yinIkkk=1;
Ms3=gamrnd(yinIkkk/ssM,ssM,[mm 1]);
parfor i=1:mm
    Par=[];
    % Par(1)=2.5*10^-5;       %sc IFN translation
    % Par(2)=1*10^-5;          %sc TNF translation
    Par(3)=1;               %Total NF-kB
    Par(4)=2;               %Total IKK
    Par(5)=Ms3(i);               %Total Ikkk

[~,Nf]=getScoreCondTot(Par,0,rpkE,myExp3,mnM,yin);
for j=1:2
nfkb(:,i,j)=Nf(j,:);
end
end

for j=1:2
for i=1:mm
    fpeak1=max(nfkb(tt<200,i,j));
    fpeak2=max(nfkb(tt>240,i,j));
    rpkNr(i,j)=fpeak2./fpeak1;
end
end

rpkNDr{1}=rpkNr(:,1);
rpkNDr{2}=rpkNr(:,2);

for j=1:2
Q(:,j,k)= get_distribution(rpkNDr{j}, edges);
d6(j,k) = jsd(P(:,j), Q(:,j,k));

end
end

else   % load the saved scan
    d6=SR.d6;
end
%%
% fig=figure("Position",[2000 100 1400 800]);
subplot(2,3,3)
lnd=[];

for k=1:4:length(ksc6)
ssM=ksc6(k);
yinIkkk=1;
xx=[.001:.001:2];
    plot(xx,gampdf(xx,yinIkkk/ssM,ssM),"LineWidth",2)
    xlim([0 2])
    ylim([0 15])
ylabel("Probability density") 
xlabel("Total IKKK/IRF/STAT expression (\muM)")
hold on
lnd=[lnd ["\theta="+num2str(ksc6(k))]];
end
axis square
hold off
legend(lnd,"Location","Best")
set(gca,"FontSize",10)

subplot(2,3,6)
plot(ksc6, d6',"LineWidth",2)
xlabel("Scale parameter")
ylabel("Jensen-Shannon distance")
legend({"PIC->TNF","CpG->LPS"},"Location","Best")
axis square
xlim([0.001 100])
ylim([0 1])
set(gca,"XScale","log","FontSize",10)
title("Total IKKK/IRF/STAT variability")

%% Shannon Distance for TLR4 and TNFR variability  -> Figure S6A (Gamma distributions) and S6B (JS distance) panel

rng(0)
tt=0:6:474;
nfkb=[];
mm=1000;
mnM=3000;
[~, edges] = histcounts(rmoutliers([rpkEDr{5}' rpkEDr{6}'],"percentiles",[1 99]), 'BinMethod', 'scott');

P=[];
P(:,1)= get_distribution(rpkEDr{5}, edges);
P(:,2)= get_distribution(rpkEDr{6}, edges);

Q=[];
d7=[];

ksc7=10.^[0:5/20:5];

Par=[];
% Par(1)=2.5*10^-5;       %sc IFN translation
% Par(2)=1*10^-5;          %sc TNF translation
Par(3)=1;               %Total NF-kB
Par(4)=2;               %Total IKK
Par(5)=1;               %Total Ikkk

if runSimulation   % simulate 1000 cells for each theta (skipped when loading saved results)
for k=1:length(ksc7)
    k
ssM=ksc7(k);

Mms=2000;
Ms2=gamrnd(Mms/ssM,ssM,[mm 1]);

parfor i=1:mm
    yin=zeros(23,1);
    yin(1) = Ms2(i);   %Receptors TNF
    yin(2) = Ms2(i);   %Receptors LPS
    yin(13) = 0;   %Nuclear NFkB
    yin(15) = 0.001;   %Free Ikba (uM)
    yin(16) = 0.0001;   %A20 (uM)

[~,Nf]=getScoreCondTot(Par,0,rpkE,myExp3,mnM,yin);
for j=1:2
nfkb(:,i,j)=Nf(j,:);
end
end

for j=1:2
for i=1:mm
    fpeak1=max(nfkb(tt<200,i,j));
    fpeak2=max(nfkb(tt>240,i,j));
    rpkNr(i,j)=fpeak2./fpeak1;
end
end

rpkNDr{1}=rpkNr(:,1);
rpkNDr{2}=rpkNr(:,2);

for j=1:2
Q(:,j,k)= get_distribution(rpkNDr{j}, edges);
d7(j,k) = jsd(P(:,j), Q(:,j,k));

end
end

else   % load the saved scan
    d7=SR.d7;
end
%%
fig=figure("Position",[20 100 1400 800]);
subplot(2,3,1)
lnd=[];
for k=1:4:length(ksc7)
ssM=ksc7(k);
Mms=2000;
xx=[0:10:8000];
    plot(xx,gampdf(xx,Mms/ssM,ssM),"LineWidth",2)
    xlim([0 6000])
    ylim([0 .006])
ylabel("Probability density") 
xlabel("Initial TLR 4 and TNF expression (\muM)")
hold on
lnd=[lnd ["\theta="+num2str(ksc6(k))]];
end
axis square
hold off
legend(lnd,"Location","Best")
set(gca,"FontSize",10)

subplot(2,3,4)
plot(ksc7, d7',"LineWidth",2)
xlabel("Scale parameter")
ylabel("Jensen-Shannon distance")
legend({"PIC->TNF","CpG->LPS"},"Location","Best")
ylim([0 1])
xlim([1 100000])
xticks([10,1000,100000])
axis square
set(gca,"XScale","log","FontSize",10)
title("Initial TLR 4 and TNF variability")





%%

%% Save the scans (only when they were simulated)
if runSimulation
    if ~isfolder(fileparts(resFile)), mkdir(fileparts(resFile)); end
    save(resFile, 'd1','d2','d3','d4','d5','d6','d7','ksc1','ksc2','ksc3','ksc4','ksc5','ksc6','ksc7','Q1', '-v7');
end

function P = get_distribution(x, edges)
    counts = histcounts(x, edges, 'Normalization', 'probability');
    P = counts + 1e-12; % avoid zeros
    P = P / sum(P);
end

function d = jsd(P, Q)
    % Ensure row vectors
    P = P(:)' + 1e-12;
    Q = Q(:)' + 1e-12;

    % Normalize
    P = P / sum(P);
    Q = Q / sum(Q);

    M = 0.5 * (P + Q);

    % KL divergence
    KL_PM = sum(P .* log2(P ./ M));
    KL_QM = sum(Q .* log2(Q ./ M));

    % Jensen-Shannon distance
    d = sqrt(0.5 * (KL_PM + KL_QM));
end