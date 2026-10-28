%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S2 related to
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
% Wang AG, Mercado-Vásquez G, Son M, et al. Corresponding author: tays@uchicago.edu
%
% Single-cell differences in the remaining 8 NF-kB features (Figure S2A, Phase I; Figure S2B, Phase II).
% Features are ranked by mutual information as in Figure 1D; the top 4 are shown in Figure 1E-H and features
% 5-12 here. Data: functions/genExtFeat.mat (n = 164 Control -> LPS, 973 CpG -> LPS, 322 Control -> TNF,
% 1795 polyI:C -> TNF cells). Black line = median. Runtime: seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions simulation_results classifier_results data/flow   % run from the main code folder: all paths are relative to the current folder
load genExtFeat.mat MIfeatE MIfeatL datvioE datvioL
fn=[" NFkB MI_{max}"," Activation Time"," Deactivation Time"," Speed"," Duration"," Oscillations",...
    " Amplitude"," Time to Peak"," Prominence"," AUC"," Early AUC"," Late AUC"];   % feature labels (as in Fig1.m)
gry=[0.6 0.6 0.6]; ltred=[0.98 0.49 0.49]; dkred=[0.8 0.15 0.15]; ltblue=[0.45 0.8 0.96]; dkblue=[0.26 0.56 0.69];

% Panel order as in the published figure (row by row within each pair of columns). These are the features
% ranked 5-12 by mutual information (checked below); their order within the figure was arranged by hand.
% Feature numbers: 1 MI_max, 2 Activation, 3 Deactivation, 4 Speed, 5 Duration, 6 Oscillations,
% 7 Amplitude, 8 Time to Peak, 9 Prominence, 10 AUC, 11 Early AUC, 12 Late AUC.
% {data, MI, x labels, colors, panel order (row 1: CpG/LPS, row 2: polyI:C/TNF), figure name}
phases={ {datvioE,MIfeatE,[{"Ctrl","CpG"};{"Ctrl","PolyIC"}],{[gry;ltred],[gry;ltblue]}, ...
          [10 12 4 11 3 8 2 6; 5 10 4 11 3 8 2 6],"Figure S2A: Phase I"}, ...
         {datvioL,MIfeatL,[{"LPS_{naive}","LPS_{cond}"};{"TNF_{naive}","TNF_{cond}"}],{[ltred;dkred],[ltblue;dkblue]}, ...
          [1 5 9 4 11 8 2 6; 9 7 10 4 11 8 2 6],"Figure S2B: Phase II"} };
for ph=1:2
    [dv,MI,xlb,col,order]=phases{ph}{1:5};
    figure("Name",phases{ph}{6},"Position",[50 50 900 800])
    for jj=1:2                                   % jj=1: CpG/LPS (columns 1-2), jj=2: polyI:C/TNF (columns 3-4)
        [~,idxp]=sort(MI(:,jj),"descend");       % rank features by mutual information
        assert(isequal(sort(idxp(5:12))',sort(order(jj,:))),"panels are not the features ranked 5-12")
        panelPos=[1 2 5 6 9 10 13 14]+2*(jj-1);  % subplot positions of the 8 panels (4 rows x 2 columns)
        for k=1:8
            subplot(4,4,panelPos(k))
            ii=order(jj,k);
            bw=median(dv{ii,jj,2},"omitnan")/3;  % bandwidth as in Fig1.m
            evalc("violin3({dv{ii,jj,1},dv{ii,jj,2}},'facecolor',col{jj},'mc',[],'medc','black','bw',bw,'facealpha',1)");   % evalc hides violin3's console messages
            xticks([1,2]); xticklabels(xlb(jj,:)); ylabel(fn{ii})
            set(gca,'YLim',[0 Inf],'FontSize',8)
        end
    end
end
