%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4 related to
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
%
% Wang AG, Mercado Vasquez JG, Son M, et al
%
% Corresponding author: tays@uchicago.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cytokine blocking experiments (Figure 4D-I). Figure 4A-C and 4J are RNA-seq analyses (R/Fig4_RNAseq.R).
% Single-cell data: table data4 (one row per cell and phase), built by functions/genBlockFeat.m.
% AUCs are normalized to the mean of the control in each replicate. Black line: mean.
% Statistics: two-sided Wilcoxon rank-sum test, Bonferroni-corrected within each panel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
addpath functions simulation_results classifier_results data/flow   % run from the main code folder: all paths are relative to the current folder
if ~isfile("single_cell_features/genBlockFeat.mat"), genBlockFeat; end
load single_cell_features/genBlockFeat.mat data4
groupcounts(data4,["Experiment","Category","Phase"])   %number of cells per condition and phase

red=[0.62 0.13 0.13; 0.98 0.5 0.45]; blue=[0.20 0.40 0.73; 0.65 0.82 1]; gry=[0.5 0.5 0.5];

countlabels(data4.Category)
%% Fig 4D
%Tracemaps of the phase I response with or without sTNFR1 (100 random cells, earliest responders at the top)
ex=["sTNFR1 (polyI:C)","sTNFR1 (polyI:C)","sTNFR1 (CpG)","sTNFR1 (CpG)"];
cats=["PolyIC -> TNF","PolyIC + sTNFR1 -> TNF","CpG -> LPS","CpG + sTNFR1 -> LPS"];
figure("Position",[50 50 900 400])
for i=1:4
    if mod(i,2)==1, rng(12), end        %same random cells as the original analysis
    zz=data4.Experiment==ex(i) & data4.Category==cats(i) & data4.Phase==1;
    nfkb=data4.NfkB(zz,:);
    nfkb=nfkb(datasample(1:size(nfkb,1),100,'Replace',false),:);
    auc=trapz(nfkb(:,1:10),2);
    [~,idx]=sort(auc,"descend");

    subplot(1,4,i)
    imagesc(0:6:234,1:100,nfkb(idx,:),[0 5])
    xlabel("Time (min)")
    yticks([])
    title(cats(i))
end
colorbar

%% Fig 4E-G
%Violin plots of the normalized AUC, control vs block
% {experiment, phase, feature, y label, colors, x labels}
pan={ {"sTNFR1 (CpG)",1,"Late AUC norm","CpG Late AUC",red,["Con","Blk"]}, {"sTNFR1 (polyI:C)",1,"Late AUC norm","PIC Late AUC",blue,["Con","Blk"]}, ... %Fig 4E
      {"sTNFR1 (CpG)",2,"AUC norm","LPS AUC",red,["Con","Blk"]},           {"sTNFR1 (polyI:C)",2,"AUC norm","TNF AUC",blue,["Con","Blk"]}, ...      %Fig 4F
      {"sIL10R",2,"AUC norm","LPS AUC",red,["CpG","CpG + IL-10 block"]} };                                                                        %Fig 4G
figure("Position",[50 50 1100 650])
fprintf('\nPanel   n (control, block)   p (two-sided rank-sum)\n')

for i=1:5
    p=pan{i};
    zz=data4.Experiment==p{1} & data4.Phase==p{2};
    g=unique(data4.Category(zz),'stable');                 %control first, block second
    v={data4.(p{3})(zz & data4.Category==g(1)), data4.(p{3})(zz & data4.Category==g(2))};
    
    subplot(2,4,i)
    violin3({max(v{1},0),max(v{2},0)},'facecolor',p{5},'mc','k','medc',[],'bw',0.3);
    pv=ranksum(v{1},v{2});
    addStar(1,2,4.2,pv)
    xticks([1 2])
    xticklabels(p{6})
    ylabel(p{4})
    ylim([0 4.6])
    set(gca,'FontSize',8)
    fprintf('%-24s %4d, %4d   %.2g\n',p{1}+" "+p{4},numel(v{1}),numel(v{2}),pv)
end

%% Fig 4H-I
%Naive, polyI:C and polyI:C + sIFNAR1 followed by TNF (H: AUC, I: late AUC), Bonferroni x3
cats=["Control -> TNF","PolyIC 1000 -> TNF","PolyIC 1000 + sIFNAR1 -> TNF"];
feat=["AUC norm","Late AUC norm"]; ylb=["TNF-\alpha AUC","TNF-\alpha late AUC"];
pairs=[1 2; 2 3; 1 3];
shown={[1 2 3],[2 3]};                                   %comparisons drawn in Fig 4H and 4I
for j=1:2
    v=cell(1,3);
    for k=1:3
        zz=data4.Experiment=="sIFNAR1" & data4.Category==cats(k) & data4.Phase==2;
        v{k}=rmmissing(data4.(feat(j))(zz));
    end

    subplot(2,4,5+j)
    violin3(cellfun(@(x) max(x,0),v,'UniformOutput',false),'facecolor',[gry;blue],'mc','k','medc',[],'bw',0.2);
    for q=1:3
        pv=min(1,3*ranksum(v{pairs(q,1)},v{pairs(q,2)}));
        if ismember(q,shown{j}), addStar(pairs(q,1),pairs(q,2),3.9+0.35*(q-1),pv), end
        fprintf('%-24s %4d, %4d   %.2g (Bonferroni x3)\n',cats(pairs(q,1))+" vs "+cats(pairs(q,2)),numel(v{pairs(q,1)}),numel(v{pairs(q,2)}),pv)
    end
    xticks(1:3)
    xticklabels(["Naive","PolyI:C","PolyI:C + IFN block"])
    xtickangle(30)
    ylabel(ylb(j))
    ylim([0 4.9])
    set(gca,'FontSize',8)
end

%%
function addStar(a,b,y,p)
%Comparison bar with * (p < 0.01) or n.s.
line([a b],[y y],'Color','k')
if p<0.01, s="*"; else, s="n.s."; end
text((a+b)/2,y+0.12,s,'HorizontalAlignment','center','FontSize',8)
end
