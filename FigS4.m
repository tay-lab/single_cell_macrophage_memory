%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S4 related to
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
%
% Wang AG, Mercado-Vásquez G, Son M, et al
%
% Corresponding author: tays@uchicago.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single-cell results of cytokine blocking with soluble receptors (Figure S4C-F): tracemaps of random cells,
% sorted by maximum (strongest responders at the top). Figure S4A-B and S4G are RNA-seq analyses (R/).
% Single-cell data: table data4 (one row per cell and phase), built by functions/genBlockFeat.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
addpath functions simulation_results classifier_results data/flow   % run from the main code folder: all paths are relative to the current folder
if ~isfile("single_cell_features/genBlockFeat.mat"), genBlockFeat; end
load single_cell_features/genBlockFeat.mat data4
tt=0:6:234;                                              %time of the 40 frames (min)

%% Fig S4C-D
%Second-stimulus response with or without sTNFR1 (100 random cells per tracemap)
% {experiment, categories (control, block), random seed, panel}
pan={ {"sTNFR1 (CpG)",["CpG -> LPS","CpG + sTNFR1 -> LPS"],12,"S4C: LPS after CpG"}, ...
      {"sTNFR1 (polyI:C)",["PolyIC -> TNF","PolyIC + sTNFR1 -> TNF"],12,"S4D: TNF after polyI:C"} };
blk=["Control","+ sTNFR1"];
figure("Name","Figure S4C-D","Position",[50 50 700 600])
for i=1:2
    p=pan{i};
    rng(p{3})                                            %fixed seed: the same random cells are plotted at every run
    for b=1:2
        zz=data4.Experiment==p{1} & data4.Category==p{2}(b) & data4.Phase==2;
        subplot(2,2,2*(i-1)+b)
        plotMap(tt,data4.NfkB(zz,:),100,[0 3])
        title(p{4}+" "+blk(b))
    end
end
colorbar

%% Fig S4E
%PolyI:C (300 ng/ml) followed by TNF, with or without sIFNAR1 (600 random cells per tracemap).
%Phase I traces are normalized to the mean maximum of each tracemap, phase II traces to the mean maximum without sIFNAR1.
cats=["PolyIC 300 -> TNF","PolyIC 300 + sIFNAR1 -> TNF"];
stim=["PolyI:C","TNF\alpha"]; blk=["Control","+ sIFNAR1"];
rng(16)
figure("Name","Figure S4E","Position",[50 50 700 600])
for b=1:2
    for ph=1:2
        zz=data4.Experiment=="sIFNAR1" & data4.Category==cats(b) & data4.Phase==ph;
        nfkb=datasample(data4.NfkB(zz,:),600,1,'Replace',false);
        if ph==1, norm=mean(max(nfkb,[],2)); elseif b==1, norm2=mean(max(nfkb,[],2)); norm=norm2; else, norm=norm2; end
        subplot(2,2,2*(b-1)+ph)
        plotMap(tt,nfkb/norm,600,[0 1.5])
        title(stim(ph)+" "+blk(b))
    end
end
colorbar

%% Fig S4F
%CpG followed by LPS, with or without sIL10R (100 random cells per tracemap)
cats=["CpG -> LPS","CpG + sIL10R -> LPS"];
stim=["CpG","LPS"]; blk=["Control","+ sIL10R"];
rng(9)
figure("Name","Figure S4F","Position",[50 50 700 600])
for b=1:2
    for ph=1:2
        zz=data4.Experiment=="sIL10R" & data4.Category==cats(b) & data4.Phase==ph;
        subplot(2,2,2*(b-1)+ph)
        plotMap(tt,data4.NfkB(zz,:),100,[0 4])
        title(stim(ph)+" "+blk(b))
    end
end
colorbar

%%
function plotMap(tt,nfkb,n,clim)
%Tracemap of n random cells (all if fewer), sorted by maximum
if size(nfkb,1)>n, nfkb=nfkb(datasample(1:size(nfkb,1),n,'Replace',false),:); end
[~,idx]=sort(max(nfkb,[],2),"descend");
imagesc(tt,1:size(nfkb,1),nfkb(idx,:),clim)
xlabel("Time (min)")
yticks([])
set(gca,'FontSize',8)
end
