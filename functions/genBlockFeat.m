function data4 = genBlockFeat()
% genBlockFeat  Single-cell table of the cytokine-blocking experiments (Figure 4D-I, Figure S4C-F).
%
%   data4 = genBlockFeat() reads the four files in single_cell_features/ and returns one table with one
%   row per cell and phase, organized like data0 (genExtFeat.mat):
%     Experiment  "sTNFR1 (polyI:C)", "sTNFR1 (CpG)", "sIL10R" or "sIFNAR1"
%     Category    condition, e.g. "PolyIC -> TNF", "PolyIC + sTNFR1 -> TNF"
%     Replicate   independent experiment (1-3)
%     Phase       1 = first stimulus (polyI:C / CpG), 2 = second stimulus (TNF / LPS)
%     NfkB        single-cell nuc/cyto NF-kB trace (40 frames, 6 min/frame)
%     Max, AUC, Late AUC (frames 21-40 = 120-240 min), Early AUC (phase 2 only)
%     AUC norm, Late AUC norm: AUCs divided by the mean of the control condition of the same experiment,
%       replicate and phase (sIFNAR1 experiment: mean of "Control -> TNF", replicates pooled)
%   The table is saved to single_cell_features/genBlockFeat.mat.
%
%   Source columns of the peakfeatures files (one row per chamber): 1/6 = phase 1/2 traces;
%   2-5 = phase 1 max, peak width, AUC, late AUC; 7-10 = phase 2 max, AUC, early AUC, late AUC;
%   11/12 = first/second stimulus code; 14 = replicate. Second-stimulus codes: 6/4 = TNF, 7/5 = LPS.

d = 'single_cell_features';
% {file, experiment, second-stimulus codes, second stimulus, {code, category} per replicate (control first)}
E = {
 'peakfeatures_PICfirst_TNFblock.mat',  "sTNFR1 (polyI:C)", [6 4], "TNF", ...
     {{2,"PolyIC";3,"PolyIC + sTNFR1"}, {2,"PolyIC";3,"PolyIC + sTNFR1"}, {0,"PolyIC";1,"PolyIC + sTNFR1"}}
 'peakfeatures_CPGfirst_TNFblock.mat',  "sTNFR1 (CpG)",     [7 5], "LPS", ...
     {{2,"CpG";3,"CpG + sTNFR1"}, {3,"CpG";4,"CpG + sTNFR1"}}
 'peakfeatures_CPGfirst_IL10block.mat', "sIL10R",           [7 5], "LPS", ...
     {{3,"CpG";21,"CpG + sIL10R"}, {2,"CpG";7,"CpG + sIL10R"}}
};
data4 = table;
for e = 1:size(E,1)
    S = load(fullfile(d,E{e,1})); fn = fieldnames(S); P = S.(fn{1});
    for r = 1:numel(E{e,5})
        grp = E{e,5}{r};
        for g = 1:size(grp,1)
            rows = cell2mat(P(:,11))==grp{g,1} & ismember(cell2mat(P(:,12)),E{e,3}) & cell2mat(P(:,14))==r;
            data4 = [data4; cellRows(P(rows,:), E{e,2}, grp{g,2}+" -> "+E{e,4}, r)]; %#ok<AGROW>
        end
    end
end
% sIFNAR1 experiment (2024): all eight first-stimulus groups, second stimulus TNF (6) or LPS (7)
S = load(fullfile(d,'peakfeatures_IFNblock.mat')); P = S.peakfeatures_IFNblock([1:62,64],:);   % row 63 is excluded
codes = [16 17 0 3 1 4 2 5];
stim2 = ["TNF","LPS"];
names = ["Control","Control + sIFNAR1","CpG","CpG + sIFNAR1","PolyIC 300","PolyIC 300 + sIFNAR1","PolyIC 1000","PolyIC 1000 + sIFNAR1"];
for s2 = [6 7]
    for g = 1:numel(codes)
        for r = 1:2
            rows = cell2mat(P(:,11))==codes(g) & cell2mat(P(:,12))==s2 & cell2mat(P(:,14))==r;
            data4 = [data4; cellRows(P(rows,:), "sIFNAR1", names(g)+" -> "+stim2(s2-5), r)]; %#ok<AGROW>
        end
    end
end
data4.Experiment = categorical(data4.Experiment); data4.Category = categorical(data4.Category);

% Normalized AUCs: control of the same experiment, replicate and phase (sIFNAR1: "Control -> TNF", pooled)
data4.("AUC norm") = nan(height(data4),1); data4.("Late AUC norm") = nan(height(data4),1);
ctrl = ["PolyIC -> TNF","CpG -> LPS","CpG -> LPS"]; ex = ["sTNFR1 (polyI:C)","sTNFR1 (CpG)","sIL10R"];
for e = 1:3
    for r = unique(data4.Replicate(data4.Experiment==ex(e)))'
        for ph = 1:2
            zz = data4.Experiment==ex(e) & data4.Replicate==r & data4.Phase==ph;
            z0 = zz & data4.Category==ctrl(e);
            data4.("AUC norm")(zz) = data4.AUC(zz)/mean(data4.AUC(z0),'omitnan');
            data4.("Late AUC norm")(zz) = data4.("Late AUC")(zz)/mean(data4.("Late AUC")(z0),'omitnan');
        end
    end
end
for ph = 1:2
    zz = data4.Experiment=="sIFNAR1" & data4.Phase==ph & endsWith(string(data4.Category),"-> TNF");
    z0 = zz & data4.Category=="Control -> TNF";
    data4.("AUC norm")(zz) = data4.AUC(zz)/mean(data4.AUC(z0),'omitnan');
    data4.("Late AUC norm")(zz) = data4.("Late AUC")(zz)/mean(data4.("Late AUC")(z0),'omitnan');
end
save(fullfile(d,'genBlockFeat.mat'),'data4')
end

function T = cellRows(C, experiment, category, rep)
% One row per cell and phase from the chamber rows C of a peakfeatures file
T = table;
for ph = 1:2
    if ph==1, tr = cell2mat(C(:,1)); f = cellfun(@(i) cell2mat(C(:,i)), {2,4,5},'UniformOutput',false); ea = nan(size(f{1}));
    else,     tr = cell2mat(C(:,6)); f = cellfun(@(i) cell2mat(C(:,i)), {7,8,10},'UniformOutput',false); ea = cell2mat(C(:,9));
    end
    n = size(tr,1);
    t = table(repmat(experiment,n,1), repmat(category,n,1), repmat(rep,n,1), repmat(ph,n,1), tr(:,1:40), ...
              f{1}, f{2}, ea, f{3}, 'VariableNames', ["Experiment","Category","Replicate","Phase","NfkB","Max","AUC","Early AUC","Late AUC"]);
    T = [T; t]; %#ok<AGROW>
end
end
