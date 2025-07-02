%% TNFblock PIC features
clc
clear
load("./single_cell_features/peakfeatures_PICfirst_TNFblock.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

%6 and 4 for TNF, 7 and 5 for LPS
conditions = peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
%For PIC
inputs = [17, 19, 17, 19, 17, 19; 0, 1, 0, 1, 0, 1; 2, 3, 2, 3, 0, 1; 4, 5, 4, 5, 0, 1 ];


maxnorms = [mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 3, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 3, 5) ), feature) ))];

figure(3) %plotting TNFblock effects
alltraces = {};
clf
%For PIC
cond = conditions(ismember( conditions(:, 1), inputs(3, :)), :);
for bb = 1:2
    temp = [];
    for cc = 1:3
        temp1 = cell2mat(  peakfeatures_PICfirst_TNFblock( ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        %Norm to no block PIC
        temp1 =  temp1./mean(cell2mat( peakfeatures_PICfirst_TNFblock(ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, 1+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ));
        temp = [temp;temp1];
    end
    length(temp)
    Violin( max( temp, 0), bb, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
        'BoxWidth', 0.1, 'ShowMean', true);
    alltraces{2, bb} = max( temp, 0);
end
set(gca, 'XLim', [0.25, 2.75], 'YLim', [-0.5, 4])
title( inputs(2, 1) );
xticks([1 2])
xticklabels({'Con','Blk'})
ylabel('S2 LPS AUC')

%% TNFblock PIC tracemap
rng(12)
load("./single_cell_features/peakfeatures_PICfirst_TNFblock.mat")
feature = 1; %feature 1 is PIC stimulus traces, feature 6 is TNF stimulus traces
conditions = peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_PICfirst_TNFblock( ([peakfeatures_PICfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
inputs = [17, 19, 17, 19, 17, 19; 0, 1, 0, 1, 0, 1; 2, 3, 2, 3, 0, 1; 4, 5, 4, 5, 0, 1, ];

figure(3) %plotting TNFblock effects
ylim1 = [0 3];

clf
counter = 1;

cond = conditions(ismember( conditions(:, 1), inputs(3, :)), :);
for bb = 1:2
    temp = [];
    tempt = [];
    for cc = 1:3
        temp1 = cell2mat(  peakfeatures_PICfirst_TNFblock( ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        tempt1 = cell2mat(  peakfeatures_PICfirst_TNFblock( ismember( cell2mat( peakfeatures_PICfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(3, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), 7) ); %3 by peak time for first interval, 7 by max
        temp = [temp;temp1];
        tempt = [tempt;tempt1];
    end
    [nr, ~] = size(temp);
    rows = datasample(1:nr, 100, 'replace', false);
    traces = temp(rows, 1:40);
    tracest = tempt(rows);

    [tracest,sortIdx]  = sort(tracest,'descend'); %ascend for time in PIC, descend for peak height in TNF
    subplot(1,2, counter)

    imagesc(traces(sortIdx, :), ylim1)
    length(temp)
    counter = counter + 1;
    axis off
    if bb == 1
        title("Control")
    else
        title("TNF block")
    end

end


%% TNFblock CpG features
clc
clear
load("./single_cell_features/peakfeatures_CPGfirst_TNFblock.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

%6 and 4 for TNF, 7 and 5 for LPS
conditions = peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);

%For CpG
inputs = [17, 19, 17, 19; 2, 3, 3, 4];
maxnorms = [mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 3, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 3, 5) ), feature) ))];

figure(3) %plotting TNFblock effects
clf
hold on
cond = conditions(ismember( conditions(:, 1), inputs(2, :)), :);
for bb = 1:2
    temp = [];
    for cc = 1:2
        temp1 = cell2mat(  peakfeatures_CPGfirst_TNFblock( ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        %Norm to no block PIC
        temp1 =  temp1./mean(cell2mat( peakfeatures_CPGfirst_TNFblock(ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, 1+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ));
        temp = [temp;temp1];
    end
    length(temp)
    if bb == 1
        tempcon = temp;
    else
        tempblk = temp;
    end
    Violin( max( temp, 0), bb, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
        'BoxWidth', 0.1, 'ShowMean', true);
    alltraces{2, bb} = max( temp, 0);
end
set(gca, 'XLim', [0.25, 2.75], 'YLim', [-0.5, 4])
title( inputs(2, 1) );
xticks([1 2])
xticklabels({'Con','Blk'})
ylabel('S2 LPS AUC')

%% TNFblock CpG heatmap
clc
clear
rng(12)
load("./single_cell_features/peakfeatures_CPGfirst_TNFblock.mat")

feature = 6; %1 for CpG, 6 for LPS
conditions = peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 7 )', 11:15);
conditions = [conditions; peakfeatures_CPGfirst_TNFblock( ([peakfeatures_CPGfirst_TNFblock{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
%For CpG
inputs = [17, 19, 17, 19; 2, 3, 3, 4];

figure(3) %plotting TNFblock effects
ylim1 = [0 3];
clf
counter = 1;
cond = conditions(ismember( conditions(:, 1), inputs(2, :)), :);
for bb = 1:2
    temp = [];
    tempt = [];
    for cc = 1:2
        temp1 = cell2mat(  peakfeatures_CPGfirst_TNFblock( ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        tempt1 = cell2mat(  peakfeatures_CPGfirst_TNFblock( ismember( cell2mat( peakfeatures_CPGfirst_TNFblock(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), 7) ); %3 for time for CpG, 7 for peak height for LPS
        temp = [temp;temp1];
        tempt = [tempt;tempt1];
    end
    [nr, ~] = size(temp);
    rows = datasample(1:nr, 100, 'replace', false);
    traces = temp(rows, 1:40);
    tracest = tempt(rows);

    [tracest,sortIdx]  = sort(tracest,'descend'); %ascend for cpg, descend for LPS
    subplot(1,2, counter)
    imagesc(traces(sortIdx, :), ylim1)
    length(temp)
    counter = counter + 1;
    axis off
    if bb == 1
        title("Control")
    else
        title("TNF block")
    end

end

%% IL10block CpG features
clc
clear

load("./single_cell_features/peakfeatures_CPGfirst_IL10block.mat")

%Feature 2-5 = Stim A max, width, AUC, and late AUC
%Feature 7-10 = Stim B max, AUC, early AUC, and late AUC
feature = 8;

%6 and 4 for TNF, 7 and 5 for LPS
conditions = peakfeatures_CPGfirst_IL10block( ([peakfeatures_CPGfirst_IL10block{:, 12}] == 6 )', 11:15);
conditions = [conditions; peakfeatures_CPGfirst_IL10block( ([peakfeatures_CPGfirst_IL10block{:, 12}] == 4 )', 11:15)];
conditions = cell2mat(conditions);
%For CpG
inputs = [17, 19, 17, 19; 3, 21, 2, 7];

maxnorms = [mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 1, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 2, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 17 & conditions(:, 4) == 3, 5) ), feature) )), ...
    mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), conditions(conditions(:, 1) == 19 & conditions(:, 4) == 3, 5) ), feature) ))];

figure(3) %plotting IL10block effects
clf
%For CpG
    hold on
    cond = conditions(ismember( conditions(:, 1), inputs(2, :)), :);
        for bb = 1:2
            temp = [];
            for cc = 1:2
                temp1 = cell2mat(  peakfeatures_CPGfirst_IL10block( ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;

                %Norm to no block PIC
                temp1 =  temp1./mean(cell2mat( peakfeatures_CPGfirst_IL10block(ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), cond(cond(:, 1) == inputs(2, 1+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ));
                temp = [temp;temp1];
            end
            length(temp)
            if bb == 1
                tempcon = temp;
            else
                tempblk = temp;
            end
            Violin( max( temp, 0), bb, 'Bandwidth', 0.3,  'MedianColor', [1,1,1], 'ViolinColor', [0 0 0], 'ViolinAlpha', .3, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
                'BoxWidth', 0.1, 'ShowMean', true);
            alltraces{2, bb} = max( temp, 0);
        end
        set(gca, 'XLim', [0.25, 2.75], 'YLim', [-0.5, 4])
        xticks([1 2])
        xticklabels({'Con','Blk'})
        ylabel('S2 TNF AUC')
%% IL10block CpG heatmap
clc
clear
rng(9)
load("./single_cell_features/peakfeatures_CPGfirst_IL10block.mat")

feature = 6; %1 for CpG, 6 for LPS
%6 and 4 for TNF, 7 and 5 for LPS
conditions = peakfeatures_CPGfirst_IL10block( ([peakfeatures_CPGfirst_IL10block{:, 12}] == 7 )', 11:15);
conditions = [conditions; peakfeatures_CPGfirst_IL10block( ([peakfeatures_CPGfirst_IL10block{:, 12}] == 5 )', 11:15)];
conditions = cell2mat(conditions);
%For CpG
inputs = [17, 19, 17, 19; 3, 21, 2, 7];

figure(3) %plotting TNFblock effects
ylim1 = [0 4];
clf
counter = 1;
cond = conditions(ismember( conditions(:, 1), inputs(2, :)), :);
for bb = 1:2
    temp = [];
    tempt = [];
    for cc = 1:2
        temp1 = cell2mat(  peakfeatures_CPGfirst_IL10block( ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), feature) ) ;
        tempt1 = cell2mat(  peakfeatures_CPGfirst_IL10block( ismember( cell2mat( peakfeatures_CPGfirst_IL10block(:, 15) ), cond(cond(:, 1) == inputs(2, bb+(cc-1)*2) & cond(:, 4) == cc, 5) ), 7) ); %2 for peak height for CpG, 7 for peak height for LPS
        temp = [temp;temp1];
        tempt = [tempt;tempt1];
    end
    [nr, ~] = size(temp);
    rows = datasample(1:nr, 100, 'replace', false);
    traces = temp(rows, 1:40);
    tracest = tempt(rows);

    [tracest,sortIdx]  = sort(tracest,'descend'); %ascend for cpg, descend for LPS
    subplot(1,2, counter)
    imagesc(traces(sortIdx, :), ylim1)
    length(temp)
    counter = counter + 1;
    axis off
    if bb == 1
        title("Control")
    else
        title("IL10 block")
    end

end
%% CPG and PIC IFN block
clc
clear
colormap parula
load("./single_cell_features/peakfeatures_IFNblock.mat")
rng(16)

peakfeatures_IFNblock = peakfeatures_IFNblock([1:62, 64], :);

lig1 = [1, 4];  %0/3 are CpG 100 nm, 1/4 are PIC 300 ng, 2/5 are PIC 1000 ng
lig2 = 6; %6 is 10ng TNF, 7 is 1 ng LPS
conditions = peakfeatures_IFNblock( :, 15);
conditions = cell2mat(conditions);


figure(1)
clf
norm2 = [];
for aa = 1:2
    conditions = 1:63;
    %conditions = cell2mat(conditions);
    conditions = conditions(cell2mat( peakfeatures_IFNblock(conditions, 11) ) == lig1(aa) & cell2mat( peakfeatures_IFNblock(conditions, 12)) == lig2);
    temp1 = [];
    temp2 = [];

    for bb = 1:length(conditions) %1:length(conditions)/2 %1:length(conditions) %1:length(conditions)/2 %(length(conditions)/2+1):length(conditions) 
        temp1 = [temp1; peakfeatures_IFNblock{conditions(bb), 1} ];
        temp2 = [temp2; peakfeatures_IFNblock{conditions(bb), 6} ];
                
    end
    temp1 = datasample(temp1, 600, 1, 'Replace', false);
    temp2 = datasample(temp2, 600, 1, 'Replace', false);
    %temp1 = datasample(temp1,  200,'Replace', false);
    maxes1 = max(temp1(:, 1:40), [], 2);
    norm1 = mean(maxes1);
    [~, idx1] = sort(maxes1, 'descend');    
    
    %temp2 = datasample(temp2,  200,'Replace', false);
    maxes2 = max(temp2(:, 1:40), [], 2);
    if aa == 1
        norm2 = mean(maxes2);
    end

    [~, idx2] = sort(maxes2, 'descend');

    subplot(2,2,aa*2-1)
    imagesc(temp1(idx1, :)./norm1, [0 1.5])

    subplot(2,2,aa*2)
    imagesc(temp2(idx2, :)./norm2, [0 1.5])

end

figure(2)
clf
for aa = 1:2
    subplot(1,2, aa)
    hold on
    for bb = 1:2
        conditions = 1:63;
        conditions = conditions(cell2mat( peakfeatures_IFNblock(conditions, 11) ) == lig1(bb) & cell2mat( peakfeatures_IFNblock(conditions, 12)) == lig2);
        temp1 = [];
        temp2 = [];

        for cc = 1:length(conditions) 
            temp1 = [temp1; peakfeatures_IFNblock{conditions(cc), 1} ];
            temp2 = [temp2; peakfeatures_IFNblock{conditions(cc), 6} ];
        end
        if aa == 1
            plot(smoothdata(median(temp1, 1),"movmean",4), 'LineWidth', 2)
        else
            plot(smoothdata(median(temp2, 1),"movmean",4), 'LineWidth', 2)
        end
                
    end
    %set(gca, 'XLim', [0, 40], "YLim", [-0.1, 1.6]);
    hold off

    
end

%% Second ligand features

lig1 = [16, 17, 2, 5]; %0/3 are CpG 100 nm, 1/4 are PIC 300 ng, 2/5 are PIC 1000 ng
lig2 = 6; %6 is 10ng TNF, 7 is 1 ng LPS
feature1 = 4;
feature2 = 10;

figure(2)
clf
for cc = 1:2
    subplot(1, 2, cc)
    hold on
    if cc == 1

        for aa = 1:4
            pos = [1, 2, 3, 4];
            conditions = 1:63;
            %conditions = cell2mat(conditions);
            conditions = conditions(cell2mat( peakfeatures_IFNblock(conditions, 11) ) == lig1(aa) & cell2mat( peakfeatures_IFNblock(conditions, 12)) == lig2);
            temp1 = [];
            temp2 = [];



            for bb = 1:length(conditions)
                temp1 = [temp1; peakfeatures_IFNblock{cell2mat( peakfeatures_IFNblock(:, 15) ) == conditions(bb), feature1} ];

            end
            temp1 = max(temp1, 0);
            if aa == 1
                set1 = temp1;
            elseif aa ==2
                set2 = temp1;
            elseif aa ==3
                set3 = temp1;
            else
                set4 = temp1;
            end


            Violin( temp1, pos(aa), 'Bandwidth', 4,  'MedianColor', [1,1,1], 'ViolinAlpha', .5, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
                'BoxWidth', 0.1, 'ShowMean', true);
        end
        set(gca, 'XLim', [0.5, 4.5], "YLim", [-1, 150]);
        hold off
    else
        for aa = 1:4
            pos = [1, 2, 3, 4];
            conditions = 1:63;
            %conditions = cell2mat(conditions);
            conditions = conditions(cell2mat( peakfeatures_IFNblock(conditions, 11) ) == lig1(aa) & cell2mat( peakfeatures_IFNblock(conditions, 12)) == lig2);
            temp1 = [];
            temp2 = [];

            for bb = 1:length(conditions) % length(conditions) %/2 gets only first replicate
                foo = peakfeatures_IFNblock{conditions(bb), feature2};
                temp2 = [temp2; foo(~isnan(foo)) ];

            end

            if aa == 1
                normauc2 = mean(temp2);
            end
            temp2 = temp2./normauc2;
%             if aa == 1
%                 normauc2 = mean(temp2);
%             elseif aa == 2
%                 ifnarauc2 = mean(temp2);
%             end
%             if aa == 1 || aa == 3
%                 temp2 = temp2./normauc2;
%             else
%                 temp2 = temp2./ifnarauc2;
%             end
            if aa == 1
                set5 = temp2;
            elseif aa ==2
                set6 = temp2;
            elseif aa ==3
                set7 = temp2;
            else
                set8 = temp2;
            end
            
            Violin( temp2, pos(aa), 'Bandwidth', 0.2,  'MedianColor', [1,1,1], 'ViolinAlpha', .5, 'BoxColor',[0 0 0], 'ShowData', false, 'Width', .32, ...
                'BoxWidth', 0.1, 'ShowMean', true);
            set(gca, 'XLim', [0.5, 4.5], "YLim", [-0.1, 5]);
            xticks([1 2 3 4 5]);
        end
    end
    hold off
end