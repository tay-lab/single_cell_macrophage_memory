function [accuracy,weights]=randomModel(data,k)
rng default;  % Set global seed (optional, for reproducibility outside parfor)
stream = RandStream('CombRecursive', 'Seed', 1);

numIterations = 5;
% Create one stream per worker using a cell array
streamCell = cell(numIterations, 1);
for i = 1:numIterations
    streamCell{i} = RandStream('CombRecursive', 'Seed',100+i);
end
streamConst = parallel.pool.Constant(streamCell);

nc=size(data,1);
r=randperm(nc,nc);
data=data(r,:);
labels=categorical(data.Cluster);
Cs=[];

classes=1:5;
weights=countlabels(labels).Percent/100;

foldSize = floor(nc / k);

parfor fold = 1:k
    RandStream.setGlobalStream(streamConst.Value{fold});  % Set stream for this iteration
    fold
    % Indices for test
    TestIdx = (1 + (fold-1)*foldSize) : fold*foldSize;
    % Ensure indices are within bounds
    TestIdx(TestIdx > nc) = [];
    trainIdx = setdiff(1:nc, TestIdx);
    dataTest=data(TestIdx,:);

    YTest = labels(TestIdx);
    YPred = categorical(randsample(classes, length(YTest), true, weights)');
    Cs(:,:,fold)=confusionmat(YTest,YPred );
end
% 

numClasses=5;
precision = zeros(numClasses,k, 1);
recall = zeros(numClasses,k, 1);
f1_score = zeros(numClasses,k, 1);
accuracy=zeros(k, 1);
for fold=1:k
    ConfMatrix=Cs(:,:,fold);
correctPredictions = sum(diag(ConfMatrix));
totalSamples = sum(ConfMatrix(:));
accuracy(fold) = correctPredictions / totalSamples;

for c = 1:numClasses
    TP = ConfMatrix(c, c);
    FP = sum(ConfMatrix(:, c)) - TP;
    FN = sum(ConfMatrix(c, :)) - TP;

    if TP + FP == 0
        precision(c,fold) = NaN;
    else
        precision(c,fold) = TP / (TP + FP);
    end

    if TP + FN == 0
        recall(c,fold) = NaN;
    else
        recall(c,fold) = TP / (TP + FN);
    end

    if isnan(precision(c,fold)) || isnan(recall(c,fold)) || (precision(c,fold) + recall(c,fold) == 0)
        f1_score(c,fold) = NaN;
    else
        f1_score(c,fold) = 2 * (precision(c,fold) * recall(c,fold)) / (precision(c,fold) + recall(c,fold));
    end

end

end

figure("Name","Random Model")
subplot(1,2,1)
boxplot(f1_score','Widths',0.5)
hold on
scatter([1:k],f1_score,300,".","blue")
ylabel("F1 score")
xlabel("Group")
ylim([0 .7])
set(gca,"FontSize",12)

subplot(1,2,2)
totalConfMatrix=sum(Cs,3);
cm=totalConfMatrix./sum(totalConfMatrix);
h=heatmap(round(100*cm,1));
h.CellLabelFormat = '%.1f%%'; % Displays one decimal place and a percent sign
h.Title = 'Group classification';
ylabel("True Group")
h.ColorbarVisible = 'off';
set(gca,"FontSize",10)
% h.XDisplayLabels = {'','','','',''};
xlabel("Predicted Group")
end