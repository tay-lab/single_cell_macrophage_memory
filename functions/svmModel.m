function [accuracy,datP]=svmModel(data,fE2,k)
%%%%%%%%%%%%%%%%%%%
rng(2);  % Set global seed (optional, for reproducibility outside parfor)
nc=size(data,1);
r=randperm(nc,nc);
data=data(r,:);
labels=categorical(data.Cluster);
predictors=zscore(table2array(data(:,fE2)));

template = templateLinear(...
    'Learner', 'SVM', ...
    'Lambda', 'auto', ...
    'BetaTolerance', 0.0001);

foldSize = floor(nc / k);

datP={};
parfor fold = 1:k
    fold
    % Indices for test
    TestIdx = (1 + (fold-1)*foldSize) : fold*foldSize;
    % Ensure indices are within bounds
    TestIdx(TestIdx > nc) = [];
    trainIdx = setdiff(1:nc, TestIdx);
    % dataTest=data(TestIdx,:);

    XTrain = predictors(trainIdx,:);
    YTrain = labels(trainIdx,:);
    XTest = predictors(TestIdx,:);
    YTest = labels(TestIdx);

    Mdl = fitcecoc(XTrain,YTrain,...
    'Learners', template, ...
    'Coding', 'onevsone');

    [~,YPred] = predict(Mdl, XTest);


    datP{fold}={Mdl,XTest,YTest,YPred};
end
% 
    classNames=[{'1'} {'2'} {'3'} {'4'} {'5'}];
    Cs=[];
    for fold=1:k
        YTest=datP{fold}{3};
        YPred=datP{fold}{4};
        classNames=categories(YTest);
        YTestSmallNet = scores2label(YPred,classNames);
    
        Cs(:,:,fold)=confusionmat(YTest,YTestSmallNet );
    end



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

figure("Name","SVM Model")
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