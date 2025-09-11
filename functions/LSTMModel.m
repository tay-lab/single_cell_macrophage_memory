function [accuracy,datP]=LSTMModel(data,fE2,k)
%%%%%%%%%%%%%%%%%%%
rng default;  % Set global seed (optional, for reproducibility outside parfor)
stream = RandStream('CombRecursive', 'Seed', 2);
numIterations = 5;
% Create one stream per worker using a cell array
streamCell = cell(numIterations, 1);
for i = 1:numIterations
    streamCell{i} = RandStream('CombRecursive', 'Seed', 2 + i);
end
streamConst = parallel.pool.Constant(streamCell);
%%%%%%%%%%%%%%%%%%%%%%%

nc=size(data,1);
r=randperm(nc,nc);
data=data(r,:);
labels=categorical(data.Cluster);
predictors=table2array(data(:,fE2));

numFeatures = 1;
numHiddenUnits = 128;
numResponses = 5;
maxEpochs = 100;
miniBatchSize = 64;

layers = [ ...
sequenceInputLayer(numFeatures)
lstmLayer(numHiddenUnits,OutputMode="last")
dropoutLayer(0.5)
selfAttentionLayer(8,64,"Name","selfattention")
fullyConnectedLayer(numResponses)
softmaxLayer
];

foldSize = floor(nc / k);

datP={};

parfor fold = 1:k
RandStream.setGlobalStream(streamConst.Value{fold});  % Set stream for this iteration
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

    muT=mean(XTrain ,"all");
    sigT=std(XTrain,0,"all");

    XTrain=(XTrain-muT)./sigT;
    XTest=(XTest-muT)./sigT;

    XTrain=num2cell(XTrain',1)';
    XTest=num2cell(XTest',1)';

    options = trainingOptions("adam", ...
    MaxEpochs=maxEpochs, ...
    MiniBatchSize=miniBatchSize, ...
    InitialLearnRate=0.001, ...
    GradientThreshold=1, ...
    Metrics = ["accuracy"],  ...
    Shuffle="never", ...
    Plots="none", ... %%"training-progress"
    Verbose=0);

    net = trainnet(XTrain,YTrain,layers,"crossentropy",options);
    YPred = minibatchpredict(net,XTest);
    XTest=sigT.*cell2mat(XTest')'+muT;

    datP{fold}={net,XTest,YTest,YPred};
end
% 
 % classNames=[{'1'} {'2'} {'3'} {'4'} {'5'}];
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

figure("Name","LSTM Model")
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