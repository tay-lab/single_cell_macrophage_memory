%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2 related to 
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
%
% Wang AG, Mercado Vasquez JG, Son M, et al
%
% Corresponding author: tays@uchicago.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Load functions and NFkB traces
addpath functions\
load genExtFeat.mat     %Load Dataclear
%% Fig 2A-D
rng default
figure
[fE2,fL2]=cleanfeat(MIfeatE,MIfeatL,fE,fL,datvioE,datvioL,0,1);
zz2=data0.Category==gname{2};
data=data0(zz2,:);

col= [ 
0 114 189;
217 83 25;
237 177 32;
126 47 142;
119 172 48;
    ]/255;

kmax=5;
ratings=table2array(data(:,fL2));
opts = statset('MaxIter',1000);
[clsT,Cnt]=kmeans(ratings,kmax,"Replicates",10,"Distance","cityblock","Options",opts);
clusters=clsT;

mn=[];
for i=1:kmax
    mn=[mn median([data.(fL2(idxp(1)))(clusters==i,:)])];
end
[~,idx]=sort(mn,'ascend');

nCC=[];     %Number of cells in cluster i
Cntp=Cnt;
for i=1:kmax
    zz=clusters==idx(i);
    nCC=[nCC nnz(zz)];
    clsT(zz)=i;
    Cntp(i,:)=Cnt(idx(i),:);
end
Cnt=Cntp;

data0.Cluster(zz2)=clsT;
data.Cluster=clsT;

subplot(4,7,[1,9.2])
gscatter(ratings(:,1),ratings(:,2),clsT,col,'.....',[16 16 16 16 16])
xlabel(fL2(1))
ylabel(fL2(2))
set(gca,'FontSize',10)

Nf=[data.("Phase1 NfkB") data.("Phase2 NfkB")];
tt=[data.("Phase1 Time") data.("Phase2 Time")];
dd=[];
nc=[];

subplot(4,7,[3.8,12])
hold on
for i=1:kmax
    zz=clsT==i;
    nf=Nf(zz,:);
    nc=[nc size(nf,1)];
    [~,idxs]=sort(trapz(nf(:,1:10),2),"descend");
   dd=[dd; nf(idxs,:)];
   plot(tt(1,:),mean(nf),"LineWidth",2,"Color",col(i,:))
end
xticks([0:120:480])
xlabel("Time (min)")
ylabel("NFkB Dynamics")
xlim([-10 480])
ylim([-.1 4])
set(gca,'FontSize',10)
hold off
legend(["1","2","3","4","5"])

subplot(4,7,[6,28])
imagesc(tt(1,:),1:size(dd,1),dd,[0 4])
yticks(cumsum(nc))
ylabel([])
set(gca, 'Xlim', [0 480]);
xticks([-10,0:120:480])
xlabel('Time (min)')
set(gca,'FontSize',10)

for i=1:4
    subplot(5,6,6*3+i)
    vioTab={};
    mpx=[];
    vT=table2array(data(clsT==1,fE2(i)));
    mvT=median(vT,"omitnan");
    mkE=[];
    for j=1:kmax
        vioTab{j}=table2array(data(clsT==j,fE2(i)));
        vv=vioTab{j};
        vv(vv<0)=0;
        vioTab{j}=vv;
        mkE(j)=median(vioTab{j},"omitnan");
        mpx=[mpx prctile(vioTab{j},96)];
    end
    violin3({vioTab{:}},'facecolor',col,'mc',[],'medc','black');
    ylabel(fE2(i))
    % set(gca,'FontSize',8,'ylim',[0 max(mpx)])
    set(gca,'FontSize',8,'ylim',[0 Inf])

    subplot(5,6,6*4+i)
    vioTab={};
    mpx=[];
    vT=table2array(data(clsT==1,fL2(i)));
    mvT=median(vT,"omitnan");
    mkE=[];
    for j=1:kmax
        vioTab{j}=table2array(data(clsT==j,fL2(i)));
        vv=vioTab{j};
        vv(vv<0)=0;
        vioTab{j}=vv;
        mkE(j)=median(vioTab{j},"omitnan");
        mpx=[mpx prctile(vioTab{j},99)];
    end
    violin3({vioTab{:}},'facecolor',col,'mc',[],'medc','black');
    ylabel(fL2(i))
    set(gca,'FontSize',8,'ylim',[0 Inf])
end
% clc
% summary(categorical(clsT))

%% LSTM model Categorical. Cross-Validation
[accuracyLSTM,datP]=LSTMModel(data,"Phase1 NfkB",5);
[accuracyLinear,~]=svmModel(data,fE2,5);
[accuracyTree,~]=treeModel(data,fE2,5);
[accuracyRandom,weights]=randomModel(data,5);
%%
classNames=[{'1'} {'2'} {'3'} {'4'} {'5'}];
k=5;
cvLabels = cell(k,1);
cvScores = cell(k,1);
for fold = 1:k
    cvLabels{fold} = datP{fold}{3};
    cvScores{fold} = datP{fold}{4};
    YTestSmallNet = scores2label(cvScores{fold},classNames);
    Cs(:,:,fold)=confusionmat(cvLabels{fold},YTestSmallNet);
end
rocObj = rocmetrics(cvLabels,cvScores,classNames);

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

figure("Name","LSTM Accuracy")
boxplot(f1_score','Widths',0.5)
hold on
scatter([1:k],f1_score,300,".","blue")
scatter([1:k],weights,300,"_","magenta","LineWidth",2)
ylabel("F1 score")
xlabel("Group")
ylim([0 .7])
set(gca,"FontSize",12)

figure
% subplot(3,3,[1 6])
totalConfMatrix=sum(Cs,3);
cm=totalConfMatrix./sum(totalConfMatrix);
h=heatmap(round(100*cm,1));
h.CellLabelFormat = '%.1f%%'; % Displays one decimal place and a percent sign
h.Title = 'Group classification';
ylabel("True Group")
h.ColorbarVisible = 'off';
set(gca,"FontSize",12)
% h.XDisplayLabels = {'','','','',''};
xlabel("Predicted Group")

% clim([0 76])
% 
% subplot(3,3,7:9)
% h=heatmap(round([100*diag(cm)'; 100*(1-diag(cm)')],1));
% xlabel("Predicted Group")
% h.ColorbarVisible = 'off';
% clim([0 76])

colormap sky

set(gca,"FontSize",12)

figure
boxplot([accuracyLSTM,accuracyLinear,accuracyTree,accuracyRandom],'Widths',0.5)
hold on
scatter([1:4],[accuracyLSTM,accuracyLinear,accuracyTree,accuracyRandom],300,".","blue")
xticklabels(["LSTM","SVM","Tree","Random"])
ylabel("Accuracy")
xlabel("Model")
set(gca,"FontSize",12)


figure("Name","LSTM ROC")
plot(rocObj,ShowConfidenceIntervals=true,ShowModelOperatingPoint=false)


% exportgraphics(gcf, 'myPlot.pdf', 'ContentType', 'vector');
