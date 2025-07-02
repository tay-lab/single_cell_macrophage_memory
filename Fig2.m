%% Generates Figure 2
addpath functions\
load genExtFeat.mat     %Load Dataclear

%% Fig 2A-D
figure
[fE2,fL2]=cleanfeat(MIfeatE,MIfeatL,fE,fL,datvioE,datvioL,0,1);
zz2=data0.Category==gname{2};
data=data0(zz2,:);

kmax=5;
col=parula(kmax);
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
legend

subplot(4,7,[6,28])
imagesc(tt(1,:),1:size(dd,1),dd,[0 6])
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
    violin({vioTab{:}},'facecolor',col,'mc',[],'medc','black');
     hold on
        plot(mkE,'o')
    title(fE2(i))
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
    violin({vioTab{:}},'facecolor',col,'mc',[],'medc','black');
    hold on
        plot(mkE,'o')
    title(fL2(i))
    set(gca,'FontSize',8,'ylim',[0 Inf])
end

%% LSTM model
rng default
seqA=data.("Phase1 NfkB");
seqB=data.("Phase2 NfkB");
hpart = cvpartition(size(seqA,1),'Holdout',.2); % Nonstratified partition

ratTrain=seqA(hpart.training,:);
ratTest=seqA(hpart.test,:);

yTrain=seqB(hpart.training,:);
yTest=seqB(hpart.test,:);

clsTTest=clsT(hpart.test);

numFeatures = size(seqA,2);
numHiddenUnits = 128;
numResponses = size(seqB,2);
maxEpochs = 200;
miniBatchSize = 60;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits,OutputMode="sequence")
    selfAttentionLayer(8,64,"Name","selfattention")
    fullyConnectedLayer(numResponses)];

options = trainingOptions("adam", ...
    MaxEpochs=maxEpochs, ...
    MiniBatchSize=miniBatchSize, ...
    InitialLearnRate=0.001, ...
    GradientThreshold=1, ...
    Shuffle="every-epoch", ...
    Metrics="rmse", ...
    Plots="training-progress", ...
    Verbose=0);

net = trainnet(num2cell(ratTrain,2),num2cell(yTrain,2),layers,"mse",options);
YPredict = cell2mat(minibatchpredict(net,num2cell(ratTest,2),MiniBatchSize=1,UniformOutput=false));

%% Fig 2G 
AA=[ratTest yTest];
BB=[ratTest YPredict];
figure 

auc=trapz(BB(:,1:10),2);
[~,idx]=sort(auc,'descend');
imagesc(BB(idx,:),[0 5])
xlabel("Time (min)")
title("Predicted NFkB traces in CpG-> LPS")

%% Fig 2H
figure
    ftime=[0:79]*6;

    plot(ftime,median(BB),"b",'LineWidth',2)
    hold on

    plot(ftime,median(AA),"--",'LineWidth',2)
    p25 = prctile(AA,25);
    p75 = prctile(AA,75);
    hIQR = fill([ftime(41:end),fliplr(ftime(41:end))],[p75(41:end),fliplr(p25(41:end))],"r","FaceAlpha",.1);

    ylabel("NFkB dynamics")
    xlabel("Time step")
    legend(["Test Data" "Predicted"],Location="best")
    ylim([0 3])
    
%% Fig 2I
idx=[85    69     5    53];
for i = 1:numel(idx)
    subplot(2,2,i)
    plot([0:79]*6-6,AA(idx(i),:),"LineWidth",2)
    hold on
    plot([40:79]*6-6,BB(idx(i),41:end),"--","LineWidth",2)
    hold off
    
    xlim([0 79]*6)
    title("Test Observation " + idx(i))
    ylim([-.5 6])
    xlabel("Time Step")
    set(gca,'FontSize',10,'FontName','sans serif')
end
legend(["Test Data" "Predicted"],Location="best")