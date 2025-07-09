%% Generates Figure 1
addpath functions\
load genExtFeat.mat

%% Fig 1A
Nf=[data0.("Phase1 NfkB") data0.("Phase2 NfkB")];            %All Traces
ti=[data0.("Phase1 Time")(1,:) data0.("Phase2 Time")(1,:)];  %Time in traces

figure
for i=1:2
    subplot(1,2,i)
    zz=data0.Category==gname{2*i};
    nfkb=Nf(zz,:);
    auc=trapz(nfkb(:,1:10),2);

    [auc,idx]=sort(auc,"descend");
    imagesc(ti,1:nnz(zz),nfkb(idx,:),[0 4])
    xlabel("Time (min)")
    title(gname{2*i})
end

%% Fig 1B

MIcomb=zeros(length(ti),1);
k=10;
for jj=1:2
    zz1=data0.Category==gname{2*jj-1};
    zz2=data0.Category==gname{2*jj};
    nc1=size(data0(zz1,:),1);
    nc2=size(data0(zz2,:),1);
parfor ii=1:size(Nf,2)  
    C=Nf;
    Rs=[C(zz1,ii); C(zz2,ii)];
    S=categorical([ones(nc1,1); 2*ones(nc2,1)]);
    MIcomb(ii,jj)=getCC(Rs,S,k);
end
end

MIcomb(MIcomb<0)=0;

figure
plot(ti,MIcomb,'LineWidth',2)
% xregion(0,240)
% xregion(240,480,"FaceColor",'r')
xticks([240:60:480])
set(gca,'Xlim',[240 480])
xlabel("Time (min)")
ylabel("Mutual Information (bits)")
set(gca,'FontSize',14,'FontName','Times New Roman')
legend({"CpG -> LPS","PolyI:C -> TNF"})

zz1=data0.Category==gname{1} | data0.Category==gname{2};
zz2=data0.Category==gname{3} | data0.Category==gname{4};
ENfMI=zeros(size(data0,1),1);
LNfMI=zeros(size(data0,1),1);
[~,tp1]=findpeaks(MIcomb(:,1),"MinPeakProminence",.1);
[~,tp2]=findpeaks(MIcomb(:,2),"MinPeakProminence",.05);
ENfMI(zz1)=Nf(zz1,tp1(1));
ENfMI(zz2)=Nf(zz2,tp2(1));
if length(tp1)>1
LNfMI(zz1)=Nf(zz1,tp1(2));
else
    LNfMI(zz1)=Nf(zz1,tp1(1));
end
if length(tp2)>1
    LNfMI(zz2)=Nf(zz2,tp2(2));
else
    LNfMI(zz2)=Nf(zz2,tp1(1));
end
data0.("Phase1 NfkB MI")=ENfMI;
data0.("Phase2 NfkB MI")=LNfMI;
%% Fig 1C

Nf1=data0.("Phase2 NfkB");
Nf2=data0.("Phase2 NfkB");

g1=["LPSnaive","LPScond","TNFnaive","TNFcond"];

MI=[];
k=10;
for ii=1:length(g1)
for jj=1:length(g1)
    zz1=data0.CategoryL2==g1(ii);
    zz2=data0.CategoryL2==g1(jj);
    nc1=size(data0(zz1,:),1);
    nc2=size(data0(zz2,:),1);

    Rs=[Nf1(zz1,1:5:size(Nf1,2)); Nf2(zz2,1:5:size(Nf2,2))];
    S=categorical([ones(nc1,1); 2*ones(nc2,1)]);
    MI(ii,jj)=getCC(Rs,S,k);
end
end

figure
colormap sky
MIm=round(triu(MI,1),2);
MIm(MIm<=0)=NaN;
h=heatmap(MIm, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.XDisplayLabels=string({"LPS_{naive}","LPS_{cond}","TNF_{naive}","TNF_{cond}"});
h.YDisplayLabels=string({"LPS_{naive}","LPS_{cond}","TNF_{naive}","TNF_{cond}"});
h.ColorLimits = [0 1];
set(gca,'XLim',[2,4],'yLim',[1,3])
set(gca,'FontSize',18,'FontName','Times New Roman')
title("Ligand specificity")

%% Fig 1D

featurenames = [
    "Phase1 NfkB MI","Phase1 Activation","Phase1 Deactivation",...
        "Phase1 Speed",'Phase1 Duration','Phase1 Fourier',...
        "Phase1 Peak","Phase1 Time to Peak","Phase1 Prominence",...
        "Phase1 AUC","Phase1 EAUC","Phase1 LAUC",...
    "Phase2 NfkB MI","Phase2 Activation","Phase2 Deactivation",...
        "Phase2 Speed",'Phase2 Duration','Phase2 Fourier',...
        "Phase2 Peak","Phase2 Time to Peak","Phase2 Prominence",...
        "Phase2 AUC","Phase2 EAUC","Phase2 LAUC",...
        ];

fE=featurenames(1:12);
fL=featurenames(13:end);

MIfeatE=zeros(length(fL),1);
MIfeatL=zeros(length(fL),1);
k=10;

datvioL={};
datvioE={};
mkL={};
mkE={};
for jj=1:2
    zz1=data0.Category==gname{2*jj-1};
    zz2=data0.Category==gname{2*jj};
    nc1=size(data0(zz1,:),1);
    nc2=size(data0(zz2,:),1);
for ii=1:length(fL)
    datvioE{ii,jj,1}=data0.(fE(ii))(zz1,:);
    datvioE{ii,jj,2}=data0.(fE(ii))(zz2,:);
    mkE{ii,jj,1}=median(data0.(fE(ii))(zz1,:),"omitnan");
    mkE{ii,jj,2}=median(data0.(fE(ii))(zz2,:),"omitnan");
    Rs=[data0.(fE(ii))(zz1,:); data0.(fE(ii))(zz2,:)];
    S=[true(nc1,1); false(nc2,1)];
    
    MIfeatE(ii,jj)=getCC(Rs,S,k);

    datvioL{ii,jj,1}=data0.(fL{ii})(zz1,:);
    datvioL{ii,jj,2}=data0.(fL{ii})(zz2,:);
    mkL{ii,jj,1}=median(data0.(fL{ii})(zz1,:),"omitnan");
    mkL{ii,jj,2}=median(data0.(fL{ii})(zz2,:),"omitnan");
    Rs=[data0.(fL{ii})(zz1,:); data0.(fL{ii})(zz2,:)];
    S=[ones(nc1,1); 2*ones(nc2,1)];

    MIfeatL(ii,jj)=getCC(Rs,S,k);
end
end

MIfeatL(MIfeatL<0)=0;
MIfeatE(MIfeatE<0)=0;

fn=[ " NFkB MI_{max}"," Activation Time"," Deactivation Time"...
     " Speed",' Duration',' Oscillations',...
     " Amplitude"," Time to Peak"," Prominence",...
     " AUC"," Early AUC"," Late AUC",...
                ];
figure
subplot(1,2,1)
[~,idxp]=sort(MIfeatE(:,1),"descend");
bh=barh(MIfeatE(idxp,:),1.5);
set(gca,'YDir','reverse','Xlim',[0 1])
% set(gca,'FontSize',12,'FontName','Times New Roman')
legend({"CpG vs Ctrl","PIC vs Ctrl"},'Location','southeast')
yticks([])
ylabel("Features")
xlabel("Mutual Information (bits)")
title("Information storaged in NFkB features Phase 1")

ytips1 = bh(1).XEndPoints;
xtips1 = bh(1).YEndPoints;
text(xtips1,ytips1,fn(idxp),'HorizontalAlignment','left',...
    'VerticalAlignment','middle')

subplot(1,2,2)
[~,idxp]=sort(MIfeatL(:,1),"descend");
bh=barh(MIfeatL(idxp,:),1.5);
legend({"CpG -> LPS vs Ctrl -> LPS","PIC -> TNF vs Ctrl -> TNF"},'Location','southeast')
set(gca,'YDir','reverse','Xlim',[0 .4])
yticks([])
ylabel("Features")
xlabel("Mutual Information (bits)")
title("Information storaged in NFkB features Phase 2")
ytips1 = bh(1).XEndPoints;
xtips1 = bh(1).YEndPoints;
text(xtips1,ytips1,fn(idxp),'HorizontalAlignment','left',...
    'VerticalAlignment','middle')

%% Fig 1E-H and SubFig 2A-B
col=hsv(2);
xlb=[{"Ctrl","CpG"}; {"Ctrl","PolyIC"}]; 
for jj=1:2
    [~,idxp]=sort(MIfeatE(:,jj),"descend");
    figure("Name","Phase 1")
    for i=1:12
        subplot(4,3,i)
        violin3({datvioE{idxp(i),jj,1},datvioE{idxp(i),jj,2}},'facecolor',col,'mc',[],'medc','black','bw',mkE{idxp(i),jj,2}/3)
        rn=ranksum(datvioE{idxp(i),jj,1},datvioE{idxp(i),jj,2});
        hold on
        xticklabels(xlb(jj,:))
        ylabel(fn{idxp(i)})
        xticks([1,2])
        set(gca,'YLim',[0 Inf])
        set(gca,'FontSize',8,'FontName','Times New Roman')
        % text(1.5,2*mkE{idxp(i),jj,2},num2str(rn)) %If Wilcoxon
    end
end

col=hsv(2);
xlb=[{"LPS_{naive}","LPS_{cond}"}; {"TNF_{naive}","TNF_{cond}"}]; 
for jj=1:2
    [~,idxp]=sort(MIfeatL(:,jj),"descend");
    figure("Name","Phase 2")
    for i=1:12
        subplot(4,3,i)
        violin3({datvioL{idxp(i),jj,1},datvioL{idxp(i),jj,2}},'facecolor',col,'mc',[],'medc','black','bw',mkL{idxp(i),jj,2}/3)
        rn=ranksum(datvioL{idxp(i),jj,1},datvioL{idxp(i),jj,2});
        hold on
        xticklabels(xlb(jj,:))
        ylabel(fn{idxp(i)})
        xticks([1,2])
        set(gca,'YLim',[0 Inf])
        set(gca,'FontSize',8,'FontName','Times New Roman')
        % text(1.5,2*mkL{idxp(i),jj,2},num2str(rn)) %If Wilcoxon
    end
end
%% Save data
%% save genExtFeat.mat data0 tAct gname groups MIfeatE MIfeatL fE fL datvioE datvioL idxp
