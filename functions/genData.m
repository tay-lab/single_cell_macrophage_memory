% This function read the raw NFkB traces stored in the file RawTraces.mat 
% and generates a table with all the extracted features in Phase1 and Phase2. 
% It saves the generated data in genExtFeat.mat

load RawTraces.mat
groups={[1:2] [3:16] [17:20] [21:41]};
gname={'Control -> LPSnaive' 'CpG -> LPScond' 'Control -> TNFnaive' 'PolyIC -> TNFcond'};

[data0,tAct]=genData(groups,gname,R);

save genExtFeat data0 tAct gname groups


function [data,threshActivation]=genData(groups,gname,R)
    aP=[];
    idelem=[];
    Rtime=[];
    threshActivation=[];
    mmx=[];
for n=1:length(groups)
    temp=[];
    tt=[];
    for aa=groups{n}
            temp =[temp; R{aa,1}];
    end
    
    ft=0:size(R{1,1},2)-1;
    ncells=size(temp,1);
    a2plot=ones(ncells,size(R{1,1},2));
    ftime=ones(ncells,size(R{1,1},2));
    aw0=[];
    
    for dd=1:ncells
        a2plot(dd,:)=temp(dd,:);
        ftime(dd,:)=ft*6;
        if anynan(a2plot(dd,:))
            aw0=[aw0 dd];
        end
    end
    a2plot(aw0,:)=[];
    ftime(aw0,:)=[];
  
    a2plot=smoothdata(a2plot',"lowess",[2 2]);

    if n==1 | n==3  
        marr=max(a2plot,2);
        threshActivation=[threshActivation; ...
            mean(a2plot(10:40,:),"all")+3*std(a2plot(10:40,:),[],"all")]
    end
    mmx=[mmx marr];
    nc=size(a2plot,2);

    aP=[aP a2plot];
    idelem=[idelem; n*ones(nc,1)];
    Rtime=[Rtime ftime'];
end

T=table(aP');
T.RTime=Rtime';
T.id=idelem;
T=renamevars(T,"Var1","Activation");
ds = arrayDatastore(T,"OutputType","same");

actid=readall(ds,"UseParallel",true).id';
features=transform(ds,@(x) extract(x,threshActivation));
data=readall(features,"UseParallel",true);
actioncats = categorical(gname)';
data.Category = actioncats(actid);

cat1={};
cat2={};
for i=1:length(gname)
    tit=gname{i};
    ss=strfind(tit,' -> ');
    cat1=[cat1 {tit(1:ss-1)}];
    cat2=[cat2 {tit(ss+4:end)}];
end

actioncats1 = categorical(cat1)';
actioncats2 = categorical(cat2)';

data.CategoryL1 = actioncats1(actid);
data.CategoryL2 = actioncats2(actid);

end