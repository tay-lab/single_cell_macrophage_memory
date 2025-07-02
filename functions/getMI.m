function [hh,propDen]=getMI(Rs,S,k,qq)
    rng default
    [zz,~]=find(~isfinite(Rs));
    Rs(zz,:)=[];
    S(zz,:)=[];
    % Rs=zscore(Rs);
    % noise=eps;
    noise=std(Rs,[],"all")/100;
    flag=true;
    while flag
        Rs=Rs+noise*(rand(size(Rs)));
        propDen=propDenKK(Rs,S,k);
        suminf=0;
        for i = 1:numel(propDen)
            suminf = suminf+sum(isinf(propDen{i}(:)));
        end
        if suminf==0
            flag = false;
        else
            noise=noise*10;
        end
    end
        nh=non_condEntropy(propDen,qq);
        cE=condEntropy(propDen,qq);
        % hh=max(0,nh-cE);
        hh=nh-cE;
end

function propDen=propDenKK(Rs,S,k)
    allS=unique(S);
    propDen=cell(length(allS),length(allS));
    for i=1:length(allS)
        zz1=S==allS(i);
        if sum(zz1)~=0
            for j=1:length(allS)
                zz2=S==allS(j);
                x=Rs(zz1,:);
                y=Rs(zz2,:);
                nj=size(y,1);
                d=size(y,2);
                vd=pi^(d/2)/gamma(d/2+1);
                if j==i || isequal(y,x)
                [~,Ds]=knnsearch(y,x,"K",k+1,"Distance","euclidean");
                dist=Ds(:,end);
                dist=dist.^d;
                pd=k./(nj*vd*dist);
                else
                [~,Ds]=knnsearch(y,x,"K",k,"Distance","euclidean");
                dist=Ds(:,end);
                dist=dist.^d;
                pd=k./(nj*vd*dist);
                end
                propDen{i,j}=pd;
            end
        end
    end
end


function HH=condEntropy(propDen,qq)
    hh=nan(length(propDen),1);
    for i=1:length(propDen)
        pd=propDen{i,i};
        pd(pd<=eps)=[]; 
        ni=length(pd);
        hh(i)=sum(log2(pd))/ni;
    end
    HH=-qq*hh;
end

function HH=non_condEntropy(propDen,qq)
    hh=nan(length(propDen),1);
    for i=1:length(propDen)
    ss=0;
    for w=1:length(propDen)
        pd=propDen{i,w};     
        ss=ss+qq(w)*pd;
    end
    ss=ss(ss>eps);
    ni=length(ss);
    hh(i)=sum(log2(ss))/ni;
    end
    HH=-qq*hh;
end