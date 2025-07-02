function [Model,tit]=getScore(Ms1,Ms2,pp,tt,nj,dt,tspan)
Par=Params;
M=Par(32);

%%TNF
w=55000;                    %Molecular weight g/mol
ConDose=10;                 %Concentration 10 ng/ml
MLigandTNF=ConDose/w;       %uM

%%LPS
w=1800;                    %Molecular weight g/mol
ConDose=1;                 %Concentration 1 ng/ml
MLigandLPS=ConDose/w;       %uM

%%PolyIC
w=150000;                     %Molecular weight g/mol
ConDose=1000;              %Concentration 1000 ng/ml extracellular
MLigandPIC=ConDose/w;       %uM

%%CpG
ConDose=300;                 %Concentration 300 nM extracellular
MLigandCpG=ConDose*10^-3;       %uM

rng default
Model={};
nn=26;

tit=["1) Ctrl->TNF","2) PIC(TNF-IFN-KO)->Ctrl","3) PIC(TNF-IFN-KO)->TNF", "4) PIC(TNF-KO+sIFNaR)->TNF",...
     "5) PIC(TNF-KO+ParaIFN)->TNF","6) PIC(TNF-KO)->TNF","7) PIC(sTNFR)->TNF","8) PIC(ParaTNF)->TNF",...
     "9) PIC(IFN-KO+sTNFR)->TNF","10) PIC(IFNa-KO)->TNF","11) PIC(sTNFR+sIFNaR)->TNF","12) PIC->TNF","13) PIC->Ctrl",...
     "14) Ctrl->LPS","15) LPS(TNF-IFN-KO)->Ctrl","16) CpG(TNF-IFN-KO)->LPS","17) CpG(TNF-KO+sIFNaR)->LPS",...
     "18) CpG(TNF-KO+Para-IFN)->LPS","19) CpG(TNF-KO)->LPS","20) CpG(sTNFR)->LPS","21) CpG(ParaTNF)->LPS",...
     "22) CpG(IFN-KO+sTNFR)->LPS","23) CpG(IFNa-KO)->LPS","24) CpG(sTNFR+sIFNaR)->LPS","25) CpG->LPS","26) CpG->Ctrl"];

%%%%%%%%%%%%%%% flag=1 (PolyIC), flag=0 CpG  %%%%%%%%%%%%%%%%%%%%
flag(1:13)=1;
flag(14:26)=0;

%%Autocrine IFN   1 2 3 4 5 6 7 8 9 10 11 12 13
scIFN =  Par(20)*[ 0 0 0 1 0 1 1 1 0 0  0  1  1   ...   %%PIC
                  0 0 0 1 0 1 1 1 0 0  0  1  1];          %%CPG

%%Paracrine IFN   1 2 3 4 5 6 7 8 9 10 11 12 13
IFNext = Par(16)*[0 0 0 0 1 1 1 1 0 0  0  1  1 ...        %%PIC
                  0 0 0 0 1 1 1 1 0 0  0  1  1 ];          %%CPG

%%Autocrine TNF   1 2 3 4 5 6 7 8 9 10 11 12 13
scTNF =  Par(19)*[0 0 0 0 0 0 1 0 1 1  1  1  1 ...        %%PIC
                  0 0 0 0 0 0 1 0 1 1  1  1  1 ];          %%CPG

%%Paracrine TNF   1 2 3 4 5 6 7 8 9 10 11 12 13
TNFext=  Par(16)*[0 0 0 0 0 0 0 1 0 1  0  1  1 ...    %%PIC
                  0 0 0 0 0 0 0 1 0 1  0  1  1 ];      %%CPG

%%First Challenge 1 2 3 4 5 6 7 8 9 10 11 12 13
Lc = MLigandPIC*[ 0 1 1 1 1 1 1 1 1 1  1  1  1];
Lc2 = MLigandCpG*[0 1 1 1 1 1 1 1 1 1  1  1  1];
Lc=[Lc Lc2];

alp=5;
bet=2000;
tnfExt=gampdf(tspan,alp,bet);
infExt=gampdf(tspan,alp,bet);

for i=1:nn
    yyT=[];
    
    m1T=[];
    m2T=[];
    m3T=[];
    m4T=[];
    m5T=[];

    IkbaT=[];
    A20T=[];
    
    ikkT=[];
    ISGT=[];
    IFNT=[];

    LT1=[];
    LT2=[];
    LT3=[];
    LT4=[];
    LT5=[];

    TT=[];

    Par(20)=scIFN(i);
    Par(19)=scTNF(i);

    parfor j=1:nj
        M1=Ms1(j);        %TLR3 Receptors
        M2=Ms2(j);        %TLR9 Receptors
        M3=M;             %IFN Receptors
        y0=zeros(23,1);

y0(1) = M;   %Receptors TNF
y0(2) = M;   %Receptors LPS

y0(13) = 0;   %Nuclear NFkB
y0(14) = 0;  %mRNA of target gene (uM)
y0(15) = 0.0001;   %Free Ikba (uM)
y0(16) = 0.0001;   %A20 (uM)
y0(17) = 0;    %Ligand TNF
y0(18) = 0;    %Ligand LPS
y0(19) = 0;    %Ligand PolyIC
y0(20) = 0;    %Ligand CpG
y0(21) = 0;    %Ligand IFNg

y0(19)=flag(i)*Lc(i);
y0(20)=(1-flag(i))*Lc(i);

    [t, y] = ode23tb(@(t,y) TNF_Model(t,y,dt,Par,M1,M2,M3,IFNext(i)*infExt,TNFext(i)*tnfExt), tspan, y0);
    T=t;
    YY=y;
  
    y0=y(end,:);
    if i~=2 && i~=13 && i~=15 && i~=26
    y0(17:20)=0;
    y0(17)=flag(i)*MLigandTNF;
    y0(18)=(1-flag(i))*MLigandLPS;
    end

    [t, y] = ode23tb(@(t,y) TNF_Model(t,y,dt,Par,M1,M2,M3,infExt*0,tnfExt*0), tspan, y0);
    T=[T; t+T(end)];
    YY=[YY; y];

    T= T/60;
    yy=YY(:,13);
    
    m1=YY(:,1)+YY(:,3);
    m2=YY(:,2)+YY(:,4);
    m3=YY(:,5);
    m4=YY(:,6);
    m5=YY(:,7);

    A1=YY(:,15);
    A2=YY(:,16);
    
    ikk=YY(:,11);
    ISG=YY(:,22);
    IFN=YY(:,23);

    L1=YY(:,17);
    L2=YY(:,18);
    L3=YY(:,19);
    L4=YY(:,20);
    L5=YY(:,21);


    yyT=[yyT; yy'];
    m1T=[m1T; m1'];
    m2T=[m2T; m2'];
    m3T=[m3T; m3'];
    m4T=[m4T; m4'];
    m5T=[m5T; m5'];

    IkbaT=[IkbaT; A1'];
    A20T=[A20T; A2'];

    ikkT=[ikkT; ikk'];
    ISGT=[ISGT; ISG'];
    IFNT=[IFNT; IFN'];

    LT1=[LT1; L1'];
    LT2=[LT2; L2'];
    LT3=[LT3; L3'];
    LT4=[LT4; L4'];
    LT5=[LT5; L5'];

    TT=[TT; T'];
    end

    tt2=minutes(TT(1,:));
    vnames=["NfkB","TNFR","TLR4","TLR3","TLR9","IFNR","IKK","IRF7","STAT","Ikba","A20",...
        "TNF","LPS","PolyIC","CpG","IFNb"];
    mod=timetable(yyT',m1T',m2T',m3T',m4T',m5T',ikkT',ISGT',IFNT',IkbaT',A20T',LT1',LT2',LT3',LT4',LT5','RowTimes',tt2','VariableNames',vnames);
    mod=synchronize(mod,minutes(tt),"spline");
        Model{i}=mod;

    if pp==1
        subplot(2,nn/2,i)
        plot(tt2,yyT,'LineWidth',2)
        hold off
        ylim([0 1])
        xlim(minutes([-10 500]))
        legend({"NFkB"},"Location","best")
        title(tit(i))
    elseif pp==2
        subplot(2,nn/2,i)
        auc=trapz(yyT(:,1:20),2);
        [~,idx]=sort(auc,"descend");
        imagesc(tt2,1:nj,yyT(idx,:),[0 1])
        xlim(minutes([-10 500]))
        ylabel([]);
        yticks([]);
        title(tit(i))
        set(gca,"FontSize",8)
    end

end

drawnow
end
