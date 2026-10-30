% Simulates one synthetic cell for the conditioned protocols (polyI:C -> TNF, CpG -> LPS). Used by FigS6.
function [s,Nf]=getScoreCondTot(Par,pp,rpkE,myExp2,M2,yin)
options = odeset('NonNegative',1:18);

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

rng(0)
s=0;
Model={};

dt=60;
tspan=dt:dt:240*60;
tt=0:6:474;

signalCyt(1,:)=[MLigandPIC 0 MLigandTNF 0];
signalCyt(2,:)=[0 MLigandCpG 0 MLigandLPS];

IFNext=zeros(1,length(tspan));
TNFext=zeros(1,length(tspan));

for i=1:2

y0=yin;

% y0(1) = M;   %Receptors TNF
% y0(2) = M;   %Receptors LPS

% y0(13) = 0;   %Nuclear NFkB
% y0(14) = 0;  %mRNA of target gene (uM)
% y0(15) = 0.0001;   %Free Ikba (uM)
% y0(16) = 0.0001;   %A20 (uM)
% y0(17) = 0;    %Ligand TNF
% y0(18) = 0;    %Ligand LPS

    % [t, y] = ode23tb(@(t,y) TNF_Model_2(t,y,dt,Par,yin(1),M2,M3,IFNext,TNFext), 0:10:1000, y0,options);
    % T=t-1000-dt;
    % YY=y;
    % y0=y(end,:);

    y0(17)=0;
    y0(18)=0;
    y0(19)=signalCyt(i,1);
    y0(20)=signalCyt(i,2);

    T=0;
    YY=y0';

    [t, y] = ode23tb(@(t,y) TNF_Model_2(t,y,dt,Par,yin(1),M2,IFNext,TNFext), tspan, y0,options);
    T=[T; t];
    YY=[YY; y];
    % T=t;
    % YY=y;

    y0=y(end,:);
    y0(17)=signalCyt(i,3);
    y0(18)=signalCyt(i,4);
    % y0(17)=0;
    % y0(18)=0;
    y0(19)=0;
    y0(20)=0;

    [t, y] = ode23tb(@(t,y) TNF_Model_2(t,y,dt,Par,yin(1),M2,IFNext,TNFext), tspan, y0,options);
    T=[T; t+T(end)+dt];
    YY=[YY; y];

    T= T/60;
    yy=YY(:,13);
    nfkb=interp1(T,yy,tt);

    ikba=YY(:,15);
    ikba=interp1(T,ikba,tt);

    a20=YY(:,16);
    a20=interp1(T,a20,tt);

    tt1=tt(tt<240);
    tt2=tt(tt>240);
    [pk1,tm1]=max(nfkb(tt<240));
    [pk2,tm2]=max(nfkb(tt>240));

    if i==1
        % s=s+(pk2-rpkE(i,3)).^2;%+(pk1-rpkE(i,1)).^2;%+(tt1(tm2)-rpkE(i,5)).^2;
    else
        s=s+(pk2/pk1-rpkE(i,2)).^2;%+(tt1(tm1)-rpkE(i,4)).^2+(tt1(tm2)-rpkE(i,5)).^2;%+(pk1-rpkE(i,1)).^2;
    end
    % myExpN=[Fourier,mxS,tmS,speed];
    % Fourier

    % s=s+sum(abs(myExpN-myExp(i,:)));

    % s=s+(mxS-mxE).^2+(tmS-tmE).^2;
    % max(nfkb)
    if pp==1
        subplot(2,2,i)
        plot(tt,nfkb,'LineWidth',2)
        
        hold on
        plot(tt,myExp2(i+4,:),'LineWidth',2)
        hold off
        ylim([0 1])
        xlim([-10 480])
        legend({"NFkB"},"Location","best")

        subplot(2,2,i+2)
        plot(tt,ikba,'LineWidth',2)
        hold on
        plot(tt,a20,'LineWidth',2)
        hold off
        ylim([0 .1])
        xlim([-10 480])
        legend({"Ikba","A20"},"Location","best")
    end

    Nf(i,:)=nfkb;
end

drawnow

end
