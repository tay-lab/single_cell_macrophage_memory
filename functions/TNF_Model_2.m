% ODE right-hand side of the NF-kB / IRF / STAT model (Methods, Eqs. 1-24). Used by getScore and getScoreCondTot.
% Source: NFkB Project Andrew\Simulation\20260414 (Apr 2026); copied 2026-10-18 without code changes.
function dy=TNF_Model_2(t,y,dt,Par,M,M2,IFNext,TNFext)
% flag=1 (PolyIC), flag=0 CpG
    [knin, Nt, Ki, klin, Kn, kt,gamma,ktl,a,ka20,IKKt,Tdeg,Kca,kb,kd,ka,ki,kp,n1,n2,kdeg,kr,ks,kb3,kb4] = Params;

    tlIFN=Par(1);
    tlTnf=Par(2);

    infExt=IFNext(floor(t/dt));
    tnfExt=TNFext(floor(t/dt));

    % kb=0;
    Nt=Par(3);
    IKKt=Par(4);
    Kca=Par(5);
    

    % gamma=Par(4);
    % kdeg=Par(5);
    % % kb3=Par(6);
    % % kb4=Par(7);
    % % kb=Par(8);
    uM2m=(6.022e17*12e-15);
    
dy=zeros(23,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Receptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(1) = -kb*y(1)*(y(17)+tnfExt)+kd*y(3)+ ktl*y(23)*uM2m;            %Receptors TNFR   (#)
dy(2) = -kb*y(2)*y(18)+kd*y(4);%+ ktlR*y(23);                       %Receptors TLR4   (#)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Active Receptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(3)=kb*y(1)*(y(17)+tnfExt)-kd*y(3);                                    %Active Receptors TNFR   (#)
dy(4)=kb*y(2)*y(18)-kd*y(4);                                             %Active Receptors TLR4   (#)
dy(5)=kb3*(M2-y(5))*y(19)-kd*y(5);                                       %Active Receptors TLR3 (#)
dy(6)=kb4*(M2-y(6))*y(20)-kd*y(6);                                       %Active Receptors TLR9 (#)
dy(7)=kb*(M-y(7))*(y(21)+infExt)-kd*y(7);                                %Active Receptors IFNAR  (#)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(8)=kr*(y(3)+y(4)+y(5)+y(6))*(Kca-y(8))*(ka20^n2/(ka20^n2+y(16)^n2))-ks*y(8);     %Active Adaptor IKKK (uM)
dy(9)=kr*(y(5)+y(6))*(Kca-y(9))-ks*y(9)-knin*y(9)^2;                                %Active IRF      (uM)
dy(10)=kr*y(7)*(Kca-y(10))-ks*y(10)-knin*y(10)^2;                                   %Active STAT      (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IKK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(11) = ka*abs(y(8)^n1/(y(8)^n1+.5^n1))*(IKKt - y(11) - y(12))*(ka20^n2/(ka20^n2+y(16)^n2)) - ki*y(11);   %Active iKKa    (uM)
dy(12) = ki*y(11)- kp*y(12);                   %Inactive iKKa  (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NFkB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(13) = knin*y(11)*(Nt-y(13))*(Ki/(Ki+y(15))) - klin*y(15)*(y(13)/(y(13)+ Kn));   %Nuclear NFkB (uM) 
dy(14) = kt*y(13).^2 - gamma*y(14);                                                %mRNA of target gene (uM)
dy(15) = ktl*y(14) - a*y(11)*(Nt - y(13))*(y(15)/(Ki+y(15)));                      %Free Ikba (uM)
dy(16) = ktl*y(14)-kdeg*y(16);                                                     %A20 (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ligands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(17) = -Tdeg*y(17) + tlTnf*y(14);                                 %Ligand TNF (uM)
dy(18) = -Tdeg*y(18);                                               %Ligand LPS (uM)
dy(19) = -Tdeg*y(19);                                               %Ligand PolyIC (uM)
dy(20) = -Tdeg*y(20);                                               %Ligand CpG (uM)
dy(21) = -Tdeg*y(21) + tlIFN*y(22);                                 %Ligand IFNAR (uM)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(22)= knin*y(9)^2-klin*y(22);                                      %Nuclear IRF (uM)
dy(23)= knin*y(10)^2-klin*y(23);                                    %Nuclear STAT (uM)

end
