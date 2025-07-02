function dy=TNF_Model(t,y,dt,Par,M1,M2,M3,IFNext,TNFext)
    k1b=Par(1);
    k2b=Par(2);
    k3b=Par(3);
    kd=Par(4);

    kr=Par(5);
    ks=Par(6);
    ka=Par(7);
    ki=Par(8);
    kp = Par(9);

    knin=Par(10); 
    klin=Par(11);  
    kpnin=Par(12); 
    kplin=Par(13); 

    kt=Par(14); 
    gamma=Par(15);
    ktl=Par(16);
    kdeg=Par(17);
    a=Par(18);

    tlTnf=Par(19);
    tlIFN=Par(20);
    ktlR=Par(21);
    Tdeg=Par(22);

    ka20=Par(23);
    Kca=Par(24);
    IKKt=Par(25);
    Kn=Par(26);
    Ki=Par(27);
    Nt=Par(28);

    n1=Par(29);
    n2=Par(30);
    Ka=Par(31);

    infExt=IFNext(floor(t/dt));
    tnfExt=TNFext(floor(t/dt));

dy=zeros(23,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Receptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(1) = -k1b*y(1)*(y(17)+tnfExt)+kd*y(3)+ ktlR*y(23);    %Receptors TNFR   (#)
dy(2) = -k1b*y(2)*y(18)+kd*y(4)+ ktlR*y(23);    %Receptors TLR4   (#)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Active Receptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(3)=k1b*y(1)*(y(17)+tnfExt)-kd*y(3);                                    %Active Receptors TNFR   (#)
dy(4)=k1b*y(2)*y(18)-kd*y(4);                                    %Active Receptors TLR4   (#)
dy(5)=k2b*(M1-y(5))*y(19)-kd*y(5);                         %Active Receptors TLR3 (#)
dy(6)=k3b*(M2-y(6))*y(20)-kd*y(6);                         %Active Receptors TLR9 (#)
dy(7)=k1b*(M3-y(7))*(y(21)+infExt)-kd*y(7);                               %Active Receptors IFNAR  (#)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(8)=kr*(y(3)+y(4)+y(5)+y(6))*(Kca-y(8))*(ka20^n2/(ka20^n2+y(16)^n2))-ks*y(8);     %Active Adaptor IKKK (uM)
dy(9)=kr*(y(5)+y(6))*(Kca-y(9))-ks*y(9)-kpnin*y(9)^2;             %Active IRF      (uM)
dy(10)=kr*y(7)*(Kca-y(10))-ks*y(10)-knin*y(10)^2;                        %Active STAT      (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IKK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(11) = ka*(y(8)^n1/(y(8)^n1+Ka^n1))*(IKKt - y(11) - y(12))*(ka20^n2/(ka20^n2+y(16)^n2)) - ki*y(11);   %Active iKKa    (uM)
dy(12) = ki*y(11)- kp*y(12);                   %Inactive iKKa  (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NFkB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(13) = knin*y(11)*(Nt-y(13))*(Ki/(Ki+y(15))) - klin*y(15)*(y(13)/(y(13)+ Kn));   %Nuclear NFkB (uM) 
dy(14) = kt*y(13).^2 - gamma*y(14);                                                %mRNA of target gene (uM)
dy(15) = ktl*y(14) - a*y(11)*(Nt - y(13))*(y(15)/(Ki+y(15)));                      %Free Ikba (uM)
dy(16) = ktl*y(14)-kdeg*y(16);                                                     %A20 (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ligands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(17) = -Tdeg*y(17) + tlTnf*y(14);                                 %Ligand TNF          (uM)
dy(18) = -Tdeg*y(18);                                             %Ligand LPS          (uM)
dy(19) = -Tdeg*y(19);                                             %Ligand PolyIC       (uM)
dy(20) = -Tdeg*y(20);                                             %Ligand CpG          (uM)
dy(21) = -Tdeg*y(21) + tlIFN*y(22);                              %Ligand IFNAR        (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(22)=kpnin*y(9)^2-kplin*y(22);                               %Nuclear IRF    (uM)
dy(23)= knin*y(10)^2-klin*y(23);                                     %Nuclear STAT    (uM)

end
