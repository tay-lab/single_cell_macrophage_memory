function Par = Params

Par(1) =.1;         %k1b Receptor Activation 1
Par(2) =.001;       %k2b Receptor Activation 2
Par(3) =.0001;      %k3b Receptor Activation 3
Par(4) =.001;       %kd  Receptor Inactivation

Par(5)=5*10^-5;     %kr Adaptor activation 
Par(6)=.01;         %ks Adaptor inactivation 
Par(7)=.0002;       %ka IKK neutral->active 
Par(8)=.0005;       %ki IKK active->inactive 
Par(9)=.002;        %kp IKK inactive->neutral

Par(10) = .01;         %NFkB, STAT nuclear import
Par(11) = .05;         %NFkB, STAT nuclear export
Par(12) = .0001;       %IRF nuclear import
Par(13) = .0005;       %IRF nuclear export

Par(14)= 1.4*10^-7;    %kt Transcription rate
Par(15)=.004;          %gamma mRNA degradation      
Par(16)=.5;            %ktl Translation rate      
Par(17)=1*10^-5;       %NFB degradation
Par(18)=0.00175;       %IkBa degradation

Par(19)=.001;          %sc TNF translation
Par(20)=5*10^-5;       %sc IFN translation
Par(21)=150;           %TNFR translation

Par(22)=2*10^-4;       %Ligand Degradation

Par(23)=0.018;         %Coefficient ka20 in hill function
Par(24)=1;             %Coefficient Kca in hill function              
Par(25)=2;             %Total IKK
Par(26)=0.029;         %Coefficient Kn in hill function
Par(27)=0.035;         %Coefficient Ki in hill function

Par(28)=1;             %Total NFkB
Par(29)=2;             %n1
Par(30)=2;             %n2
Par(31)=.5;            %Coefficient Ka in hill function
Par(32)=2000;

end