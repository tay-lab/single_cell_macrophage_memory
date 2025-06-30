function feat=extract(data,threshActivation)
    t=data.RTime';
    Y=data.Activation';
    id=data.id;
    if id<=2 
        threshActivation=threshActivation(1);
    else
        threshActivation=threshActivation(2);
    end

    time2peak=NaN(2,1);
    peak=NaN(2,1);
    prominence=NaN(2,1);
    tlast=NaN(2,1);
    tfirst=NaN(2,1);
    speed=NaN(2,1);
    duration=zeros(2,1);
    AUC=NaN(2,1);
    Fourier=NaN(2,1);
    EAUC=NaN(2,1);
    LAUC=NaN(2,1);

    feat=table;

    z={};
    z{1}=t>=0 & t<240;
    z{2}=t>=240 & t<=480;
    
    tt1=t(z{1});
    YY1=Y(z{1});    
    tt2=t(z{2});
    YY2=Y(z{2});

    feat=addvars(feat,tt1','NewVariableNames',"Phase1 Time");
    feat=addvars(feat,tt2','NewVariableNames',"Phase2 Time");
    feat=addvars(feat,YY1','NewVariableNames',"Phase1 NfkB");
    feat=addvars(feat,YY2','NewVariableNames',"Phase2 NfkB");

    for i=1:2
        tt=t(z{i});
        tt=tt-tt(1);
        YY=Y(z{i});
        YY=YY-YY(1);
        dt=6;

        dt=dt*.1;
        xt=tt(1):dt:tt(end);
        YY = interp1(tt,YY,xt,'spline');
        tt=xt;

%%%%%%%%%% Speed %%%%%%%%%%%%%%
        dv=diff(YY)/dt;
        speed(i)=max(dv);

%%%%%%%%%% Duration %%%%%%%%%%%%%%
        zz=YY>=threshActivation;
        afirst=find(zz==1,1, 'first');
        if ~isempty(afirst)
            tfirst(i)=tt(afirst);
            alast=find(zz(afirst:end)==0,1, 'first');
            if ~isempty(alast)
                tlast(i)=tt(afirst-1+alast);
            else
                tlast(i)=tt(end);
            end
            duration(i)=nnz(zz)*dt;
        end
        
        AUC(i)=trapz(tt,YY);

        zz=tt<=60;
        EAUC(i)=trapz(tt(zz),YY(zz)); 
        LAUC(i)=AUC(i)-EAUC(i); 

%%%%%%%%%% Fourier %%%%%%%%%%%%%%
        L = length(tt);             % Length of signal
        Fs = 1/dt;            % Sampling frequency  
        f = Fs*(0:floor((L/2)))/L;
        ff = fft(YY);

        P2 = abs(ff/L);
        P1 = P2(1:floor(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        [~,locVal]=findpeaks(P1,f,'SortStr','descend','NPeaks',1);
        tf=1./locVal;
        Fourier(i)=tf;
        
%%%%%%%%%% Peaks %%%%%%%%%%%%%%
        [peakVal,locVal,~,p]=findpeaks(YY,tt,'SortStr','descend','NPeaks',1);
        if ~isempty(peakVal)
            peak(i)=peakVal;
            time2peak(i)=locVal;
            prominence(i)=p;
        end
    end
    
    featurenames = [
        "Phase1 EAUC","Phase1 LAUC","Phase1 AUC",...
        "Phase1 Speed",'Phase1 Duration','Phase1 Fourier'...
        "Phase1 Peak","Phase1 Time to Peak","Phase1 Prominence",...
        "Phase1 Activation","Phase1 Deactivation",...
        "Phase2 EAUC","Phase2 LAUC","Phase2 AUC",...
        "Phase2 Speed",'Phase2 Duration','Phase2 Fourier',...
        "Phase2 Peak","Phase2 Time to Peak","Phase2 Prominence",...
        "Phase2 Activation","Phase2 Deactivation"
        ];

    feat = addvars(feat, ...
        EAUC(1),LAUC(1),AUC(1),...
        speed(1),duration(1),Fourier(1,:),...
        peak(1,:),time2peak(1,:),prominence(1,:), ...
        tfirst(1),tlast(1),...
        EAUC(2),LAUC(2),AUC(2),...
        speed(2),duration(2),Fourier(2,:),...
        peak(2,:),time2peak(2,:),prominence(2,:), ...
        tfirst(2),tlast(2),...
        'NewVariableNames',featurenames);
end
