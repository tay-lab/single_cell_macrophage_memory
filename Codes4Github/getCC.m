function [I,qqmax] =getCC(Rs,S,k)
    dimq=length(unique(S));
    fun = @(qq) -getMI(Rs,S,k,qq);

    Aeq=ones(1,dimq);
    beq=1;
    lb =zeros(1,dimq);
    ub =ones(1,dimq);
    qq0=ones(1,dimq)/dimq;
    
    options = optimoptions('fmincon','Algorithm','interior-point','Display','none');
    [qqmax,I] = fmincon(fun,qq0,[],[],Aeq,beq,lb,ub,[],options);
    I=-I;
end

