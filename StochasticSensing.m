function [Env R Hit]=StochasticSensing1(rL, rH, beta, p, l, N, method,Hx)
    %2022_06: Goal - represent simulated growth rates assuming optimal
    %decision-making.
    
    %Assume state variable m tracks number of H hits
    %Assume starting in low state so m/l = 0/l.
    
    %Inputs: rL:    growth in low state
    %       rH:     growth in high state
    %       beta:   benefit of committing to fixed state (beta>1)
    %       p:      probability of rH state
    %       l:      memory of population
    %       N:      Number of iterations
    %       method: Describes how Env is driven
    %           NaN => driven randomly with p
    %           vector => driven with artificial vector.  Vector must be of
    %           length N.
    %       Hx:     Indicates if desire to use deterministic prior Hx
    %           1 => adds beta*rL.*[1:...:1] to beginning of R.
    
    %Outputs:
    %       Env:    The 0/1 environment sequence
    %       R:      The temporal growth decision
    %       Hit:    The 0/1 indicator of decision accuracy.

%Calculate pI:
pI=(beta-1)*rL/((beta-1)*rL+rH);
    
if length(Hx)>1 %Explicit history
    m=Hx;
elseif length(Hx)==1
    m=zeros(1,l); %Start with constant 0 history.
else
    disp('problem categorizing Hx')
end
%Select environment:
Hit=nan(1,N);
if isnan(method(1))==1
    Env=rand(1,N)<p;
elseif isnan(method(1))==0
    Env=method;
end
R=nan(1,N);
for t=1:N
    %Current estimate
    m;
    if sum(m)/l<pI %Estimating RL
        R(t)=beta*rL*(1-Env(t));
        Hit(t)=Env(t)==0; %1 if Env is 0
    else %Estimating RH
        R(t)=rL*(1-Env(t)) + rH*Env(t);
        Hit(t)=Env(t)==1; %1 if Env is 1
    end
    %Update m, keeping l fixed.
    m=[m(2:end) Env(t)];
end
if length(Hx)==1
    %include average R from previous state
    RPrior=beta*rL*ones(1,l);
    R=[RPrior R];
    Env=[zeros(1,l) Env];
elseif length(Hx)>1
    RPrior=nan(1,length(Hx));
    for j=1:length(Hx)
        if Hx(j)==1
            RPrior(j)=rH;
        elseif Hx(j)==0
            RPrior(j)=beta*rL;
        else
            disp('Problem with Hx values')
        end
    end
    R=[RPrior R];
else
    disp('Problem categorizing Hx')
end

end