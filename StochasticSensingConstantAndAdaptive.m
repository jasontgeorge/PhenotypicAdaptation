function [RConst HitConst RVariable HitVariable Env LSeq] =...
            StochasticSensing1ConstantAndAdaptiveLPIProxv2(...
                rL, rH, beta, p, l, lMin, lMax, N, method,Hx,adaptation)
 %2022_07: This changes memory based on proximity of p_est to pI.
    %State variable m tracks number of H hits
    %Starting in low state so m/l = 0/l.
    
    %Inputs: rL:    growth in low state
    %       rH:     growth in high state
    %       beta:   benefit of committing to fixed state (beta>1)
    %       p:      probability of rH state
    %       l:      INITIAL memory of population
    %       lMin:   Minimal memory size
    %       lMax:   Maximal memory size
    %       N:      Number of iterations
    %       method: Describes how Env is driven
    %           NaN => driven randomly with p
    %           vector => driven with artificial vector.  Vector must be of
    %           length N.
    %       Hx:     Indicates if desire to use deterministic prior Hx
    %           1 => history of 1s
    %           0 => history of 0s
    %       adaptation: 'Variance'  => sequential variance comparison
    %                   'Proximity' => proximity comparison
    
    %Outputs:
    %       Env:    The 0/1 environment sequence
    %       R:      The temporal growth decision
    %       Hit:    The 0/1 indicator of decision accuracy.
    %       LSeq:   Sequence of memory sizes at each time.

    l0=l;
    %Calculate pI:
    pI=(beta-1)*rL/((beta-1)*rL+rH);
    LSeq=nan(1,N); %LSeq(1)=l;
    
    if length(Hx)>1 %Explicit history
        m=Hx;
    elseif Hx==1
        m=ones(1,l); %Start with constant 1 history.
    elseif Hx==0
        m=zeros(1,l); %Start with constant 0 history.    
    else
        disp('problem categorizing Hx')
    end
    %Select environment:
    HitConst=nan(1,N);
    HitVariable=nan(1,N);
    if isnan(method(1))==1
        Env=rand(1,N)<p;
    elseif isnan(method(1))==0
        Env=method;
    end
    RConst=nan(1,N);
    RVariable=nan(1,N);
    
    mConst=m;
    mVariable=m;
    for t=1:N
        LSeq(t)=l;
        %Constant
        if sum(mConst)/l0<pI %Estimating RL
            RConst(t)=beta*rL*(1-Env(t));
            HitConst(t)=Env(t)==0; %1 if Env is 0
        else %Estimating RH
            RConst(t)=rL*(1-Env(t)) + rH*Env(t);
            HitConst(t)=Env(t)==1; %1 if Env is 1
        end
        %Update m, keeping l fixed.
        mConst=[mConst(2:end) Env(t)];
        
        %Variable
        if sum(mVariable)/l<pI %Estimating RL
            RVariable(t)=beta*rL*(1-Env(t));
            HitVariable(t)=Env(t)==0; %1 if Env is 0
        else %Estimating RH
            RVariable(t)=rL*(1-Env(t)) + rH*Env(t);
            HitVariable(t)=Env(t)==1; %1 if Env is 1
        end
        
       
        if strcmp(adaptation,'Proximity')==1 || strcmp(adaptation,'proximity')==1
            %Update l: Proximity function
            %LINEAR:
             lnPlus1=round( lMax-(lMax-lMin)*abs(sum(mVariable)/l-pI) );
             if lnPlus1>l %increase memory by 1
                if l==lMax %Already capped
                    %l unchanged
                    mVariable=[mVariable(2:end) Env(t)];
                else %l<lMax not capped
                    l=l+1;
                    mVariable=[mVariable(1:end) Env(t)];
                end
            elseif lnPlus1==l %Memory stays constant
                mVariable=[mVariable(2:end) Env(t)];
            elseif lnPlus1<l
                if l==lMin
                    mVariable=[mVariable(2:end) Env(t)];
                else %l > lMin
                     diff = 1;
                     l=l-1;
                    mVariable=[mVariable(2+diff:end) Env(t)];
                end
            else
                disp('problem categorizing ln vs. lnPlus1');
             end
            
        elseif strcmp(adaptation,'Variance')==1 || strcmp(adaptation,'variance')==1
           
            mVariablePlus1=[mVariable(2:end) Env(t)];
            Varm=sum(mVariable)/l*(1-sum(mVariable)/l);
            VarmPlus1=sum(mVariablePlus1)/l*(1-sum(mVariablePlus1)/l);
            
            pcutoff=1;
            if Varm < VarmPlus1
                Random=rand;
                if Random<pcutoff
                    if l==lMax
                        mVariable=[mVariable(2:end) Env(t)];
                    else %l<lMax not capped
                        l=l+1;
                        mVariable=[mVariable(1:end) Env(t)];
                    end
                else %Memory stays constant
                    mVariable=[mVariable(2:end) Env(t)];
                end
            elseif Varm > VarmPlus1 %New variance is smaller
                Random=rand;
                if Random<pcutoff
                    if l==lMin
                        mVariable=[mVariable(2:end) Env(t)];
                    else
                        l=l-1;
                        mVariable=[mVariable(3:end) Env(t)];
                    end
                else %Memory stays constant
                    mVariable=[mVariable(2:end) Env(t)];
                end
            elseif Varm==0 && VarmPlus1==0
                Random=rand;
                if Random<pcutoff
                    if l==lMin
                        mVariable=[mVariable(2:end) Env(t)];
                    else
                        l=l-1;
                        mVariable=[mVariable(3:end) Env(t)];
                    end
                else
                   mVariable=[mVariable(2:end) Env(t)];
                end
            else
                mVariable=[mVariable(2:end) Env(t)];
            end
            if l~=length(mVariable)
                disp('PROBLEM!')
                break
            end
        end
    end
        if length(Hx)==1
            RPrior=beta*rL*ones(1,l0);
            RConst=[RPrior RConst];
            RVariable=[RPrior RVariable];
            Env=[zeros(1,l0) Env];
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
            RConst=[RPrior RConst];
            RVariable=[RPrior RVariable];
        else
            disp('Problem categorizing Hx')
        end
end
