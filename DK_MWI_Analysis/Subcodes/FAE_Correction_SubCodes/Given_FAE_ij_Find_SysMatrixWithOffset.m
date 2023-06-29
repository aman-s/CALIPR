
function H = Given_FAE_ij_Find_SysMatrixWithOffset(DeltaVox, alphaEx_NominalRad, alphaRef_NominalRad, TENew, esp, T1, T2_Scale, OffsetTag)
    EchoSampled =  int16(TENew./esp);
    [alphaEX, alphaRE] = Calculate_Exc_n_Ref_FA(alphaEx_NominalRad, alphaRef_NominalRad, DeltaVox);
    
    R1 = 1./T1;
    R2vec = 1./T2_Scale;
    tau = esp;
    n = max(EchoSampled);
    %Converting n to double for sparse matrix used later
    n = double(n);
    
    nRates = length(T2_Scale);
    tau=tau/2;

    H0 = zeros(n,nRates);
        
    % RF mixing matrix
    T0 = [cos(alphaRE/2)^2, sin(alphaRE/2)^2, sin(alphaRE); ...
        sin(alphaRE/2)^2, cos(alphaRE/2)^2, -sin(alphaRE); ...
        -0.5*sin(alphaRE), 0.5*sin(alphaRE), cos(alphaRE)];

    TArray = cell(1,n);
    TArray(:) = {sparse(T0)};
    T = blkdiag(1,TArray{:});

    % Selection matrix to move all traverse states up one coherence level
    S = sparse(3*n+1,3*n+1);
    S(1,3)=1; S(2,1)=1;
    S(3,6)=1; S(4,4)=1;
    for o=2:n
        offset1=( (o-1) - 1)*3 + 2;
        offset2=( (o+1) - 1)*3 + 3;
        %display(['offset1: ', num2str(offset1), '; ','offset2: ', num2str(offset2)]);
        if offset1<=(3*n+1)
            S(3*o-1,offset1)=1;  % F_k <- F_{k-1}
        end
        if offset2<=(3*n+1)
            S(3*o,offset2)=1;  % F_-k <- F_{-k-1}
        end
        S(3*o+1,3*o+1)=1;  % Z_order
    end

    for iRate=1:nRates
        % Relaxation matrix
        R2=R2vec(iRate);
        R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);

        RArray = cell(1,n);
        RArray(:) = {sparse(R0)};
        R = blkdiag(exp(-tau*R2),RArray{:});

        % Precession and relaxation matrix
        P = (R*S);

        % Matrix representing the inter-echo duration  
        E = (P*T*P);

        % Recursively apply matrix to get echo amplitudes
        x = zeros(size(R,1),1);
        x(1)=1*sin(alphaEX);
        %x(4)=1*cos(alphaEX);
      
        %Effect of first echo:
        for iEcho=1:n
            x=E*x;  
            H0(iEcho,iRate) = x(1);
        end
    end
    
    %Return only for echo sampled
    H = H0(EchoSampled, :);
    if (strcmp(OffsetTag, 'Yes') || strcmp(OffsetTag, 'yes'))
        H = [H ones(size(H,1), 1)]; % Add one more row for taking into account the baseline offset.
    end
end

function [alpha_Excite, alpha_Refocus] = Calculate_Exc_n_Ref_FA(alphaEx_NominalRad, alphaRef_NominalRad, Delta_FAE)
    alpha_Excite    = (1 + Delta_FAE)*alphaEx_NominalRad;
    alpha_Refocus   = (1 + Delta_FAE)*alphaRef_NominalRad;
end