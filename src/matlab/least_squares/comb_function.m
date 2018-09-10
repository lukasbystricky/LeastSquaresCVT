function [N,K,H,R,T,MTC] = comb_function(N,K,H,R,T,MTC )
%This function is based on the algorithim nexcom from the book
% Combinitoral Algorithims for computers and calculators Second Edition
%It generates the next compisition of an intger N into K parts
%N the integer
%K the number of parts we wish to split N
%MTC True if this is not the last composition, otherwise false
%R the Ith part of the output composition
%H,T placeholders

if (MTC == true)
    if (T>1)
        H = 0;
    end
    
    H = H+1;
    T=R(H);
    R(H) = 0;
    R(1) = T-1;
    R(H+1) = R(H+1) + 1;
    
    if (R(K) ==  N)
        MTC = false;
        return
    else
        MTC = true;
        return
    end
end

if(MTC == false)
    R(1) = N;
    T = N;
    H = 0;
    
    for i = 2:K
        R(i) = 0;
    end
    
    if (R(K) ==  N)
        MTC = false;
        return
    else
        MTC = true;
        return
    end
    
end

end

