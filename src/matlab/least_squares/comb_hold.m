function [A] = comb_hold(N,D)
%Holds all the possible combinations of a number from 1 to N split into D
%Parts
%i.e. N = 2 D = 3 output will be 
%
%1 0 0
%0 1 0
%0 0 1
%2 0 0
%1 1 0
%0 2 0
%1 0 1
%0 1 1
%0 0 2
%D is the dimension
%N is the degree the Vandermonde matrix will go up to

A = zeros(i4_choose(N+D,D)-1,D);
count = 1;

for i=1:N
    p = i4_choose(i+D,D) - i4_choose(i+D-1,D);
    R = zeros(p,D);
    R = comb_gen(i,D);
    R;
    
    for k =1:p
        for j =1:D
            A(count,j) = R(k,j);
        end
        count = count + 1;
    end
end


end

