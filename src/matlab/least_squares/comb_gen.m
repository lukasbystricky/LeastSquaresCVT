function [A] = comb_gen(N,D)
%Call comb_function to return all possible combinations of N into K parts
%N is the number
%K is the number of parts we wish to split N into

p = i4_choose(N+D,D) - i4_choose(N+D-1,D); %determine the amount of memory needed
A = zeros(p,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Set Default Values%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MTC = false; 
T = 1000000;
H = 1000000;
R = zeros(1,D);
for i=1:p
    %i
    [N,D,H,R,T,MTC] = comb_function(N,D,H,R,T,MTC);
    A(i,:) = R(1,:);
end
%A

