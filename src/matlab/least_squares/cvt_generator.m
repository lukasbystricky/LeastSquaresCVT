function [gens, fCond, cond_flag] = cvt_generator(gens, Phi)
nBigIter = 0;

cond_flag = 0;

fCond = cond(create_gramian(gens, Phi));
%disp(fCond)
if fCond < 3
    cond_flag = 1;
    return;
end

prevWeights = zeros(size(gens, 1));
[gens, prevWeights] = kmeans(gens, Phi, 5000, prevWeights);
fCond = cond(create_gramian(gens, Phi));
disp(fCond)
if fCond < 3
    return;
end
    
for i=1:nBigIter
    [gens, prevWeights] = kmeans(gens, Phi, 100, prevWeights);
    fCond = cond(create_gramian(gens, Phi));
    %disp(fCond)
    if fCond < 3
        return;
    end
end

end