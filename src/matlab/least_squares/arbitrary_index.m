function [index_set] = arbitrary_index(M,D)
% Construct a randomly generated index set of cardinality M, dimension D

index_set = zeros(1, D);

for i=1:D  %ensures that every dimension has at least one index
    this_index = zeros(1,D);
    this_index(i) = 1;
    index_set = [index_set ; this_index];
    
end

S = 1;
while i4_choose(S + D, D) < M
    S = S + 1;
end

test_set = comb_hold(S,D);

index_count = 1 + D;

while index_count < M
    vs = randperm(size(test_set,1));
    
    for i=vs
        bad = 0;
        this_index = test_set(i, :);
        downwards = downward_closed( this_index );
        temp_set = [index_set ; this_index];
        for v_bar=downwards.'
            if ~is_in_index(temp_set, v_bar)
                bad = 1;
                break
            end
        end
        
        if ~bad && ~is_in_index( index_set, this_index.' )
            index_set = [index_set ; this_index];
            test_set(i, :) = [];
            index_count = index_count + 1;
            break
        end
    end
end

end