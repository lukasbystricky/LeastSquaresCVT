function [index_set] = index_set_rule(search_space,rule,s)
%Returns an index set based on given rule
%td- Total Degree space
%hc- Hyperbolic Cross
%s- order of rule
%%%For example for a given index set \nu in dimesions d the hc will be defined by
%%% \prod_{k=1}^{d} (\nu_k + 1) \leq s
%%%% s in the above formula is the order

[n,d] = size(search_space); %n is number of possible index sets we are searching over,d is dimension of stochastic space
index_set = [];

if(strcmp(rule,'hc'))
%if(rule == 'hc')
    for i=1:n
        val = 1;
        for j = 1:d
            val = val*(search_space(i,j)+1);
        end
        if(val <= s+1)
            index_set = [index_set,search_space(i,:)']; %#ok<AGROW>
        end
    end
end

if(strcmp(rule,'td'))
%if(rule == 'td')
    for i=1:n
        val = 0;
        for j = 1:d
            val = val + search_space(i,j);
        end
        if(val <= s)
            index_set = [index_set,search_space(i,:)']; %#ok<AGROW>
        end
    end
end

if(strcmp(rule,'td_aniso'))
%if(rule == 'td')
    for i=1:n
        val = sum(search_space(i,:) .* (1:2:2*d));
        if(val <= s)
            index_set = [index_set,search_space(i,:)']; %#ok<AGROW>
        end
    end
end

if(strcmp(rule,'tensor'))
%if(rule == 'tensor')
    for i=1:n
        if all(search_space(i,:) <= s)
            index_set = [index_set,search_space(i,:)']; %#ok<AGROW>
        end
    end
end

index_set = index_set';

end

