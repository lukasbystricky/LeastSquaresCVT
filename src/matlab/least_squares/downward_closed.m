function [index_set] = downward_closed(v)
% returns all indexes u <= v, such that u_i <= v_i for each i

D = size(v, 2);

index_set = [0:v(1)]';
for i=2:D
    c = [];
    for j=0:v(i)
       for this_row=index_set.'
           c = [c; this_row(:).' j];
       end       
    end
    index_set = c;
end

end