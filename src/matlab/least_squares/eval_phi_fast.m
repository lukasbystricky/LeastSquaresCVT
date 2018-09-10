function [phi] = eval_phi_fast(x, index_set)
%x all points being passed in

[npts, dx] = size(x);
[index_size, d] = size(index_set);

assert(dx == d);

max_PD = max(index_set); %Maximum polynomial degree for each dimension

phi = ones(npts, index_size);
for j = 1:d
     eval_arr = eval_leg(max_PD(j), x(:,j));
     phi_j = eval_arr(:, index_set(:,j)+1);
     phi = phi .* phi_j;
end

% for i = 2:index_size
%      index_hold = repmat(index_set(i,:),npts,1);
%      phi(:,i) = prod(legendreP(index_hold,x),2);
% end

end