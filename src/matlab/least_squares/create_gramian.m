function [ A ] = create_gramian( gens, Phi )

[n,~] = size(gens);
[m,d] = size(Phi.index_set);
s = max(Phi.index_set(:));

legs = zeros(n, s+1, d);
Pweight = @(x) sum(Phi.value(x).^2, 2);
ws = zeros(n, 1);

for i=1:d
    legs(:,:,i) = eval_leg(s, gens(:,i));
end

for i=1:n
    ws(i) = m/Pweight(gens(i,:))/n;
end

A = zeros(m, m);

for i=1:m
    for j=1:m
        nui = Phi.index_set(i,:);
        nuj = Phi.index_set(j,:);
        for ii=1:n
            term = ws(ii);
            for jj=1:d
                term = term * legs(ii, nui(jj)+1, jj) * legs(ii, nuj(jj)+1, jj);
            end
            A(i,j) = A(i,j) + term;
        end
    end
end
end
