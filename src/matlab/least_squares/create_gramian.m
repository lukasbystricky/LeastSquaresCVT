function [ G ] = create_gramian( gens, Phi )
%% create_gramian creates an approximate Gramian matrix
%
% Creates an approximation to the Gramian matrix needed to approximate a
% function over the hypercube [-1,1]^d. Given a set of n generators in R^d, 
% we % construct an approximate Gramian matrix G that can be used to create 
% the best approximation to a function in the space of polynomials of a 
% given order. The basis functions will be defined by a multi-index set 
% A_{ij}, such that \ell_i(x) = PROD( \phi_(A_{ij})(x_j) ). The multi-index 
% set has cardinality m, i.e. it defines m basis functions. The function 
% \phi_k(x_j) is a one-dimensional polynomial of degree k. 
%
% PARAMETERS:
% INPUTS:
% gens, REAL (n,d), the generators used to construct the matrix
% Phi, struct containing the following information:
%   - index set, INTEGER (m,d), the index set needed to construct the basis
%   - basis_card, INTEGER, the cardinality of the index set (i.e. m)
%   - value, FUNCTION HANDLE, a function handle that specifies how to
%   compute the basis, for example @(X) eval_phi_fast(X, Phi.index_set)
%
% OUTPUT:
% G, real(m,m), the Gramian matrix

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

G = zeros(m, m);

for i=1:m
    for j=1:m
        nui = Phi.index_set(i,:);
        nuj = Phi.index_set(j,:);
        for ii=1:n
            term = ws(ii);
            for jj=1:d
                term = term * legs(ii, nui(jj)+1, jj) * legs(ii,...
                    nuj(jj)+1, jj);
            end
            G(i,j) = G(i,j) + term;
        end
    end
end
end
