function [sample_points] = gen_sequential_sample_points(d, rule, order, num_sample)
rng('shuffle')
%Driver file for least squares optimization
%d - dimension of hypercube
%rule - sparse tensor basis rule
%order of the polynomial space we want
switch rule
    case {'tensor'}
        poly_order = d * order;
    case ('td')
        poly_order = order;
        full_index_set = comb_hold(poly_order,d); %Returns full tensor product index set up to order search_val 
        index_set = index_set_rule(full_index_set,rule,order); %Returns index set based on given rule and possible weights
        %account for the fact that the index set returned did not include the constant function
        index_set = [zeros(1, d); index_set];
    % 'dc' means downward closed, would return an arbitrary index with 'order'
    % cardinality, and dimension d. 'arbitrary_index' automatically adds zeros
    case ('dc')
        index_set = arbitrary_index(order,d); %what would be called normally
        
end

%M_full = i4_choose(poly_order+d, d) - 1;
%commented out so that the "order" argument can be used with the arbitrary
%   cardinality without it thinking the index set will be too big
%if d * M_full > 10 * 1000000
%    error('attempt to construct d=%i dimensional index set with n=%i entries.', d, M_full)
%end



basis_card = size(index_set, 1);

% commented out for the same reason as above
%if d * basis_card > 10 * 1000
%    error('attempt to construct d=%i dimensional quadrature with n=%i polynomials.', d, M_full)
%end

% encode underlying domain (here unit cube)
Omega = struct();
Omega.type = 'cube';
Omega.d = d;
Omega.box = repmat([-1;1], 1, d); % domain Omega =[-1,1]^d
Omega.vol = 2^d;  % volume of domain

% setup structure for sensitivities (put gradient etc. also here later)
% extensible for different polynomials
% hardcode index set for now (TODO think about efficency later, i.e., to reduce the number of evaluations)
Phi = struct();
Phi.index_set = index_set;
Phi.basis_card = basis_card;
Phi.value = @(X) eval_phi_fast(X, index_set);


% generate some Migliorati samples
%N_sample = 3 * basis_card; %TODO implement exact formula that guarantee cond(G) <= 3
%N_sample = 4 * basis_card * log(basis_card); %TODO implement exact formula that guarantee cond(G) <= 3
Pweight = @(x) sum(Phi.value(x).^2, 2) / basis_card;
gens = sequential_sampling_uniform(num_sample,Phi);
cvt_points = cvt_generator(gens, Phi);

disp(cond(create_gramian(gens, Phi)));

scatter(cvt_points(:,1), cvt_points(:,2))

end