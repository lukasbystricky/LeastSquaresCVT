function [sample_points] = sequential_sampling_uniform(n,Phi)
%n - number of desired sample points
%We are using the uninform density with a Legendre Polynomial basis
[basis_card,d] = size(Phi.index_set); %Size of index set and dimension of final random variable
%sample_points = zeros(n,d); % Array to store final random points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start with q = 1 case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Linfty = sqrt(2*max(Phi.index_set(:,1)) + 1); %Maximum value of the polynmial space i.e. for f(x) < Mg(x) in rejection sampling this value is M
a = -1;
b = 1;
tot_samples = 2*n*max(Phi.index_set(:,1))^2; %We will need to consider 2*d*n*basis_card^2 sample points see mig/cohen paper

sample_points = zeros(n,1);
for i = 1:n
    flag = 1;
    while(flag == 1)
        r = a + (b-a).*rand(tot_samples,1);
        rej_samples = (1/basis_card)/(Linfty)*sum(eval_phi_fast(r,Phi.index_set(:,1)).^2,2);% this is f(x)/(Mg(x) somethings simplified ask for detials if not clear
        u = rand(tot_samples,1);
        accept = u < rej_samples;
        fin_accept = find(accept,1); %Find first points which was accepted
        if(length(fin_accept) == 1)
            sample_points(i) = r(fin_accept);
            flag = 0;
        end
    end
end




    for q = 2:d
        %%%%Calculate all the marginals in one go%%%%%%
        marginals = bsxfun(@rdivide,eval_phi_fast(sample_points,Phi.index_set(:,1:q-1)).^2,sum(eval_phi_fast(sample_points,Phi.index_set(:,1:q-1)).^2,2));
        check = sum(marginals,2); %sanity check this should be equal to 1 by definition
        %%%Begin Rejection Sampling%%%%%%
        %%%%I don't think we can avoid using a for loop over the total number of desired sample points here
        %%%% Each of the initial points has its own associated marginal so will
        %%%% need its individual rejection sample
        sample_points_temp = zeros(n,1);
        for i = 1:n
            flag = 1;
            while(flag == 1)
                tot_samples = 2*max(Phi.index_set(:,q))^2;
                r = a + (b-a).*rand(tot_samples,1); %random samples on [-1,1] for rejection sampling
                Linfty = sqrt(2*max(Phi.index_set(:,q)) + 1);
                rej_samples = 1/(Linfty)*eval_phi_fast(r,Phi.index_set(:,q)).^2 * marginals(i,:)';% this is f(x)/(Mg(x) somethings simplified ask for detials if not clear
                u = rand(tot_samples,1);
                accept = u < rej_samples;
                fin_accept = find(accept,1); %Find first points which was accepted
                if(length(fin_accept) == 1)
                    sample_points_temp(i) = r(fin_accept);
                    flag = 0;
                end
            end
        end
        sample_points = [sample_points sample_points_temp];

        %marginal(:) = prod(eval_phi_fast(sample_points(i,q-1)
        %for j = 1:basis_card
        %    marginals(j) = prod(eval_phi_fast(sample_points(i,q-1),Phi.index_set(:,1)).^2);
        %end
        %r = a + (b-a).*rand(tot_samples,1); %random samples on [-1,1] for rejection sampling
        %rej_samples = (1/basis_card)/(Linfty)*sum(eval_phi_fast(r,Phi.index_set(:,q).^2),2);% this is f(x)/(Mg(x) somethings simplified ask for detials if not clear
%%%%Now we need to determine the conditional probability%%%%%%
%     for j = 1:q
%         eval_arr = eval_leg(phi.index_set(:,j), sample_points(:,j));
   
    end






