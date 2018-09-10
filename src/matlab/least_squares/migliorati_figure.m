function [] = migliorati_figure( dim, trials )
nRange = 100;
mRange = 30;

cond_matrix_mig = zeros(mRange,nRange);
cond_matrix_cvt = zeros(mRange,nRange);

for n=1:nRange
    bad = 0;
    for m=1:mRange
        disp([n+2 m+2])
        for k=1:trials
            Phi.index_set = arbitrary_index(m+2,dim);
            Phi.basis_card = m+2;
            Phi.value = @(X) eval_phi_fast(X, Phi.index_set);
            
            gens = sequential_sampling_uniform(n+2,Phi);
            
            [fgens, fCond, cond_flag] = cvt_generator(gens, Phi);
            
            if fCond < 3
                cond_matrix_mig(m, n) = cond_matrix_mig(m, n) + 1; 
                if cond_flag
                    cond_matrix_cvt(m, n) = cond_matrix_cvt(m, n) + 1;     
                end
            end
        end
        if cond_matrix_cvt(m, n) == 0
            bad = bad + 1;
            if bad >= 3
                break
            end
        end
    end
end

figure(1)
[~,h1] = contourf(cond_matrix_mig,trials);
set(h1,'LineColor','none')

figure(2)
[~,h2] = contourf(cond_matrix_cvt,trials);
set(h2,'LineColor','none')

save('ConditionMatrixTrialMig.mat', 'cond_matrix_mig');
save('ConditionMatrixTrialCVT.mat', 'cond_matrix_cvt');
end