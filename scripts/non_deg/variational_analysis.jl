function eigenvalues_jackknife_samples(corr)
    sample = delete1_resample(corr)
    nops, T, nconf = size(sample)[2:4]
    eigvals_jk = zeros(eltype(sample),(nops,T,nconf))
    for t in 1:T, s in 1:nconf
        t0 = 1
        # smaller values correspond to a faster decay, and thus correspond to a larger masses
        # use sortby to sort the eigenvalues by ascending eigen-energy of the meson state
        eigvals_jk[:,t,s] = eigen(sample[:,:,t,s],sample[:,:,t0,s],sortby=x -> -abs(x)).values
    end
    return eigvals_jk
end
function eigenvalues_jackknife(corr)
    eigvals_jk = eigenvalues_jackknife_samples(corr)
    eigvals, Δeigvals = apply_jackknife(eigvals_jk;dims=3)
    return eigvals, Δeigvals
end
# generate a resample of the original correlator matrix
function delete1_resample(corr_matrix)
    nops,T,nconf = size(corr_matrix)[2:end]
    samples = similar(corr_matrix)
    # temporary array for jackknife sampling
    temp = zeros(eltype(corr_matrix),(nops,nops,T,nconf-1))
    for index in 1:nconf    
        for i in 1:index-1
            temp[:,:,:,i] = corr_matrix[:,:,:,i]
        end
        for i in 1+index:nconf
            temp[:,:,:,i-1] = corr_matrix[:,:,:,i]
        end
        # perform average after deleting one sample
        samples[:,:,:,index] = dropdims(mean(temp,dims=4),dims=4)
    end
    return samples
end
# apply jackknife resampling along dimension dims
function apply_jackknife(obs::AbstractArray;dims::Integer)
    N  = size(obs)[dims]
    O  = dropdims(mean(obs;dims);dims)
    ΔO = dropdims(sqrt(N-1)*std(obs;dims,corrected=false);dims)
    return O, ΔO
end
function apply_jackknife(obs::AbstractVector)
    N  = length(obs)
    O  = mean(obs)
    ΔO = sqrt(N-1)*std(obs,corrected=false)
    return O, ΔO
end
