function error_distribution(corr,cov,ncut;nexp2,sign=+1,constant=false)
    n     = 10_000
    cov   = (cov+cov')/2
    # If the matrix is not positive definite then there must be a numerical
    # problem. For now just ignore the off-diagonal elements, i.e. the
    # correlation
    if !isposdef(cov)
        @show size(cov)
        cov = LinearAlgebra.diagm([cov[i,i] for i in 1:size(cov)[1]])
    end
    mass  = zeros(n)
    decay = zeros(n)
    w = inv(cov)
    dist = MvNormal(corr,cov)
    m = 0
    for i in 1:n
        c = rand(dist)
        if constant
            f = fit_corr_LSD(c,cov,ncut;sign)[1]
        else
            f = fit_corr(c,cov,ncut;nexp2,sign)[1]
        end
        if !f.converged
            i = i -1
            m = m +1
            continue
        end
        mass[i] = f.param[1]
        decay[i] = f.param[2]
    end
    return mass, decay
end

function data_binning(s,nbins=200)
    counts = zeros(nbins)
    bins = collect(range(minimum(s),maximum(s),length=nbins))
    index = zeros(Int,length(s))
    for i in eachindex(s)
        if s[i] < bins[end]
            index = findfirst(x-> x >= s[i],bins)
            counts[index] += 1
        end
    end
    # normalize with trapezoid rule
    h = bins[2]-bins[1]
    norm = h*(sum(counts) - counts[1]/2 - counts[end]/2 )
    return counts/norm, bins
end
function fit_hist2normal(bins,hist,median)
    @. fitnormal(x,p) = pdf(Normal(p[1],abs(p[2])),x)
    p = [median,median/10]
    f = curve_fit(fitnormal,bins,hist,p)
    μ,σ = f.param[1],abs(f.param[2])
    return μ,σ
end
function decay_mass_histogram(corr,cov_m,ncut;nexp2,sign=+1)
    mass, decay = error_distribution(corr,cov_m,ncut;nexp2,sign)
    median_mass = median(mass)
    median_decay = median(decay)
    hist_mass,  bins_mass  = data_binning(filter(x -> x < 5median_mass,  mass ))
    hist_decay, bins_decay = data_binning(filter(x -> x < 5median_decay, decay))
    m,Δm = fit_hist2normal(bins_mass, hist_mass , median_mass )
    f,Δf = fit_hist2normal(bins_decay,hist_decay, median_decay)
    return m, Δm, f, Δf
end
#function decay_mass_histogram(file::AbstractString,key,typeUD;sign=+1,nexp2,ncut,kws...)
#    corr,cov_v,cov_m = average_correlator(file,key,typeUD;kws...)
#    L = latticesize(file)[2]
#    _rescale_corrs!(corr,cov_v,cov_m,L)
#    m, Δm, f, Δf = decay_mass_histogram(corr,cov_m,ncut;nexp2,sign)
#    return m, Δm, f, Δf
#end
#function decay_mass_histogram(file::AbstractString,key,typeU,typeD;sign=+1,ncut,nexp2,kws...)
#    corr,cov_v,cov_m = average_correlator(file,key,typeU,typeD;kws...)
#    L = latticesize(file)[2]
#    _rescale_corrs!(corr,cov_v,cov_m,L)
#    m, Δm, f, Δf = decay_mass_histogram(corr,cov_m,ncut;nexp2,sign)
#    return m, Δm, f, Δf
#end
function decay_mass_histogram(hdf5file,hdf5group,key,typeUD;sign=+1,nexp2,ncut,kws...)
    corr,cov_v,cov_m = average_correlator(hdf5file,hdf5group,key,typeUD;kws...)
    L = latticesize(hdf5file,hdf5group)[2]
    _rescale_corrs!(corr,cov_v,cov_m,L)
    m, Δm, f, Δf = decay_mass_histogram(corr,cov_m,ncut;nexp2,sign)
    return m, Δm, f, Δf
end
function decay_mass_histogram(hdf5file,hdf5group,key,typeU,typeD;sign=+1,ncut,nexp2,kws...)
    corr,cov_v,cov_m = average_correlator(hdf5file,hdf5group,key,typeU,typeD;kws...)
    L = latticesize(hdf5file,hdf5group)[2]
    _rescale_corrs!(corr,cov_v,cov_m,L)
    m, Δm, f, Δf = decay_mass_histogram(corr,cov_m,ncut;nexp2,sign)
    return m, Δm, f, Δf
end