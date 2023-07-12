function effectivemass(c;step=1)
    N = length(c)
    m = zeros(N)
    for i in 1:N
        r = c[i]/c[mod1(i+step,N)]
        m[i] = log(abs(r))
    end
    return abs.(m)
end
function effectivemass_err(c,Δc)
    N  = length(c)
    Δm = zeros(N)
    for i in 1:N-1
        Δm[i] = sqrt((Δc[i]/c[i])^2 + (Δc[i+1]/c[i+1])^2)
    end
    Δm[N] = sqrt((Δc[N]/c[N])^2 + (Δc[1]/c[1])^2)
    return Δm
end
function effectivemass(c,Δc)
    return effectivemass(c), effectivemass_err(c,Δc)
end
function effectivemass_cosh(c)
    T = length(c)
    t = 1:T
    mid  = div(T,2) +1 # 1-based indexing
    return @. abs(acosh(c/c[mid])/(mid-t))
end
function effectivemass_cosh_err(c,Δc;norm=1)
    T = length(c)
    mid  = div(T,2) + 1 # 1-based indexing
    acoshderiv(x) = 1/sqrt(x^2 - 1)
    Δm = similar(c)
    for t in 1:T
        err1 = acoshderiv(c[t]/c[mid])*Δc[t]/c[mid]
        err2 = acoshderiv(c[t]/c[mid])*c[t]*Δc[mid]/c[mid]^2
        Δm[t] = LinearAlgebra.norm((err1,err2),norm)/abs(mid-t)
    end
    return Δm
end
function effectivemass_sinh(c)
    T = length(c)
    t = 1:T
    mid  = div(T,2) +1 # 1-based indexing
    return @. abs(asinh(c/c[mid])/(mid-t))
end
function effectivemass_sinh_err(c,Δc;norm=1)
    T = length(c)
    mid  = div(T,2) + 1 # 1-based indexing
    asingderiv(x) = 1/sqrt(x^2 + 1)
    Δm = similar(c)
    for t in 1:T
        err1 = asingderiv(c[t]/c[mid])*Δc[t]/c[mid]
        err2 = asingderiv(c[t]/c[mid])*c[t]*Δc[mid]/c[mid]^2
        Δm[t] = LinearAlgebra.norm((err1,err2),norm)/abs(mid-t)
    end
    return Δm
end
function effectivemass(file::String,type;key="g5",kws...)
    c,Δc2 = average_correlator(file,key,type;kws...)[1:2]
    m = effectivemass(c)
    Δm = effectivemass_err(c,sqrt.(Δc2))
    return m, Δm
end
function effectivemass_cosh(file::String,type;key="g5",kws...)
    c,Δc2 = average_correlator(file,key,type;kws...)[1:2]
    m = effectivemass_cosh(c)
    Δm = effectivemass_cosh_err(c,sqrt.(Δc2))
    return m, Δm
end
function effectivemass(file::String,typeU,typeD;key="g5",kws...)
    cU,Δc2U = average_correlator(file,key,typeU;kws...)[1:2]
    cD,Δc2D = average_correlator(file,key,typeD;kws...)[1:2]
    m  = effectivemass((cD+cU)/2)
    Δm = effectivemass_err((cD+cU)/2,sqrt.(( Δc2D+Δc2U)/2))
    return m, Δm
end
# implicit methods
# see equation (10) in arXiv:1607.06654 [hep-lat]
# (Notation in Gattringer/Lang is misleading!)
function _meff_at_t(c,t,T;sign=+1)
    # non-implicit mass as initial guess
    m0 = 1 # log(abs(c[t]/c[mod1(t+1,T)]))
    t0 = mod1(t-1,T)
    r0 = abs(c[t0]/c[t])
    # correlator at large times (dropped overall factor)
    cor_lt(m,T,t) = exp(-m*t) + sign*exp(-m*(t-T/2))
    # function to fit the effective mass
    g(m,T,t,t0) = cor_lt(m,T,t0)/cor_lt(m,T,t)
    # Use the more simpler algorithms from the Roots.jl package
    # find_zero() has more overhead and fails if the algorithm does not converged
    # Here we just use two simple, derivative free methods. If they do not converge
    # they return NaN. If that is the case then we try a slightly more robust algorithm.
    m = Roots.secant_method(x->g(x,T,t,t0)-r0,m0;maxevals=5000)
    if isnan(m)
       m = Roots.dfree(x->g(x,T,t,t0)-r0,m0)
    end
    return m
end
function implicit_meff(c;sign=+1)
    T = length(c)
    m = zeros(eltype(c),div(T,2))
    for t in 1:div(T,2)
        m[t] = _meff_at_t(c,t,T;sign=sign)
    end
    return m
end
function implicit_meff(CC,ΔCC;sign=+1)
    m1 = implicit_meff(CC;sign=sign)
    m2 = implicit_meff(CC + ΔCC;sign=sign)
    m3 = implicit_meff(CC - ΔCC;sign=sign)
    Δm = @. abs(m3-m2)/2
    return m1, Δm
end
function implicit_meff_jackknife(fileM,typeM,key;sign=+1,kws...)
    corrs = correlators(fileM,typeM,key;kws...)
    return implicit_meff_jackknife(corrs;sign=+1)
end
function implicit_meff_jackknife(corrs;sign=+1)
    T, N = size(corrs)
    # create arrays for decay constant
    corrs_delete1 = zeros(T,N-1)
    meff = zeros(N,T÷2)
    # set up jack-knife (one deletion)
    for i in 1:N
        for t in 1:T
            for j in 1:N
               (j < i) && (corrs_delete1[t,j]   = corrs[t,j])
               (j > i) && (corrs_delete1[t,j-1] = corrs[t,j])
            end
        end
        # perform averaging for fitting weights
        C = reshape(mean(corrs_delete1,dims=2),T)
        meff[i,:] = implicit_meff(C;sign)
    end
    return apply_jackknife(meff)
end
function effectivemass_cosh_jackknife(file::String,type,key;kws...)
    corrs = correlators(file,type,key)
    return implicit_meff_cosh(corrs)
end
function implicit_meff_cosh(corrs)
    T, N = size(corrs)
    # create arrays for decay constant
    corrs_delete1 = zeros(T,N-1)
    meff = zeros(N,T÷2)
    # set up jack-knife (one deletion)
    for i in 1:N
        for t in 1:T
            for j in 1:N
               (j < i) && (corrs_delete1[t,j]   = corrs[t,j])
               (j > i) && (corrs_delete1[t,j-1] = corrs[t,j])
            end
        end
        # perform averaging for fitting weights
        C = reshape(mean(corrs_delete1,dims=2),T)
        meff[i,:] = effectivemass_cosh(C)[1:T÷2]
    end
    return apply_jackknife(meff)
end
