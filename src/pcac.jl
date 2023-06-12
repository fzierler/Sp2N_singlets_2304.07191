function awi_mass(file::AbstractString,typeUD::AbstractString;therm,step,autocor,ncut)
    AP = correlators(file,typeUD,"g0g5_g5_re")[:,therm+1:step:end]
    g5 = correlators(file,typeUD,"g5")[:,therm+1:step:end]
    return awi_mass(AP,g5;autocor,ncut)
end
function awi_mass(hdf5file,hdf5group,typeUD;therm,step,autocor,ncut)
    AP = correlators(hdf5file,hdf5group,typeUD,"g0g5_g5_re")[:,therm+1:step:end]
    g5 = correlators(hdf5file,hdf5group,typeUD,"g5")[:,therm+1:step:end]
    return awi_mass(AP,g5;autocor,ncut)
end
function awi_mass(AP,g5;autocor,ncut)
    T, N = size(AP)
    corr = zeros(T-3,N)
    for i in 1:N
        corr[:,i] = awi_correlator(AP[:,i],g5[:,i])
    end
    AWI, Δ2AWI = average_correlator(corr;autocor=autocor)[1:2]
    return awi_fit(AWI,sqrt.(Δ2AWI),ncut)
end
function awi_correlator(AP,g5)
    T = length(AP) # cuts for differentiation
    c1 = 2:T-2 # disregard the t=0 datapoint
    c2 = 3:T-1
    c3 = 4:T
    dtAP = @. (AP[c3]-AP[c1])/2 # perform differentiation (central difference)
    AWI = @. dtAP/g5[c2] / 2 # assemble final AWI-correlator
    return AWI
end
@. constfit(x,p) = p[1] + 0*x
function awi_fit(AWI,ΔAWI,ncut)
    # we have thrown away two datapoints at t=0,1 and one at t=T
    c = first(ncut)
    t0 = c-2:length(AWI)-c+1
    # initial guess at T/2
    p0 = AWI[div(length(AWI),2)]
    fit  = curve_fit(constfit,t0,AWI[t0],inv.(ΔAWI[t0]).^2,[p0])
    fitP = curve_fit(constfit,t0,AWI[t0].+ΔAWI[t0],inv.(ΔAWI[t0]).^2,[p0])
    fitM = curve_fit(constfit,t0,AWI[t0].-ΔAWI[t0],inv.(ΔAWI[t0]).^2,[p0])
    # extract mass and estimate error
    mAWI = fit.param[1]
    ΔmAWI1 = stderror(fit)[1]
    ΔmAWI2 = abs(fitP.param[1]-fitM.param[1])/2
    # fudge factor: not a completely symmetric plateau
    ΔmAWI = 2max(ΔmAWI1,ΔmAWI2)
    return mAWI, ΔmAWI
end
function nondeg_pcac_mass(mAWI,ΔmAWI,posdeg)
    mdeg, Δmdeg = mAWI[posdeg], ΔmAWI[posdeg]
    md  = @. 2mAWI - mdeg
    Δmd = @. sqrt(4ΔmAWI^2 + Δmdeg^2)
    return md, Δmd
end
