using Plots
using HiRepAnalysis
using LinearAlgebra
using LaTeXStrings
using DelimitedFiles
pgfplotsx(legendfontsize=20,labelfontsize=24,tickfontsize=18,titlefontsize=22,ms=4,framestyle=:box,legend=:outerright,size=(800,500))

function errorbars_semilog(x,Δx)
    lower = similar(Δx)
    for i in eachindex(x)
        lower[i] = ifelse(x[i] > Δx[i],Δx[i],x[i]-1E-15)
    end
    bars = (lower,Δx)
    return bars
end
parse_fit_intervall(s) = isa(s,Integer) ? s : Tuple(parse.(Int,split(s[2:end-1],',')))
function read_prm(file,i)
    prm = readdlm(file,';')
    odir = prm[i,1]
    hits = prm[i,2]
    fitη,fitπ,fitσ,fita0,fitρ = parse_fit_intervall.(prm[i,3:7])
    fileM, fileS = joinpath(odir,"out_spectrum.hdf5"), joinpath(odir,prm[i,9])
    group,name,vsub,deriv,gs_sub,constant = prm[i,10:15]
    return odir, hits, fitη, fitπ, fitσ, fita0, fitρ ,fileM, fileS, group, name, vsub, deriv, gs_sub, constant
end
function singlet_analysis(hdf5file,hdf5groupM,hdf5groupS,hits,vsub,key,cut,gs_sub,sigma)
    type  = "DISCON_SEMWALL SINGLET"
    typeM = "DEFAULT_SEMWALL TRIPLET"
    T, L = latticesize(hdf5file,hdf5groupM)[1:2]
    rescale = (L^3)^2 /L^3
    # disconnected contributions from code in /Spectrum/
    C_dis_MC = disconnected_eta_MC(hdf5file,hdf5groupS,type,hits;maxhits=hits,rescale,vsub,key="$(key)_disc_re")
    C_con_MC = correlators(hdf5file,hdf5groupM,typeM,key)
    HiRepAnalysis._rescale_corrs!(C_con_MC, L)
    return HiRepAnalysis.singlet_jackknife_corr(C_con_MC,C_dis_MC,cut;gs_sub,sigma)
end
function constant_plot()
    prmfile = "input/parameters/param_deriv.csv"

    # Compare eta vs pion correlator in presence of a constant
    plt_corr_pi_eta = plot()
    odir,hits,fitη,fitπ,fitσ,fita0,fitρ,fileM,fileS,group,name,vsub,deriv,gs_sub,constant = read_prm(prmfile,11)
    hdf5groupM = splitext(fileM)[1]
    hdf5groupS = splitext(fileS)[1]

    Nconf = nconfigs(hdf5file,hdf5groupM)
    m = quarkmasses(hdf5file,hdf5groupM)
    G = gaugegroup(hdf5file,hdf5groupM)
    β = couplingβ(hdf5file,hdf5groupM)
    T, L = latticesize(hdf5file,hdf5groupM)[1:2]
    typeM, type = "DEFAULT_SEMWALL TRIPLET", "DISCON_SEMWALL SINGLET"
    rescale = (L^3)^2 /L^3

    # connected pion correlator
    Cπ, ΔCπ2 = average_correlator(hdf5file,hdf5groupM,"g5",typeM)[1:2]
    HiRepAnalysis._rescale_corrs!(Cπ, ΔCπ2, L)
    ΔCπ = sqrt.(ΔCπ2)
    # connected eta' correlator
    gs_sub = false # do not perform ground state subtraction
    C, ΔC, Cd, ΔCd = singlet_analysis(hdf5file,hdf5groupM,hdf5groupS,hits,vsub,"g5",fitπ,gs_sub,false)
    #title for plotting
    title = L"$ %$(T)\times %$(L)^3, \beta=%$β, %$G, m_q=%$(first(m))$, $n_{\rm src} = %$hits$"
    label = L"C_{\eta'}(t)"
    plot!(plt_corr_pi_eta,title=title,legend=:top,yaxis=:log10)
    scatter!(plt_corr_pi_eta,2Cπ,yerr=errorbars_semilog(2Cπ,2ΔCπ),label=L"2C_{\pi}(t)")
    scatter!(plt_corr_pi_eta,C,yerr=errorbars_semilog(C,ΔC),label=label,markershape=:diamond,ms=5)
    plot!(plt_corr_pi_eta,yticks=10.0.^(-4:0),xticks=2:2:T,xlabel=L"t")

    # Compare eta' correlator for different fixed Q
    plt_corr_eta_fQ = plot()
    odir,hits,fitη,fitπ,fitσ,fita0,fitρ,fileM,fileS,group,name,vsub,deriv,gs_sub,constant = read_prm(prmfile,23)
    hdf5groupM = splitext(fileM)[1]
    hdf5groupS = splitext(fileS)[1]

    absQ = true
    markers  = deleteat!(Plots.supported_markers(),(1,2,4,7,9,13,14,16,17,18))
    for Q in [1,2,3,4]
        hdf5fileQ    = absQ ? "output/data_fixed_absQ.hdf5" : "output/data_fixed_Q.hdf5"
        hdf5groupM_Q = joinpath(hdf5groupM,"Q$(Float64(Q))")
        hdf5groupS_Q = joinpath(hdf5groupS,"Q$(Float64(Q))")

        Nconf = nconfigs(hdf5fileQ,hdf5groupM_Q)
        m = quarkmasses(hdf5fileQ,hdf5groupM_Q)
        G = gaugegroup(hdf5fileQ,hdf5groupM_Q)
        β = couplingβ(hdf5fileQ,hdf5groupM_Q)
        T, L = latticesize(hdf5fileQ,hdf5groupM_Q)[1:2]

        # parameter of configurations
        typeM, type = "DEFAULT_SEMWALL TRIPLET", "DISCON_SEMWALL SINGLET"
        rescale = (L^3)^2 /L^3

        # connected eta' correlator
        gs_sub = false # do not perform ground state subtraction
        C, ΔC, Cd, ΔCd = singlet_analysis(hdf5fileQ,hdf5groupM_Q,hdf5groupS_Q,hits,vsub,"g5",fitπ,gs_sub,false)
        #title for plotting
        title = L"$ %$(T)\times %$(L)^3, \beta=%$β, %$G, m_q=%$(first(m))$, $n_{\rm src} = %$hits$"
        if absQ
            label = L"C_{\eta'}(t), |Q|=%$(Q), N_{\rm cfg} = %$Nconf"
        else
            label = L"C_{\eta'}(t), Q=%$(Q), N_{\rm cfg} = %$Nconf"
        end
        plot!(plt_corr_eta_fQ,title=title,legend=:top,yaxis=:log10)
        scatter!(plt_corr_eta_fQ,C,yerr=errorbars_semilog(C,ΔC),label=label,markershape=markers[Q],ms=5)
        plot!(ylims=(max(1E-5,minimum(C)/10),ylims(plt_corr_eta_fQ)[2]))
        plot!(plt_corr_eta_fQ,yticks=10.0.^(-4:0),xticks=2:2:T,xlabel=L"t")
    end
    savefig(plt_corr_eta_fQ,"output/figures/fixedQ_correlator.pdf")
    savefig(plt_corr_pi_eta,"output/figures/Constant_in_eta_shifted.pdf")
end
constant_plot()
