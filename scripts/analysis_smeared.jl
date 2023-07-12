using HiRepAnalysis
using Plots
using Statistics
using LinearAlgebra
using LaTeXStrings
pgfplotsx(legendfontsize=18,labelfontsize=20,tickfontsize=18,titlefontsize=20,ms=5,framestyle=:box,legend=:outerright)

hdf5groupN  = "runsSp4/Lt32Ls16beta6.9m1-0.90m2-0.90/out_spectrum/"
hdf5groupD  = "runsSp4/Lt32Ls16beta6.9m1-0.90m2-0.90/out_spectrum_discon/"
hdf5groupS1 = "runsSp4/Lt32Ls16beta6.9m1-0.90m2-0.90/out_spectrum_smeared_N60N0/"
hdf5groupS2 = "runsSp4/Lt32Ls16beta6.9m1-0.90m2-0.90/out_spectrum_smeared_N60N60/"
typeN = "DEFAULT_SEMWALL TRIPLET"

T = latticesize(hdf5file,hdf5groupN)[1]
L = latticesize(hdf5file,hdf5groupN)[2]
key  = "g5"

# We have fewer measurements for with smearing than measurements in total
# We restrict ourselves only to those ensembles where we have smearinf
n = nconfigs(hdf5file,hdf5groupS1)

# obtain connected part without smearing
kws = (therm=1,step=1,autocor=false)
corr  = correlators(hdf5file,hdf5groupN,typeN,key)[:,1:n]
C  = vec(mean(corr,dims=2))
ΔC = vec(std(corr,dims=2))./sqrt(n)

# obtain connected part
typeS = "source_N60_sink_N0 TRIPLET"
kws = (therm=1,step=1,autocor=false)
corr_smear0 = correlators(hdf5file,hdf5groupS1,typeS,key)

typeS = "source_N60_sink_N60 TRIPLET"
kws = (therm=1,step=1,autocor=false)
corr_smearN = correlators(hdf5file,hdf5groupS2,typeS,key)

ratio = similar(corr_smear0)
for i in eachindex(ratio)
    ratio[i] =  corr_smear0[i] .^2 ./ corr_smearN[i]
end
ratio
r  = vec(mean(ratio,dims=2))
Δr = vec(std(ratio,dims=2))/sqrt(n)


# rescale with same factors as in _rescale_corrs
LF =  L^3/2
@. r  *= LF
@. C  *= LF
@. Δr *= LF
@. ΔC *= LF

# disconnected contributions from code in /Spectrum/
typeD = "DISCON_SEMWALL SINGLET"
hits  = 128
rescale = (L^3)^2 /L^3
D, ΔD = disconnected_eta(hdf5file,hdf5groupD,typeD,hits;rescale=rescale,maxhits=hits,maxconf=n)

# obtain improved connected correlator
kws = (ncut=8,error=:hist,nexp2=false)
CI, ΔCI = groundstate_correlator(hdf5file,hdf5groupN,key,typeN;kws...)
mπ, Δmπ, fπ, Δfπ = HiRepAnalysis.meson_mass_decay(hdf5file,hdf5groupN,key,typeN;kws...)[1:4]
mρ, Δmρ = HiRepAnalysis.meson_mass_decay(hdf5file,hdf5groupN,"g1",typeN;kws...)[1:2]
sπ  = errorstring(mπ,Δmπ)
sρ  = errorstring(mρ,Δmρ)
sfπ = errorstring(fπ,Δfπ)

# Improved singlet correlator
Im  = CI - 2D
ΔIm = @. sqrt(ΔCI^2 + 4ΔD^2)

# Smeared singlet correlator
CS  = r - 2D
ΔCS = @. sqrt(Δr^2 + 4ΔD^2)

# direct singlet correlator
CD  =  C -  2D
ΔCD = @. sqrt(ΔC^2 + 4ΔD^2)

# effectivemasses
meffC,  ΔmeffC  = implicit_meff(C,  ΔC)
meffIm, ΔmeffIm = implicit_meff(Im, ΔIm)
meffCS, ΔmeffCS = implicit_meff(CS, ΔCS)
meffCD, ΔmeffCD = implicit_meff(CD, ΔCD)

# extract effective mass from smeared
fitint = (5,8)
mη, Δmη, fη, Δfη = HiRepAnalysis.decay_mass_histogram(CS,diagm(ΔCS.^2),fitint;nexp2=false)[1:4]
mη, Δmη, fη, Δfη = HiRepAnalysis.decay_mass_histogram(Im,diagm(ΔIm.^2),fitint;nexp2=false)[1:4]
sη = errorstring(mη,Δmη)
fitη = first(fitint):last(fitint)

# title for the plot
m = quarkmasses(hdf5file,hdf5groupN)
G = gaugegroup(hdf5file,hdf5groupN)
β = couplingβ(hdf5file,hdf5groupN)
title=L"$ %$(T)\times %$(L)^3, \beta=%$β, m_q=%$(first(m))$"
name="$(T)x$(L)^3 $G, m=$(first(m)), smeared"
state = L"\eta'"

t = 2:11
#scatter(t,meffC[t],yerr=ΔmeffC[t],label="non-singlet")
scatter(t,meffCD[t],yerr=ΔmeffCD[t],markershape=:circle;label=L"$~~$wall source: no subtraction, no smearing",legend_cell_align = "left",extra_kwargs = :subplot)
scatter!(t,meffIm[t],yerr=ΔmeffIm[t],markershape=:rect,label=L"$~~$wall source: excited state subtraction",legend_cell_align = "left",extra_kwargs = :subplot)
scatter!(t,meffCS[t],yerr=ΔmeffCS[t],markershape=:pentagon,ms=5,label=L"$~~$smeared connected piece",legend_cell_align = "left",extra_kwargs = :subplot)
plot!(t,mπ*ones(length(t)),ribbon=Δmπ, label=L"$~~m_\pi$",markershape=:none,legend_cell_align = "left",extra_kwargs = :subplot)
plot!(t,mρ*ones(length(t)),ribbon=Δmρ,label=L"$~~m_\rho$",markershape=:none,legend_cell_align = "left",extra_kwargs = :subplot)
plot!(fitη,mη*ones(length(fitη)),ribbon=Δmη, label=L"$~~m_{\eta'}$",markershape=:none,legend_cell_align = "left",extra_kwargs = :subplot)
plot!(ylims=(0.4,1.),xticks=1:T)
plot!(title=title,xlabel=L"t",ylabel="effective mass")
savefig("output/figures/smeared.pdf")
plot!(legend=:outerright)
