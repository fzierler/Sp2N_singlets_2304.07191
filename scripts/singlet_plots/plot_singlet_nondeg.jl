using Plots
using LaTeXStrings
using Colors
using DelimitedFiles
using HiRepAnalysis
pgfplotsx(legendfontsize=20,labelfontsize=20,tickfontsize=14,titlefontsize=20,ms=8,framestyle=:box,size=(500,500))

parametric_colour(name,β;N=40,n1=10,n2=10)=reverse(colormap(name,N; mid=0.55, logscale=false, b=0.1 ,w=0.5))[n1:end-n2][(Int ∘ round)((2β*10)%(N-n1-n2)+1)]
blue(β) = parametric_colour("Blues",β)
green(β) = parametric_colour("Greens",β)
orange(β) = parametric_colour("Oranges",β)

data = readdlm("output/data/data_non_deg.csv", ';',header=true)[1]
# AWI mass at degeneracy
mAWIdeg, ΔmAWIdeg = data[end,30],data[end,31]
# rest of the data
β    = data[:,1]
o = 10
mπ0  = data[:,6+o]
Δmπ0 = data[:,7+o]
mπf  = data[:,10+o]
Δmπf = data[:,11+o]
mρ0  = data[:,12+o]
Δmρ0 = data[:,13+o]
mρf  = data[:,14+o]
Δmρf = data[:,15+o]
mη   = data[:,18+o]
Δmη  = data[:,19+o]
mAWI = data[:,20+o]
ΔmAWI= data[:,21+o]
# obtain AWI ratio
rq  = @. 2mAWI/mAWIdeg - 1
Δrq = @. 2sqrt((ΔmAWI/mAWIdeg)^2 + (ΔmAWIdeg*mAWI/mAWIdeg^2)^2)

Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
Δproduct(x,y,Δx,Δy) = sqrt((Δx * y)^2 + (Δy*x)^2)

rπf = mπf ./ mπ0
rρ0 = mρ0 ./ mπ0
rρf = mρf ./ mπ0
rη = mη ./ mπ0
Δrπf = Δratio.(mπf,mπ0,Δmπf,Δmπ0)
Δrρ0 = Δratio.(mρ0,mπ0,Δmρ0,Δmπ0)
Δrρf = Δratio.(mρf,mπ0,Δmρf,Δmπ0)
Δrη = Δratio.(mη,mπ0,Δmη,Δmπ0)

x1, Δx1 = rq, Δrq
x2, Δx2 = mπf ./ mπ0, Δratio.(mπf,mπ0,Δmπf,Δmπ0)

plt = plot(legend=:outerright,ylabel=L" m_{\rm meson}/m_{\pi^0}",title=L"Sp(4),~ m_u \neq m_d")
plt_pi   = plot(plt,xlabel=L" m_{\pi^\pm} / m_{\pi^0}")
plt_pcac = plot(plt,xlabel=L"PCAC mass ratio $m_d/m_u$")
plt_direct =  plot(legend=:outerright,ylabel=L" a m_{\rm meson}",title=L"Sp(4),~ m_u \neq m_d",xlabel=L"am_{\pi^\pm}")

# obtain minimal and maximal pion-ratio
xcont1 = [1,maximum(x1+2Δx1)]
scatter!(plt_pcac,x1,xerr=Δx1,rπf,yerr=Δrπf,label=L"{\pi^\pm}",markershape=:rect,ms=7)
plot!(plt_pcac,xcont1,ones(2),label=L"{\pi^0}",lw=3)
scatter!(plt_pcac,x1,xerr=Δx1,rη ,yerr=Δrη,label=L"{\eta'}",markershape=:circle,ms=7)
scatter!(plt_pcac,x1,xerr=Δx1,rρf,yerr=Δrρf,label=L"{\rho^\pm}",markershape=:dtriangle,ms=8)
scatter!(plt_pcac,x1,xerr=Δx1,rρ0,yerr=Δrρ0,label=L"{\rho^0}",markershape=:diamond,ms=8)

xcont2 = [1,maximum(x2+2Δx2)]
plot!(plt_pi,xcont2,xcont2,label=L"{\pi^\pm}",lw=3)
plot!(plt_pi,xcont2,ones(2),label=L"{\pi^0}",lw=3)
scatter!(plt_pi,x2,xerr=Δx2,rη ,yerr=Δrη ,label=L"{\eta'}",markershape=:circle,ms=7)
scatter!(plt_pi,x2,xerr=Δx2,rρf,yerr=Δrρf,label=L"{\rho^\pm}",markershape=:dtriangle,ms=8)
scatter!(plt_pi,x2,xerr=Δx2,rρ0,yerr=Δrρ0,label=L"{\rho^0}",markershape=:diamond,ms=8)

extrema_mπf = [extrema(mπf)...]
plot!(plt_direct,extrema_mπf,extrema_mπf,label=L"{\pi^\pm}",lw=3)
scatter!(plt_direct, mπf, xerr=Δmπf, mπ0, yerr=Δmπ0, label=L"{\pi^0}",markershape=:circle,ms=7)
scatter!(plt_direct, mπf, xerr=Δmπf, mη , yerr=Δmη , label=L"{\eta'}",markershape=:rect,ms=7)
scatter!(plt_direct, mπf, xerr=Δmπf, mρ0, yerr=Δmρ0, label=L"{\rho^0}",markershape=:dtriangle,ms=8)
scatter!(plt_direct, mπf, xerr=Δmπf, mρf, yerr=Δmρf, label=L"{\rho^\pm}",markershape=:diamond,ms=8)

savefig(plt_pcac,"output/figures/Sp4_eta_nondeg_AWI.pdf")
savefig(plt_pi,"output/figures/Sp4_eta_nondeg.pdf")
savefig(plt_direct,"output/figures/Sp4_eta_nondeg_not_a_ratio.pdf")
plt_pi
plt_pcac
plt_direct