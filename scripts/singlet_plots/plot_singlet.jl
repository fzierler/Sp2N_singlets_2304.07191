using Plots
using LaTeXStrings
using Colors
using DelimitedFiles
using HiRepAnalysis
pgfplotsx(legendfontsize=16,labelfontsize=20,tickfontsize=14,titlefontsize=20,framestyle=:box,size=(700,500))

# set up parametric colours for plotting different shades
parametric_colour(name,β;N=40,n1=10,n2=10)=reverse(colormap(name,N; mid=0.55, logscale=false, b=0.1 ,w=0.5))[n1:end-n2][(Int ∘ round)((2β*10)%(N-n1-n2)+1)]
blue(β) = parametric_colour("Blues",β)
green(β) = parametric_colour("Greens",β)
orange(β) = parametric_colour("Oranges",β)

#hardcode markershapes for FIGURE V 
function bw_shape(β) 
    β≈6.9 && return :circle 
    β≈7.2 && return :rect 
    β≈2   && return :pentagon
end
bw_color(β) = ifelse(Base.isgreater(β,6.8),:white ,:black)

data = readdlm("output/data/data_deriv.csv",';',skipstart=1)
β,mq,L,T = ( data[:,i] for i in 1:4 )
mπ,Δmπ,mρ,Δmρ,ma0,Δma0,mη,Δmη,mσ,Δmσ = ( data[:,i] for i in 12:21 )

# give everything in terms of pion mass
x   = @. mπ / mρ
mρπ = @. mρ / mπ
mηπ = @. mη / mπ
mσπ = @. mσ / mπ
mηρ = @. mη / mρ
# conservative error estimate
Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
Δx   = Δratio.(mπ,mρ,Δmπ,Δmρ)
Δmρπ = Δratio.(mρ,mπ,Δmρ,Δmπ)
Δmηπ = Δratio.(mη,mπ,Δmη,Δmπ)
Δmσπ = Δratio.(mσ,mπ,Δmσ,Δmπ)
Δmηρ = Δratio.(mη,mρ,Δmη,Δmρ)

labelη = L"\eta'"
labelσ = L"\sigma" 

# get indices of comparable sets of ensembles_degenerate
iβ1 = findall(x -> isapprox(6.9,x),β)
iβ2 = findall(x -> isapprox(7.2,x),β)
iSU2   = findall(x -> isapprox(2.0,x),β)

# labels for common axes
y_mesons = L" m_{\rm meson} / m_\pi "
x_ratio  = L"m_\pi / m_\rho"

# create initial plota
plt_σ_Sp4 = plot(legend=:topright,title = L"Sp(4) : ~ m_\pi L > 6",xlabel=x_ratio,ylabel=y_mesons)
plt_η_Sp4 = plot(legend=:topright,title = L"Sp(4) : ~ m_\pi L > 6",xlabel=x_ratio,ylabel=y_mesons)
plt_η_SU2 = plot(legend=:topright,title = L"SU(2) : ~ m_\pi L > 6",xlabel=x_ratio,ylabel=y_mesons)
plt_η_v_ρ_Sp4 = plot(legend=:topright,title = L"Sp(4) : ~ m_\pi L > 6",xlabel=x_ratio,ylabel=L" m_{\eta'} / m_\rho ")
# set up data for solid lines
min_x = minimum(hcat(x[iβ1],x[iβ1]))
max_x = maximum(hcat(x[iβ1],x[iβ1]))
xd = 0.95*min_x:0.001:max_x*1.05

# plot sigma meson for Sp(4)
plot!(plt_σ_Sp4,xd,inv.(xd),label=L"\rho",color=green(iβ1[1]),lw=3)
scatter!(plt_σ_Sp4,x[iβ1],xerr=Δx[iβ1],mσπ[iβ1],yerr=Δmσπ[iβ1],label=L"\sigma~(\beta = %$(first(β[iβ1])))",markershape=:circle,ms=8,markercolor=orange(first(β[iβ1])))   
scatter!(plt_σ_Sp4,x[iβ2],xerr=Δx[iβ2],mσπ[iβ2],yerr=Δmσπ[iβ2],label=L"\sigma~(\beta = %$(first(β[iβ2])))",markershape=:diamond,ms=10,markercolor=orange(first(β[iβ2])))   
plot!(plt_σ_Sp4,xd,ones(length(xd)),label=L"\pi",color=blue(iβ1[1]),lw=3)

# plot eta' meson for Sp(4)
plot!(plt_η_Sp4,xd,inv.(xd),label=L"\rho",color=green(iβ1[1]),lw=3)
scatter!(plt_η_Sp4,x[iβ1],xerr=Δx[iβ1],mηπ[iβ1],yerr=Δmηπ[iβ1],label=L"\eta'~(\beta = %$(first(β[iβ1])))",markershape=:circle,ms=8,markercolor=orange(first(β[iβ1])))   
scatter!(plt_η_Sp4,x[iβ2],xerr=Δx[iβ2],mηπ[iβ2],yerr=Δmηπ[iβ2],label=L"\eta'~(\beta = %$(first(β[iβ2])))",markershape=:diamond,ms=10,markercolor=orange(first(β[iβ2])))   
plot!(plt_η_Sp4,xd,ones(length(xd)),label=L"\pi",color=blue(iβ1[1]),lw=3)

# plot eta' meson vs vector meson for Sp(4)
scatter!(plt_η_v_ρ_Sp4,x[iβ1],xerr=Δx[iβ1],mηρ[iβ1],yerr=Δmηρ[iβ1],label=L"\beta = %$(first(β[iβ1]))",markershape=:circle,ms=8,markercolor=orange(first(β[iβ1])))   
scatter!(plt_η_v_ρ_Sp4,x[iβ2],xerr=Δx[iβ2],mηρ[iβ2],yerr=Δmηρ[iβ2],label=L"\beta = %$(first(β[iβ2]))",markershape=:diamond,ms=10,markercolor=orange(first(β[iβ2])))   
plot!(plt_η_v_ρ_Sp4,[0.6,1.0],ones(2),label="",ls=:dash,colour=:black)

# plot eta' meson for SU(2)
plot!(plt_η_SU2,xd,inv.(xd),label=L"\rho",color=green(iβ1[1]),lw=3)
scatter!(plt_η_SU2,x[iSU2],xerr=Δx[iSU2],mηπ[iSU2],yerr=Δmηπ[iSU2],label=L"\eta'~(\beta = %$(first(β[iSU2])))",markershape=:circle,ms=8,markercolor=orange(first(β[iβ1])))   
plot!(plt_η_SU2,xd,ones(length(xd)),label=L"\pi",color=blue(iβ1[1]),lw=3)

# set appropriate y-limits
plot!(plt_σ_Sp4,xlims=(0.64,0.9))
plot!(plt_η_Sp4,xlims=(0.67,0.9))
plot!(plt_η_SU2,xlims=(0.73,0.9))
plot!(plt_η_v_ρ_Sp4,xlims=(0.67,0.9),ylims=(0.8,1.1))

savefig(plt_η_Sp4,"output/figures/eta_Sp4.pdf")
savefig(plt_η_SU2,"output/figures/eta_SU2.pdf")
savefig(plt_σ_Sp4,"output/figures/sigma_Sp4.pdf")
savefig(plt_η_v_ρ_Sp4,"output/figures/Sp4_eta_in_rho.pdf")