using Plots
using LaTeXStrings
using Colors
using DelimitedFiles
using HiRepAnalysis
pgfplotsx(legendfontsize=16,labelfontsize=20,tickfontsize=14,titlefontsize=20,ms=8,framestyle=:box,size=(700,500))

data = readdlm("output/data_flow.csv", ';',header=true)[1]
data = replace(data,"-"=>NaN)

parametric_colour(name,β;N=40,n1=10,n2=10)=reverse(colormap(name,N; mid=0.55, logscale=false, b=0.1 ,w=0.5))[n1:end-n2][(Int ∘ round)((2β*10)%(N-n1-n2)+1)]
blue(β) = parametric_colour("Blues",β)
green(β) = parametric_colour("Greens",β)
orange(β) = parametric_colour("Oranges",β)

β    = data[:,1]
mπ   = data[:,12]
Δmπ  = data[:,13]
mρ   = data[:,14]
Δmρ  = data[:,15]
mη   = data[:,18]
Δmη  = data[:,19]
w0_c  = data[:,28]
Δw0_c = data[:,29]

Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
Δproduct(x,y,Δx,Δy) = sqrt((Δx * y)^2 + (Δy*x)^2)

scale  = w0_c
Δscale = Δw0_c
mπ   = @. mπ * scale
mρ   = @. mρ * scale
mη   = @. mη * scale
Δmπ  = @. Δproduct(mπ,scale,Δmπ,Δscale)
Δmρ  = @. Δproduct(mρ,scale,Δmρ,Δscale)
Δmη  = @. Δproduct(mη,scale,Δmη,Δscale)

i1, i2 = 1:6, 7:18
plt  = plot(title=latexstring(L"Sp(4) : ~ m_\pi L > 6"),ylabel=L" m_{\rm meson}  w_0 ",xlabel=L"m_\pi w_0")
scatter!(plt,mπ[i1],mρ[i1],xerr=Δmπ[i1],yerr=Δmρ[i1],label=L"\rho (\beta  = %$(first(β[i1])))", markershape=:circle,markercolor=green(first(β[i1])))
scatter!(plt,mπ[i2],mρ[i2],xerr=Δmπ[i2],yerr=Δmρ[i2],label=L"\rho (\beta  = %$(first(β[i2])))", markershape=:circle,markercolor=green(first(β[i2])))
scatter!(plt,mπ[i1],mη[i1],xerr=Δmπ[i1],yerr=Δmη[i1],label=L"\eta (\beta' = %$(first(β[i1])))",markershape=:diamond,markercolor=orange(first(β[i1])),ms=10)
scatter!(plt,mπ[i2],mη[i2],xerr=Δmπ[i2],yerr=Δmη[i2],label=L"\eta (\beta' = %$(first(β[i2])))",markershape=:diamond,markercolor=orange(first(β[i2])),ms=10)
# draw reference line of pions
lims_x, lims_y = xlims(plt), ylims(plt)
plot!(plt,[0.2,0.7],[0.2,0.7],label=L"\pi",color=blue(first(β[i1])),xlims=lims_x,ylims=lims_y)
plot!(plt,legend=:bottomright)
savefig("output/figures/eta_Sp4_w0.pdf")

