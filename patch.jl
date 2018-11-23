function patch(verts, faces ,fc)
#https://github.com/krcools/CompScienceMeshes.jl/blob/master/examples/plotlyjs_patches.jl
    v = verts
    c = faces

    x = v[:,1]; y = v[:,2]; z = v[:,3]
    i = c[:,1].-1; j = c[:,2].-1; k = c[:,3].-1

	#=
    if fcr == nothing
        #a = [barytocart(chart(Γ,cells(Γ,i)),[1,1]/3)[3] for i in 1:numcells(Γ)]
        #a = [cartesian(center(chart(Γ,cell)))[3] for cell in cells(Γ)]
    else
        a = fcr
    end

    m, M = extrema(a)
    if isapprox(m, M)
        n = ones(Integer, a)
    else
        n = floor.(Integer, (a.-m)./(M.-m)*(length(cm).-1)).+1
    end
    fc = [cm[i] for i in n]
	=#
	
    s = PlotlyJS.mesh3d(;
        x=x,
        y=y,
        z=z,
        i=i,
        j=j,
        k=k,
        facecolor=fc,
        colorbar=PlotlyJS.attr(title="z"),
    )

    #PlotlyJS.plot([s])
    return(s)
end