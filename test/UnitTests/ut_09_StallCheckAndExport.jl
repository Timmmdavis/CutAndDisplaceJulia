# Unit test 9: the "has the crack stopped?" check and the crack-shape export (no CGAL needed).
#  - CheckIfInsidePreviousBoundary: true when every rim point of the new mesh lies within a
#    tube (radius 10% of the crack size) around the old rim, false once the crack has grown
#    past it. This is what stops the propagation loop. Checked with the same penny, a penny 5%
#    larger (inside the tube), and one 50% larger (outside).
#  - ExportCrackMesh: writes both crack walls as one closed .off surface, each wall moved by half
#    the displacement. Checks: twice the points and triangles; at every point the gap between the
#    walls equals the opening there (average of the touching triangles, 0 on the rim); the walls
#    are mirror images; reading the file back, the volume it encloses is the volume of that
#    opening, and close to Σ Area·Dn; shear moves the walls apart by Dss along the crack.
#  - ConnectedComponentsReader: turns CGAL's 0-based component labels into 1-based ones.
# Figure: the three rims against the tube; a cross-section of the exported walls against Dn.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Stall check and crack export" begin

	(Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint)=LoadSmallPennyMesh()
	(Area,HalfPerimeter)=TriAreas(P1,P2,P3)
	max_target_edge_length=maximum(HalfPerimeter)*(2/3)
	Scaled(s)=(P1.*s,P2.*s,P3.*s,MidPoint.*s)

	#Stall check
	Inside=Dict{Float64,Bool}()
	for s in (1.0,1.05,1.5)
		(A1,A2,A3,AM)=Scaled(s)
		Inside[s]=CutAndDisplaceJulia.CheckIfInsidePreviousBoundary(MidPoint,P1,P2,P3,FaceNormalVector,A1,A2,A3,AM,max_target_edge_length)
	end
	@test Inside[1.0]==true
	@test Inside[1.05]==true
	@test Inside[1.5]==false

	#Export of the opened crack
	r=hypot.(MidPoint[:,1],MidPoint[:,2])
	a=maximum(hypot.(Points[:,2],Points[:,3]))
	Dn=0.1.*sqrt.(max.(1 .-(r./a).^2,0))
	n=size(P1,1)
	(Tris,Pts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
	np=size(Pts,1)
	Dir=mktempdir()
	f=cd(Dir) do
		CutAndDisplaceJulia.ExportCrackMesh(P1,P2,P3,Dn,zeros(n),zeros(n),"OpenedCrack")
	end
	(PtsW,TrisW)=CutAndDisplaceJulia.OFFReader(f)
	@test size(TrisW,1)==2n
	@test size(PtsW,1)==2np
	Pos=PtsW[1:np,2:4]; Neg=PtsW[np+1:end,2:4]
	#Expected opening at each point, worked out here independently: mean Dn of the triangles
	#touching the point, 0 on the rim
	RimPts=Set(vcat([[Tuple(round.(p,digits=10)),Tuple(round.(q,digits=10))] for (p,q) in BoundaryEdgesBruteForce(P1,P2,P3)]...))
	DnPt=zeros(np); Cnt=zeros(np)
	for i=1:n, k=1:3
		DnPt[Tris[i,k]]+=Dn[i]; Cnt[Tris[i,k]]+=1
	end
	DnPt./=Cnt
	OnRim=[Tuple(round.(Pts[k,2:4],digits=10)) in RimPts for k=1:np]
	DnPt[OnRim].=0
	Gap=Pos[:,3].-Neg[:,3]
	@test Pos[:,1:2]≈Pts[:,2:3] && Neg[:,1:2]≈Pts[:,2:3]       #pure opening: walls move along z only
	@test Pos[:,3]≈-Neg[:,3]                                    #mirror images about z = 0
	@test maximum(abs.(Gap.-DnPt))<1e-12                        #gap = opening (not 2x)
	#Read back: volume enclosed by the closed surface (divergence theorem)
	EnclosedVolume=abs(sum(dot(PtsW[TrisW[i,1],2:4],cross(PtsW[TrisW[i,2],2:4],PtsW[TrisW[i,3],2:4])) for i=1:2n))/6
	VolumeOfPointOpening=sum(Area[i]*mean(DnPt[Tris[i,:]]) for i=1:n)  #linear opening between points
	@test EnclosedVolume ≈ VolumeOfPointOpening rtol=1e-10
	@test EnclosedVolume ≈ sum(Area.*Dn) rtol=0.2                     #point averaging smooths the rim
	#Shear: walls slide apart by Dss along the crack (in the x-y plane here)
	DssTest=fill(0.05,n)
	f2=cd(Dir) do
		CutAndDisplaceJulia.ExportCrackMesh(P1,P2,P3,zeros(n),DssTest,zeros(n),"ShearedCrack")
	end
	(PtsS,~)=CutAndDisplaceJulia.OFFReader(f2)
	Slide=sqrt.(sum((PtsS[1:np,2:4].-PtsS[np+1:end,2:4]).^2,dims=2))
	@test maximum(abs.(vec(Slide)[.!OnRim].-0.05))<1e-12
	@test maximum(abs.(PtsS[:,4]))<1e-12                        #stays in the crack plane

	#Connected components file (as written by BuildCGAL.ConnectedComponentsCGAL)
	g=joinpath(Dir,"ConnectedComponents.txt")
	write(g,"0\n0\n1\n1\n1\n")
	@test CutAndDisplaceJulia.ConnectedComponentsReader(g)==[1,1,2,2,2]

	if MakeFigures
		fig=CM.Figure(size=(1500,620))
		ax1=CM.Axis(fig[1,1],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="CheckIfInsidePreviousBoundary\nold rim (black) with its tube (grey band, ±10% of crack size)\nnew rims: same → $(Inside[1.0]), ×1.05 → $(Inside[1.05]), ×1.5 → $(Inside[1.5])")
		Rim=BoundaryEdgesBruteForce(P1,P2,P3)
		Tube=0.1*2a
		t=range(0,2π,length=300)
		CM.band!(ax1,CM.Point2f.((a+Tube).*cos.(t),(a+Tube).*sin.(t)),CM.Point2f.((a-Tube).*cos.(t),(a-Tube).*sin.(t)),color=(:gray,0.25))
		for (s,c) in ((1.0,:black),(1.05,:blue),(1.5,:red))
			for (p,q) in Rim
				CM.lines!(ax1,s.*[p[1],q[1]],s.*[p[2],q[2]],color=c,linewidth=2)
			end
			CM.lines!(ax1,[NaN],[NaN],color=c,linewidth=2,label="×$s: inside = $(Inside[s])")
		end
		CM.axislegend(ax1,position=:rb,labelsize=11)
		rp=hypot.(Pts[:,2],Pts[:,3])
		ax2=CM.Axis(fig[1,2],xlabel="distance from centre r",ylabel="gap between the exported walls",
			title="ExportCrackMesh: gap between the two walls at each point\nline: the opening Dn given (2·Dn dashed: the old, doubled export)\nvolume read back $(round(EnclosedVolume,sigdigits=4)) vs Σ Area·Dn $(round(sum(Area.*Dn),sigdigits=4))")
		rr=range(0,a,length=200)
		CM.lines!(ax2,rr,0.1.*sqrt.(max.(1 .-(rr./a).^2,0)),color=:black,label="Dn (opening given)")
		CM.lines!(ax2,rr,0.2.*sqrt.(max.(1 .-(rr./a).^2,0)),color=:black,linestyle=:dash,label="2·Dn (old export)")
		CM.scatter!(ax2,rp,Gap,color=:blue,markersize=7,label="exported gap")
		CM.axislegend(ax2,position=:lb,labelsize=11)
		CM.Label(fig[0,:],"Unit test 9: the stop (stall) check and exporting the opened crack",fontsize=20)
		SaveTestFigure(fig,"ut_09_StallCheckAndExport")
	end
end
