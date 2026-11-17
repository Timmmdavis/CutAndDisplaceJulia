# Unit test 8: the per-step crack figures (CutAndDisplaceJuliaPlots: CrackViews.jl).
# Needs CutAndDisplaceJuliaPlots; skipped inside Pkg.test of CutAndDisplaceJulia.
# A made-up "run" of three steps (a vertical penny, radius ~500 m, rising and growing) is written in
# the same text format the propagation test writes (MeshWithDisp / MeshEdges .txt files). Checks:
#  - FindSteps finds the steps of the run in step order; ReadStep reads back what was written
#  - ViewFrame gets the dip right (vertical crack, and the same crack rotated to strike 60°, dip 60°)
#  - CrackWalls opens the walls by exactly Exag·Dn along the normal, closed at the rim
#  - PlotCrackState draws, SaveFig writes a PNG, RenderRun renders every step
# Figure: the input (opening, tip K/Kc) next to the PlotCrackState figure it produced, for the
# vertical crack and for the rotated one.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

if !isdefined(Main,:CutAndDisplaceJuliaPlots)
	@info "CutAndDisplaceJuliaPlots not loaded: skipping the plotting test (run RunUnitTests.jl to include it)"
else
@testset "Plotting: crack views" begin

	CV=CutAndDisplaceJuliaPlots
	(Points,Triangles,Q1,Q2,Q3,FaceNormalVector,MidPoint)=LoadSmallPennyMesh()
	#flat penny (x,y,0) -> vertical penny in the x-z plane, radius R, centre depth Zc
	ToWorld(P,R,Zc)=[R.*P[:,1] zeros(size(P,1)) R.*P[:,2].+Zc]
	OpeningAt(P1,P2,P3,R,Zc)=begin
		M=(P1.+P2.+P3)./3
		r=hypot.(M[:,1],M[:,3].-Zc)./R
		0.5.*sqrt.(max.(1 .-r.^2,0))                    #elliptical opening, 0.5 m at the centre
	end

	#Write the made-up run: steps 1, 2 and 10 (to check numeric, not text, order)
	Dir=mktempdir()
	Kc=5e7
	Written=Dict{Int,Any}()
	for (i,R,Zc) in ((10,600.0,-1500.0),(1,400.0,-2000.0),(2,450.0,-1900.0))
		P1=ToWorld(Q1,R,Zc); P2=ToWorld(Q2,R,Zc); P3=ToWorld(Q3,R,Zc)
		Dn=OpeningAt(P1,P2,P3,R,Zc)
		Rim=BoundaryEdgesBruteForce(P1,P2,P3)
		TipXYZ=reduce(vcat,[((a.+b)./2)' for (a,b) in Rim])
		TipK=Kc.*(1 .+0.5.*(TipXYZ[:,3].-Zc)./R)               #K/Kc from 0.5 (bottom) to 1.5 (top)
		open(joinpath(Dir,"$i-2-MeshWithDisp-ForDrawing-7.txt"),"w") do io
			println(io,"P1x P1y P1z P2x P2y P2z P3x P3y P3z Dn Dss Dds")
			for k=1:size(P1,1)
				println(io,join(string.([P1[k,:];P2[k,:];P3[k,:];Dn[k];0.0;0.0])," "))
			end
		end
		open(joinpath(Dir,"$i-2-MeshEdges-ForDrawing-7.txt"),"w") do io
			println(io,"X Y Z StrainEnergy KCrit")
			for k=1:size(TipXYZ,1)
				println(io,join(string.([TipXYZ[k,:];TipK[k];Kc])," "))
			end
		end
		Written[i]=(P1=P1,P2=P2,P3=P3,Dn=Dn,TipXYZ=TipXYZ,TipK=TipK)
	end

	#Reading the run back
	Steps=CV.FindSteps(Dir)
	@test [s.i for s in Steps]==[1,2,10]
	d=CV.ReadStep(Steps[3])
	w=Written[10]
	@test d.P1≈w.P1 && d.P2≈w.P2 && d.P3≈w.P3
	@test d.Dn≈w.Dn
	@test d.TipXYZ≈w.TipXYZ && d.TipK≈w.TipK
	@test d.Kc==Kc

	#Strike and dip of the crack
	F=CV.ViewFrame(d.P1,d.P2,d.P3)
	@test F.dip ≈ 90 atol=1e-6
	Rot=[cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0 0 1]*[1 0 0; 0 cosd(30) -sind(30); 0 sind(30) cosd(30)]
	c=vec(mean(vcat(d.P1,d.P2,d.P3),dims=1))
	rot(P)=((P.-c')*Rot').+c'
	FR=CV.ViewFrame(rot(d.P1),rot(d.P2),rot(d.P3))
	@test FR.dip ≈ 60 atol=1e-6
	@test F.az ≈ 90 atol=1e-6                 #strike
	@test FR.az ≈ 60 atol=1e-6

	#Crack walls: separation = Exag*Dn (vertex average), zero on the rim
	Exag=100.0
	(ids,X,W1,W2)=CV.CrackWalls(d.P1,d.P2,d.P3,d.Dn,Exag,F.n̂)
	(ids2,vd)=CV.VertexDn(d.P1,d.P2,d.P3,d.Dn)
	Sep=sqrt.(sum((W1.-W2).^2,dims=2))
	@test vec(Sep)≈Exag.*vd
	RimV=unique(vcat([[e[1],e[2]] for e in CV.BoundaryEdges(ids)]...))
	@test all(Sep[RimV].==0)
	@test all(abs.((W1.-W2)*F.n̂).≈vec(Sep))                       #along the crack normal

	#Drawing and saving
	fig=CV.PlotCrackState(d.P1,d.P2,d.P3,d.Dn;TipXYZ=d.TipXYZ,TipK=d.TipK,Kc=d.Kc,meta=(step=10,V=1.0e5))
	@test fig isa CM.Figure
	Png=joinpath(Dir,"step10.png")
	CV.SaveFig(Png,fig)
	@test isfile(Png) && filesize(Png)>20_000
	figR=CV.PlotCrackState(rot(d.P1),rot(d.P2),rot(d.P3),d.Dn;TipXYZ=rot(d.TipXYZ),TipK=d.TipK,Kc=d.Kc,meta=(step="rotated",))
	PngR=joinpath(Dir,"rotated.png")
	CV.SaveFig(PngR,figR)
	@test isfile(PngR)
	times=CV.RenderRun(Dir)
	@test length(times)==3
	@test sort(readdir(joinpath(Dir,"frames")))==["frame_0001.png","frame_0002.png","frame_0010.png"]

	if MakeFigures
		fig=CM.Figure(size=(1900,1500))
		for (col,(name,P1,P2,P3,TipXYZ,Png,Fr)) in enumerate((("vertical crack (strike 90°, dip 90°)",d.P1,d.P2,d.P3,d.TipXYZ,Png,F),
		                                              ("same crack rotated (strike 60°, dip 60°)",rot(d.P1),rot(d.P2),rot(d.P3),rot(d.TipXYZ),PngR,FR)))
			#input drawn in the crack's own plane coordinates
			c0=vec(mean(vcat(P1,P2,P3),dims=1))
			U(P)=(P.-c0')*Fr.û; V(P)=(P.-c0')*Fr.v̂
			A1=[U(P1) V(P1) zeros(size(P1,1))]; A2=[U(P2) V(P2) zeros(size(P2,1))]; A3=[U(P3) V(P3) zeros(size(P3,1))]
			g=CM.GridLayout(fig[1,col])
			ax=CM.Axis(g[1,1],aspect=CM.DataAspect(),xlabel="along strike [m]",ylabel="up dip [m]",
				title="INPUT: $name, in the crack plane\nfill: opening Dn, dots: tip K/Kc\nViewFrame finds dip $(round(Fr.dip,digits=1))°")
			m=FillTris!(ax,A1,A2,A3,d.Dn,colormap=:viridis,colorrange=(0,maximum(d.Dn)))
			sc=CM.scatter!(ax,U(TipXYZ),V(TipXYZ),color=d.TipK./d.Kc,colormap=:coolwarm,colorrange=(0.5,1.5),markersize=12)
			CM.Colorbar(g[1,2],m,label="Dn [m]")
			CM.Colorbar(g[1,3],sc,label="K/Kc")
			axi=CM.Axis(fig[2,col],aspect=CM.DataAspect(),title="OUTPUT: PlotCrackState → SaveFig (the per-step figure)",titlesize=18)
			CM.image!(axi,rotr90(CM.load(Png)))
			CM.hidedecorations!(axi); CM.hidespines!(axi)
		end
		CM.rowsize!(fig.layout,1,CM.Relative(0.3))
		CM.Label(fig[0,:],"Unit test 8: crack figures (CutAndDisplaceJuliaPlots) from a made-up run: inputs on top, the figures PlotCrackState made below",fontsize=20)
		SaveTestFigure(fig,"ut_08_PlottingCrackViews")
	end
end
end
