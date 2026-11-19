# Unit test 11: graded remeshing (GradedRemesh), the :graded option of testProp.
# Small tip elements (hf) along the crack front, one equilateral rim triangle per front edge,
# growing to h inside. Tested on a flat disc and on a dish (spherical cap) so curvature is covered.
#  Part 1 (no CGAL needed, runs in Pkg.test):
#  - GradedRemeshResampleFront: front points every ~hf, the original front points kept, all of
#    them on the true rim circle (the spline through the coarse rim, not its straight edges).
#  - GradedRemeshMLSProject: points on the dish's flat facets go back onto the sphere.
#  - GradedRemeshInsertPoint!: splitting a triangle keeps the area and adds two triangles.
#  Part 2 (needs BuildCGAL with the graded_remeshing tool; skipped otherwise):
#  - GradedRemesh: every front edge has its own equilateral rim triangle, no bad triangles
#    anywhere, the front and the surface stay on the true shape, and CleanEdgeTris +
#    IsosceliseEdgeTris then change nothing (so it fits the propagation loop as is).
#  - The same on a wavy, kinked front remeshed as the propagation loop does (FrontCleanup=true):
#    own rim triangle and inner point per front edge, front kept, clean-up changes nothing.
# Figure: the flat and dish meshes from part 2 (rim triangles in red) and the front spacing.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

#Coarse disc of radius 1 made of nr rings, n points on the outer ring, optionally lifted to
#z(x,y). Neighbouring rings are zipped together by angle. Returns Points,Triangles.
function GradedTestDisc(n,nr;z=(x,y)->0.0)
	X=[[0.0,0.0]]; Rings=Vector{Vector{Int}}(); Ang=Vector{Vector{Float64}}()
	for k=1:nr
		m=max(6,round(Int,n*k/nr)); r=k/nr
		t=[2π*(j+0.5*isodd(k))/m for j=0:m-1]
		push!(Rings,collect(length(X)+1:length(X)+m)); push!(Ang,t)
		append!(X,[[r*cos(v),r*sin(v)] for v in t])
	end
	T=[[1,Rings[1][j],Rings[1][mod1(j+1,length(Rings[1]))]] for j in eachindex(Rings[1])]
	for k=1:nr-1   #zip ring k (inner, A) and ring k+1 (outer, B), anticlockwise
		A=Rings[k]; B=Rings[k+1]; na=length(A); nb=length(B)
		ta(i)=Ang[k][mod1(i+1,na)]+2π*fld(i,na); tb(j)=Ang[k+1][mod1(j+1,nb)]+2π*fld(j,nb)
		i=0; j=0
		#start the outer ring at the first point after the inner ring's first point
		while tb(j)<ta(0); j+=1; end
		j0=j
		while i<na || j<j0+nb
			if i<na && (j==j0+nb || ta(i+1)<tb(j+1))
				push!(T,[A[mod1(i+1,na)],B[mod1(j,nb)],A[mod1(i+2,na)]]); i+=1
			else
				push!(T,[A[mod1(i+1,na)],B[mod1(j,nb)],B[mod1(j+1,nb)]]); j+=1
			end
		end
	end
	P=[collect(1:length(X)) permutedims(reduce(hcat,X)) [z(p[1],p[2]) for p in X]]
	return P,permutedims(reduce(hcat,T))
end
GradedRemeshArea(P,T)=TriAreas(CutAndDisplaceJulia.CreateP1P2P3(T,P)...)[1]
DishZ(x,y)=1.5-sqrt(1.5^2-x^2-y^2)          #spherical cap, radius 1.5, rim radius 1
DishError(X)=maximum(abs(norm(X[i,:]-[0,0,1.5])-1.5) for i in 1:size(X,1))
TriMinAngles(X,T)=[minimum(begin a,b,c=X[t[i],:],X[t[mod1(i+1,3)],:],X[t[mod1(i+2,3)],:]
	acosd(clamp(dot(b-a,c-a)/norm(b-a)/norm(c-a),-1,1)) end for i=1:3) for t in eachrow(T)]

@testset "Graded remeshing" begin

	h=0.2; hf=0.05
	(Pd,Td)=GradedTestDisc(32,5)                #flat disc, rim edges ≈0.2
	(Pc,Tc)=GradedTestDisc(32,5;z=DishZ)        #the same on a dish
	@test all(GradedRemeshArea(Pd,Td).>0)       #(a valid test mesh)

	#Part 1a: front resampling
	for (P,T,curved) in ((Pd,Td,false),(Pc,Tc,true))
		X0=P[:,2:4]
		(X,Tn,Loops)=CutAndDisplaceJulia.GradedRemeshResampleFront(X0,T,hf)
		@test length(Loops)==1
		F=reduce(hcat,[X[i] for i in Loops[1]])'
		Lf=norm.(eachrow(F-circshift(F,-1)))
		@test 0.75hf<minimum(Lf) && maximum(Lf)<1.25hf
		@test all(i in Loops[1] for i in CutAndDisplaceJulia.GradedRemeshBorderLoops(T)[1])
		#on the true rim (radius 1 in map view; on the sphere for the dish). The coarse rim's
		#straight edges are up to 1-cos(π/32) ≈ 5e-3 inside the circle.
		@test maximum(abs.(hypot.(F[:,1],F[:,2]).-1))<5e-4
		curved && @test DishError(F)<5e-4
		@test length(Tn)==size(T,1)+size(F,1)-length(CutAndDisplaceJulia.GradedRemeshBorderLoops(T)[1])
	end

	#Part 1b: smooth-surface projection on the dish
	X0=Pc[:,2:4]; N0=CutAndDisplaceJulia.GradedRemeshVertexNormals(X0,Tc)
	R=2.5*CutAndDisplaceJulia.GradedRemeshMeanEdgeLength(X0,Tc)
	Cent=[vec(mean(X0[t,:],dims=1)) for t in eachrow(Tc)]          #on the flat facets
	Proj=[CutAndDisplaceJulia.GradedRemeshMLSProject(c,X0,N0,R) for c in Cent]
	@test DishError(reduce(hcat,Cent)')>5e-3                       #facets are off the sphere...
	@test DishError(reduce(hcat,Proj)')<1e-3                       #...the projection is on it
	Xd=Pd[:,2:4]; Nd=CutAndDisplaceJulia.GradedRemeshVertexNormals(Xd,Td)
	@test all(abs(CutAndDisplaceJulia.GradedRemeshMLSProject([x,y,0.0],Xd,Nd,R)[3])<1e-12
	          for (x,y) in ((0.1,0.2),(-0.5,0.3),(0.8,-0.4)))      #a plane stays a plane

	#Part 1c: inserting a point
	X=[Vector(r) for r in eachrow(Xd)]; T=[Vector(r) for r in eachrow(Td)]
	A0=sum(GradedRemeshArea(Pd,Td))
	CutAndDisplaceJulia.GradedRemeshInsertPoint!(X,T,[0.31,0.17,0.0])
	P2=[collect(1:length(X)) reduce(hcat,X)']; T2=permutedims(reduce(hcat,T))
	@test size(T2,1)==size(Td,1)+2
	@test sum(GradedRemeshArea(P2,T2))≈A0 rtol=1e-12

	#Part 2: the whole graded remesh (CGAL)
	HasTool=isdefined(Main,:BuildCGAL) && isdefined(BuildCGAL,:GradedRemeshingCGAL) &&
	        isfile(joinpath(BuildCGAL.tools_dir(),Sys.iswindows() ? "graded_remeshing.exe" : "graded_remeshing"))
	if !HasTool
		@info "BuildCGAL with the graded_remeshing tool not loaded: skipping the CGAL part of the graded remeshing test"
	else
		Results=[]
		for (name,P,T,curved) in (("Flat disc",Pd,Td,false),("Dish",Pc,Tc,true))
			(Pg,Tg)=CutAndDisplaceJulia.GradedRemesh(P,T,h,hf,BuildCGAL.GradedRemeshingCGAL;growth=1.3)
			X=Pg[:,2:4]
			L=CutAndDisplaceJulia.GradedRemeshBorderLoops(Tg)[1]
			FrontEdges=[(L[k],L[mod1(k+1,length(L))]) for k in eachindex(L)]
			RimTri=[findfirst(r->(a in r)&&(b in r),collect(eachrow(Tg))) for (a,b) in FrontEdges]
			@test allunique(RimTri)                                   #one rim triangle per front edge...
			Inner=[first(setdiff(Tg[t,:],[a,b])) for ((a,b),t) in zip(FrontEdges,RimTri)]
			@test allunique(Inner)                                    #...each with its own inner point
			MinA=TriMinAngles(X,Tg)
			@test all(MinA[RimTri].>58)                               #equilateral tip triangles
			@test minimum(MinA)>20                                    #no poor triangles anywhere
			Lf=[norm(X[a,:]-X[b,:]) for (a,b) in FrontEdges]
			@test 0.75hf<minimum(Lf) && maximum(Lf)<1.25hf
			@test maximum(abs.(hypot.(X[L,1],X[L,2]).-1))<5e-4       #front on the true rim
			curved ? (@test DishError(X)<1e-3) : (@test all(X[:,3].==0))
			@test sum(GradedRemeshArea(Pg,Tg))≈sum(GradedRemeshArea(P,T)) rtol=0.02
			#CleanEdgeTris + IsosceliseEdgeTris (as in the propagation loop) change nothing
			(Q1,Q2,Q3,Tc2,Pc2,MP,FNV)=CutAndDisplaceJulia.CleanEdgeTris(Pg,Tg)
			(S1,S2,S3,Ti2,Pi2,MP2,FNV2)=CutAndDisplaceJulia.IsosceliseEdgeTris(MP,Q1,Q2,Q3,Tc2,Pc2,FNV)
			@test size(Ti2,1)==size(Tg,1)
			#(IsosceliseEdgeTris works in each triangle's own plane, the apex points sit on the
			# curved surface, so on the dish it may nudge them by a tiny fraction of hf)
			@test maximum(abs.(sort(vec(Pi2[:,2:4]))-sort(vec(Pg[:,2:4]))))<0.01hf
			push!(Results,(name,Pg,Tg,RimTri,Lf))
		end

		#A wavy, kinked front like the propagation loop's (the disc's rim pushed in and out), remeshed as the
		#loop does (FrontCleanup=true). Every front edge must still get its own rim triangle and inner point,
		#the front must stay on the input front (the spline through it only rounds its kinks), and the rim
		#clean-up must then change nothing, so no front point is dropped. (Before this was fixed, 8 front
		#edges here shared inner points, the clean-up removed 8 triangles and the front moved by 0.63hf.)
		(Pw,Tw)=GradedTestDisc(30,4)
		θw=atan.(Pw[:,3],Pw[:,2]); fw=1 .+ 0.12 .* hypot.(Pw[:,2],Pw[:,3]) .* (sin.(7θw) .+ 0.4 .* sin.(13θw))
		Pw[:,2].*=fw; Pw[:,3].*=fw
		hw=0.2; hfw=0.1
		(Pg,Tg)=CutAndDisplaceJulia.GradedRemesh(Pw,Tw,hw,hfw,BuildCGAL.GradedRemeshingCGAL;growth=1.3,FrontCleanup=true)
		X=Pg[:,2:4]
		L=CutAndDisplaceJulia.GradedRemeshBorderLoops(Tg)[1]
		FrontEdges=[(L[k],L[mod1(k+1,length(L))]) for k in eachindex(L)]
		RimTri=[findfirst(r->(a in r)&&(b in r),collect(eachrow(Tg))) for (a,b) in FrontEdges]
		@test allunique(RimTri)
		@test allunique([first(setdiff(Tg[t,:],[a,b])) for ((a,b),t) in zip(FrontEdges,RimTri)])
		Lw=CutAndDisplaceJulia.GradedRemeshBorderLoops(Tw)[1]; Xw=Pw[:,2:4]
		SegDist(p,a,b)=(t=clamp(dot(p-a,b-a)/dot(b-a,b-a),0,1); norm(p-(a+t*(b-a))))
		FrontOff=maximum(minimum(SegDist(X[i,:],Xw[Lw[k],:],Xw[Lw[mod1(k+1,length(Lw))],:]) for k in eachindex(Lw)) for i in L)
		@test FrontOff<0.4hfw                                         #measured 0.2hf
		@test minimum(TriMinAngles(X,Tg))>10                          #measured 13°: the sharpest kinks
		(Q1,Q2,Q3,Tc2,Pc2,MP,FNV)=CutAndDisplaceJulia.CleanEdgeTris(Pg,Tg)
		(S1,S2,S3,Ti2,Pi2,MP2,FNV2)=CutAndDisplaceJulia.IsosceliseEdgeTris(MP,Q1,Q2,Q3,Tc2,Pc2,FNV)
		#(at the sharpest kinks IsosceliseEdgeTris may split a rim triangle rather than squeeze its neighbours, so the
		# triangle count can grow; what must hold is that the front is untouched and every front edge then has an
		# isosceles rim triangle with an inner point of its own)
		Lc=CutAndDisplaceJulia.GradedRemeshBorderLoops(Ti2)[1]; Xc=Pi2[:,2:4]
		@test length(Lc)==length(L)
		@test maximum(minimum(norm(Xc[i,:]-X[j,:]) for j in L) for i in Lc)<1e-9    #the same front points
		FE2=[(Lc[k],Lc[mod1(k+1,length(Lc))]) for k in eachindex(Lc)]
		RT2=[findfirst(r->(a in r)&&(b in r),collect(eachrow(Ti2))) for (a,b) in FE2]
		In2=[first(setdiff(Ti2[t,:],[a,b])) for ((a,b),t) in zip(FE2,RT2)]
		@test allunique(In2)
		@test all(abs(norm(Xc[c,:]-Xc[a,:])-norm(Xc[c,:]-Xc[b,:]))<1e-6*norm(Xc[b,:]-Xc[a,:]) for ((a,b),c) in zip(FE2,In2))
		push!(Results,("Wavy front (as in the loop)",Pg,Tg,RimTri,[norm(X[a,:]-X[b,:]) for (a,b) in FrontEdges]))

		if MakeFigures
			fig=CM.Figure(size=(2000,560))
			for (k,(name,Pg,Tg,RimTri,Lf)) in enumerate(Results)
				(A,B,C)=CutAndDisplaceJulia.CreateP1P2P3(Tg,Pg)
				wavy=k==3
				ax=CM.Axis(fig[1,k],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
					title=wavy ? "$name: $(size(Tg,1)) triangles, h=$hw, hf=$hfw" : "$name: $(size(Tg,1)) triangles, h=$h at the centre, hf=$hf on the front")
				DrawTriEdges!(ax,A,B,C)
				DrawTriEdges!(ax,A[RimTri,:],B[RimTri,:],C[RimTri,:];color=:red,linewidth=0.8)
				wavy ? (CM.xlims!(ax,-1.25,1.25); CM.ylims!(ax,-1.25,1.25)) : (CM.xlims!(ax,0.45,1.05); CM.ylims!(ax,-0.3,0.3))
			end
			ax=CM.Axis(fig[1,4],title="Front edge lengths, disc and dish (target hf)",xlabel="edge length",ylabel="count")
			for (k,(name,Pg,Tg,RimTri,Lf)) in enumerate(Results[1:2])
				CM.hist!(ax,Lf;bins=range(0.7hf,1.3hf,length=25),label=name,color=(k==1 ? :steelblue : :orange,0.6))
			end
			CM.vlines!(ax,[hf];color=:black,linestyle=:dash)
			CM.axislegend(ax)
			CM.Label(fig[0,:],"GradedRemesh: one rim triangle (red) per front edge, equilateral on smooth fronts; sizes grow by 1.3 per row inwards (disc and dish: close-up of the right-hand rim)",fontsize=16)
			SaveTestFigure(fig,"ut_11_GradedRemesh")
		end
	end
end
