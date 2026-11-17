# Unit test 10: dip/azimuth/strike of each triangle, and finding each triangle's neighbours.
# A unit penny bent into a trough, z = 0.5·x², so the dip changes across it:
#  - DipAndAzimuth: dip = angle of the triangle from horizontal, azimuth = compass direction it
#    dips towards (x = east, y = north), strike 90° anticlockwise from that. Checked against
#    the same angles worked out here from the normals, and the answer must not depend on which
#    way the normals point. On this trough: east half dips west (270°), west half dips east (90°).
#  - FindConnectedTriangles: the triangles sharing an edge with each triangle. Checked against
#    ConnectedConstraints (used in the meshing), and rim triangles have fewer than 3.
# Figure: the trough coloured by dip, by azimuth, and by number of neighbours.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Dip/azimuth and neighbours" begin

	(Points,Triangles,~,~,~,~,~)=LoadSmallPennyMesh()
	Points[:,4]=0.5.*Points[:,2].^2                         #bend into a trough along y
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
	(Triangles,Points)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)   #shared corners → shared points
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
	n=size(P1,1)

	#Dip, azimuth, strike
	(Dip,Azimuth,Strike)=CutAndDisplaceJulia.DipAndAzimuth(Triangles,Points,FaceNormalVector)
	Up=FaceNormalVector.*sign.(FaceNormalVector[:,3])               #normals made to point up
	DipHere=acosd.(clamp.(Up[:,3],-1,1))
	AzHere=mod.(atand.(Up[:,1],Up[:,2]),360)                        #horizontal part of an upward normal points down-dip
	@test maximum(abs.(Dip.-DipHere))<1e-8
	Steep=DipHere.>1                                                #azimuth is undefined on flat triangles
	AzDiff(a,b)=abs(mod(a-b+180,360)-180)
	@test maximum(AzDiff.(Azimuth[Steep],AzHere[Steep]))<1e-6
	@test maximum(AzDiff.(Strike[Steep],mod.(Azimuth[Steep].-90,360)))<1e-6
	East=(MidPoint[:,1].>0.2); West=(MidPoint[:,1].<-0.2)
	#(flat triangles on a curved surface tilt a little along y too: within 10°)
	@test all(AzDiff.(Azimuth[East],270).<10) && all(AzDiff.(Azimuth[West],90).<10)
	@test maximum(abs.(Dip.-atand.(abs.(MidPoint[:,1]))))<3        #≈ slope of z = 0.5x² at the triangle
	#Flipping every normal gives the same answer
	(Dip2,Azimuth2,Strike2)=CutAndDisplaceJulia.DipAndAzimuth(Triangles,Points,-FaceNormalVector)
	@test Dip2≈Dip
	@test maximum(AzDiff.(Azimuth2[Steep],Azimuth[Steep]))<1e-6

	#Neighbours
	(Dist,SortedTriangles,Noconnections)=CutAndDisplaceJulia.FindConnectedTriangles(Triangles,MidPoint)
	(SortedTrisCC,~)=CutAndDisplaceJulia.ConnectedConstraints(P1,P2,P3,MidPoint)
	Mine=[Set(filter(!=(0),SortedTriangles[i,[2,4,6]])) for i=1:n]
	Theirs=[Set(filter(!=(0),SortedTrisCC[i,:])) for i=1:n]
	@test Mine==Theirs
	@test Noconnections[:,2]==length.(Mine)
	@test all(i in Mine[j] for i=1:n for j in Mine[i])            #neighbours are mutual
	(F1,F2,F3)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	OnRim=F1.FreeFlg .| F2.FreeFlg .| F3.FreeFlg
	@test all(Noconnections[OnRim,2].<3) && all(Noconnections[.!OnRim,2].==3)
	@test all((Dist.>0).==(SortedTriangles[:,[2,4,6]].>0))

	if MakeFigures
		fig=CM.Figure(size=(1700,560))
		ax1=CM.Axis(fig[1,1],aspect=CM.DataAspect(),xlabel="x (east)",ylabel="y (north)",
			title="DipAndAzimuth: dip [°]\npenny bent into a trough z = 0.5x²\n(steepest at the east and west edges)")
		m1=FillTris!(ax1,P1,P2,P3,Dip,colormap=:viridis,colorrange=(0,45))
		DrawTriEdges!(ax1,P1,P2,P3,color=(:white,0.4),linewidth=0.4)
		CM.Colorbar(fig[1,2],m1,label="dip [°]")
		ax2=CM.Axis(fig[1,3],aspect=CM.DataAspect(),xlabel="x (east)",ylabel="y (north)",
			title="DipAndAzimuth: dip direction (azimuth)\neast half dips west (270°), west half east (90°)\narrows: direction of dip")
		m2=FillTris!(ax2,P1,P2,P3,Azimuth,colormap=:twilight,colorrange=(0,360))
		k=1:3:n
		CM.arrows2d!(ax2,MidPoint[k,1],MidPoint[k,2],0.08.*sind.(Azimuth[k]).*Steep[k],0.08.*cosd.(Azimuth[k]).*Steep[k],color=:black,shaftwidth=1,tipwidth=5,tiplength=5)
		CM.Colorbar(fig[1,4],m2,label="azimuth [°] (0 = north, 90 = east)",ticks=0:90:360)
		ax3=CM.Axis(fig[1,5],aspect=CM.DataAspect(),xlabel="x (east)",ylabel="y (north)",
			title="FindConnectedTriangles: neighbours per triangle\nsame as ConnectedConstraints: $(Mine==Theirs)\nrim triangles have fewer than 3")
		m3=FillTris!(ax3,P1,P2,P3,Float64.(Noconnections[:,2]),colormap=CM.cgrad([:red,:orange,:lightblue],3,categorical=true),colorrange=(0.5,3.5))
		DrawTriEdges!(ax3,P1,P2,P3,color=(:black,0.3),linewidth=0.4)
		CM.Colorbar(fig[1,6],m3,label="number of neighbours",ticks=1:3)
		CM.Label(fig[0,:],"Unit test 10: dip, azimuth and strike of each triangle, and each triangle's neighbours",fontsize=20)
		SaveTestFigure(fig,"ut_10_DipAzimuthAndNeighbours")
	end
end
