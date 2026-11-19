# Builds the meshes used by the elliptical / gradient stress-intensity tests and saves them in
# this folder, so the tests themselves do not need BuildCGAL (they run in Pkg.test and on GitHub).
# Only re-run this if you want to change those meshes. Needs BuildCGAL:
#	julia --project=<env with CutAndDisplaceJulia and BuildCGAL> MakeSIFTestMeshes.jl
# Writes:
#	EllipseMesh_a2_b1_1500Faces.off       ellipse, semi-axes 2 (x) and 1 (y), in the x-y plane
#	EllipseMesh_a1_b3_650Faces.off        ellipse, semi-axes 1 (x) and 3 (y), in the x-y plane
#	VertPennyMesh_r1_1500Faces.off        unit-radius penny in the x-z plane (vertical)
# The recipe is the one the tests used before: remesh to the target size with CGAL, clean and
# isoscelise the rim triangles, snap the rim back onto the exact ellipse/circle, repeat.
using CutAndDisplaceJulia, BuildCGAL

function RemeshCleanIsoscelise(Points,Triangles,target_edge_length)
	f=CutAndDisplaceJulia.OFFExport(Points,Triangles,size(Triangles,1),size(Points,1),"BeforePolygonRemshing")
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(BuildCGAL.PolygonRemeshingCGAL(f,target_edge_length))
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
	return Points,Triangles
end

#Snap the rim onto the ellipse x²/a² + y²/b² = 1 (via the unit circle)
function SnapRimToEllipse(Points,Triangles,a,b)
	Points[:,2]=Points[:,2]./a; Points[:,3]=Points[:,3]./b
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
	(P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles)=CutAndDisplaceJulia.PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,1,1,2)
	Points[:,2]=Points[:,2].*a; Points[:,3]=Points[:,3].*b
	return Points,Triangles
end

function MakeEllipse(a,b,NoTris)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.off"))
	Points[:,2]=Points[:,2]*a; Points[:,3]=Points[:,3]*b
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
	(target_edge_length,~)=CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,NoTris)
	for k=1:2
		(Points,Triangles)=RemeshCleanIsoscelise(Points,Triangles,target_edge_length)
		(Points,Triangles)=SnapRimToEllipse(Points,Triangles,a,b)
	end
	(Points,Triangles)=RemeshCleanIsoscelise(Points,Triangles,target_edge_length)
	return Points,Triangles
end

function MakeVerticalPenny(NoTris)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.off"))
	(Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,3],Points[:,4],0.0,0.0,cosd(90),sind(90))
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
	(target_edge_length,~)=CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,NoTris)
	Snap(Points,Triangles)=begin
		(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
		(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
		(P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles)=CutAndDisplaceJulia.PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,1.0,1,3)
		(Points,Triangles)
	end
	for k=1:3
		(Points,Triangles)=RemeshCleanIsoscelise(Points,Triangles,target_edge_length)
		(Points,Triangles)=Snap(Points,Triangles)
	end
	return Points,Triangles
end

cd(mktempdir()) do
	for (name,(Points,Triangles)) in (("EllipseMesh_a2_b1_1500Faces",MakeEllipse(2.0,1.0,1500)),
	                                  ("EllipseMesh_a1_b3_650Faces",MakeEllipse(1.0,3.0,650)),
	                                  ("VertPennyMesh_r1_1500Faces",MakeVerticalPenny(1500)))
		f=CutAndDisplaceJulia.OFFExport(Points,Triangles,size(Triangles,1),size(Points,1),name)
		cp(f,joinpath(@__DIR__,"$name.off");force=true)
		println("wrote $name.off: $(size(Triangles,1)) triangles")
	end
end
