#test remeshing

meshpath="C:\\Users\\timmm\\AppData\\Local\\Julia-1.1.0\\bin\\remeshed.off"
#Reload	
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(meshpath)
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);




@enter CutAndDisplaceJulia.CleanAndIsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)