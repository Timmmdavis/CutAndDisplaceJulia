function CreateTrianglesPointsFromP1P2P3(P1,P2,P3)

Points=zeros(Int(length(P1)),3)
Points[1:3:end,:]=copy(P1);
Points[2:3:end,:]=copy(P2);
Points[3:3:end,:]=copy(P3);
n=length(Points)
Points=[1:n/3 Points]
Triangles=fill(0,Int(n/9),3)
Triangles[:,1]=1:3:n/3;
Triangles[:,2]=2:3:n/3;
Triangles[:,3]=3:3:n/3;


return Points,Triangles
end 