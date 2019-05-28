function testp2(UniqueEdges,LeadingPoint,TrailingPoint,InnerPoint,LeadingPoints,TrailingPoints,InnerPoints,scene,Idx,P1,P2,P3)

#For each edge loop
for i=1:length(UniqueEdges)
	#Run through list of this edge (Unique edges is a unit range)
	for j=UniqueEdges[i] 

		LeadingPointOld  =LeadingPoint
		TrailingPointOld =TrailingPoint
		InnerPointOld    =InnerPoint
		IdxOld=Idx

		#Extract the points on the current bit of the edge
		(LeadingPoint,~) =GrabPoint(LeadingPoints,P1,P2,P3,j)
		(TrailingPoint,~)=GrabPoint(TrailingPoints,P1,P2,P3,j)
		(InnerPoint,~) =GrabPoint(InnerPoints,P1,P2,P3,j)
		
		@info LeadingPoint
		scatter!(scene,[LeadingPoint[1] LeadingPoint[2] LeadingPoint[3]],markersize = 70)

	end
end

end