# Graded remeshing of a crack surface: small triangles at the crack front (tipline), growing
# row by row to the normal size inside. Works for flat and curved (dish shaped, warped) cracks.
#
#   (Points,Triangles)=GradedRemesh(Points,Triangles,h,hf,BuildCGAL.GradedRemeshingCGAL;growth=1.3)
#
#   Points,Triangles  the crack mesh (as from OFFReader); the mesh border is taken as the front
#   h                 target edge length inside the crack
#   hf                edge length along the front (the tip elements), hf <= h
#   GradedRemesher    the CGAL step, BuildCGAL.GradedRemeshingCGAL (passed in so that
#                     CutAndDisplaceJulia itself does not depend on BuildCGAL)
#   growth            size ratio from one row of triangles to the next (1.3 gives triangles with
#                     all angles above ~24 degrees)
#   FrontCleanup      true: first a uniform remesh at h (the same CGAL remeshing as
#                     PolygonRemeshingCGAL: short front edges are collapsed, front vertices slide
#                     along the front), so the spline goes through evenly spaced points. Use this
#                     in the propagation loop: the mesh from AdvancingFrontCGAL has unevenly spaced
#                     front points (and small zigzags) that a spline through every point turns
#                     into spikes. The front points are then no longer honoured exactly (as in the
#                     uniform loop). Default false: every front point is kept.
#
# What it does:
#  1. Front: new front points every ~hf along a smooth closed curve (centripetal Catmull-Rom
#     spline) through the ORIGINAL front points. The original front points are kept, so the
#     points made by PropagateFracture are honoured exactly, and the tipline gets no corners at
#     the old vertices.
#  2. Tip triangles: for every front edge a point is placed at the apex of the equilateral
#     triangle on that edge, on the crack surface. Each front edge then has its own inner point,
#     which is the rim triangle shape the stress intensity formula in StressIntensity3D needs
#     (what CleanEdgeTris + IsosceliseEdgeTris produce on a uniform mesh).
#  3. CGAL (cgal_tools/src/graded_remeshing.cpp in BuildCGAL): isotropic remeshing with target
#     edge length min(h, hf + c*distance-from-front), keeping the front edges and the apex points
#     fixed, then angle smoothing. New vertices are put back onto the input surface. Where CGAL leaves
#     two front edges' triangles sharing an inner vertex, one gets a new inner vertex of its own
#     (GradedRemeshOwnRimVertices), so no front point is ever dropped by the rim clean-up.
#  4. Curved cracks: CGAL leaves new vertices on the input's flat facets. Every inner vertex is
#     moved onto a smooth surface (a quadric fitted by moving least squares) through the ORIGINAL
#     vertices. On a flat crack this changes nothing.
# Running CleanEdgeTris + IsosceliseEdgeTris afterwards leaves the mesh unchanged.
#
# Notes for curved cracks:
#  - "Distance from the front" in (3) is a straight-line distance. On a curved sheet that is
#    never longer than the distance along the sheet, so the mesh can only come out finer than
#    asked, never coarser.
#  - The surface used in (4) is fitted within ~2.5 input edge lengths of each point, so the
#    input mesh must resolve the curvature (as it does for the remeshed meshes in the loop).
function GradedRemesh(Points,Triangles,h,hf,GradedRemesher;growth=1.3,FrontCleanup=false)

	hf<=h || error("GradedRemesh: the front size hf ($hf) must not be larger than h ($h)")
	X0=Points[:,2:4]                         #the original surface
	N0=GradedRemeshVertexNormals(X0,Triangles)
	R=2.5*GradedRemeshMeanEdgeLength(X0,Triangles)

	#0. Optional clean-up of the front: uniform remesh at h (as PolygonRemeshingCGAL)
	if FrontCleanup
		(Points,Triangles)=GradedRemeshRunCGAL(GradedRemesher,Points,Triangles,h,h,Int[];
		                                       ProtectFront=false,smooth=0)
	end

	#1. Front points every ~hf on the smooth front through the original front points
	(X,T,Loops)=GradedRemeshResampleFront(Points[:,2:4],Triangles,hf)

	#2. One inner point per front edge at the apex of its equilateral triangle, on the surface.
	#   Skipped where it would come within 0.6hf of a point already placed (very tight bends).
	Placed=[X[i] for L in Loops for i in L]
	Apex=Vector{Vector{Float64}}()
	for L in Loops, k in eachindex(L)
		a=X[L[k]]; b=X[L[mod1(k+1,length(L))]]; e=b-a
		mid=(a+b)/2
		m=normalize(cross(GradedRemeshSurfaceNormal(mid,X0,N0,R),e))  #in the surface, ⟂ to the edge
		#make sure m points into the crack: towards the third vertex of the edge's triangle
		t=T[findfirst(t->(L[k] in t)&&(L[mod1(k+1,length(L))] in t),T)]
		c=X[first(setdiff(t,[L[k],L[mod1(k+1,length(L))]]))]
		dot(m,c-mid)<0 && (m=-m)
		p=GradedRemeshMLSProject(mid+sqrt(3)/2*norm(e)*m,X0,N0,R)
		minimum(norm(p-q) for q in Placed)<0.6hf && continue
		push!(Apex,p); push!(Placed,p)
	end
	Fixed=[GradedRemeshInsertPoint!(X,T,p) for p in Apex]
	#The front points are held fixed too: CGAL protects the front edges from being split or collapsed,
	#but without this it can still collapse an inner edge from a front point onto a fixed apex point,
	#which re-routes the front through that apex (moving the front by up to ~0.7hf).
	append!(Fixed,[i for L in Loops for i in L])

	#3. CGAL graded remeshing with the front and the apex points held fixed
	Pin=[collect(1:length(X)) permutedims(reduce(hcat,X))]
	Tin=permutedims(reduce(hcat,T))
	(P,Tout)=GradedRemeshRunCGAL(GradedRemesher,Pin,Tin,h,hf,Fixed;growth=growth,ProtectFront=true,
	                             band=sqrt(3)/2*hf)

	#3b. CGAL keeps the front edges and the apex points, but not the edges between them, so a few front
	#    edges' triangles can end up sharing an inner vertex. Give each its own (the front is not touched).
	(P,Tout,NewApex)=GradedRemeshOwnRimVertices(P,Tout)
	append!(Apex,NewApex)

	#4. Every vertex that is not on the front and not an apex point goes onto the smooth surface
	Keep=falses(size(P,1))
	for L in GradedRemeshBorderLoops(Tout), i in L; Keep[i]=true; end
	for p in Apex   #the apex points come back unmoved (to the 13 digits of the .off files); new ones are in place
		d=vec(sum((P[:,2:4].-p').^2,dims=2)); j=argmin(d)
		d[j]<1e-16 && (Keep[j]=true)
	end
	for i in findall(.!Keep)
		P[i,2:4]=GradedRemeshMLSProject(P[i,2:4],X0,N0,R)
	end
	return P,Tout
end

# Every front edge gets a triangle with an inner vertex of its own (what StressIntensity3D needs, and
# what CleanEdgeTris + IsosceliseEdgeTris would otherwise make by merging the triangles and dropping
# the front point between them). Two cases are repaired, without changing any front point or edge:
#  - a triangle with two front edges (a sharp tip of the front): each front edge gets a new vertex;
#  - front edges whose triangles share an inner vertex: all but the most regular get a new vertex.
# A new vertex goes on the perpendicular bisector of its front edge, inside the triangle, at the height
# (up to the equilateral one) that gives the new triangles the largest smallest angle.
# Returns the points, triangles and the new vertices.
function GradedRemeshOwnRimVertices(P,T)
	X=[Vector(r) for r in eachrow(P[:,2:4])]; T=[Vector(r) for r in eachrow(T)]
	cnt=Dict{Tuple{Int,Int},Int}()
	for t in T, k=1:3; e=minmax(t[k],t[mod1(k+1,3)]); cnt[e]=get(cnt,e,0)+1; end
	Front(x,y)=cnt[minmax(x,y)]==1
	New=Vector{Vector{Float64}}()
	Nrm(p,q,r)=cross(q-p,r-p)
	MinAngle(p,q,r)=minimum(acosd(clamp(dot(u,v)/(norm(u)*norm(v)),-1,1)) for (u,v) in ((q-p,r-p),(p-q,r-q),(p-r,q-r)))
	#candidate points for front edge a→b inside triangle (a,b,c): on the bisector, up to the equilateral
	#height and up to 0.9 of the way to the triangle's far side
	function Candidates(a,b,c)
		e=X[b]-X[a]; L=norm(e); u=e/L; mid=(X[a]+X[b])/2
		w=X[c]-mid; m=normalize(w-dot(w,u)*u)
		cx=dot(w,u); cy=dot(w,m)
		sb=cx<L/2 ? (L/2)/(L/2-cx) : Inf; sa=cx>-L/2 ? (L/2)/(cx+L/2) : Inf
		Hmax=min(sqrt(3)/2*L,0.9*cy*min(sb<=1 ? sb : Inf, sa<=1 ? sa : Inf))
		return [mid+f*Hmax*m for f in range(0.15,1.0,length=12)]
	end
	#best layout: the candidate(s) whose triangles keep the orientation and have the largest smallest angle
	function Best(cands,tris,n0)
		best=-Inf; arg=nothing
		for q in cands
			pos(v)=v isa Int ? X[v] : q[v.k]
			all(dot(Nrm(pos.(t)...),n0)>0 for t in tris) || continue
			a=minimum(MinAngle(pos.(t)...) for t in tris)
			a>best && (best=a; arg=q)
		end
		return arg
	end
	Add!(q)=(push!(X,q); push!(New,q); length(X))
	#1. tips of the front: triangle a→b→d with front edges a→b and b→d
	for i in 1:length(T)
		t=T[i]; k=findfirst(k->Front(t[k],t[mod1(k+1,3)]) && Front(t[mod1(k+1,3)],t[mod1(k+2,3)]),1:3)
		(k===nothing || Front(t[mod1(k+2,3)],t[k])) && continue      #no tip, or a lone triangle
		a,b,d=t[k],t[mod1(k+1,3)],t[mod1(k+2,3)]
		q1,q2=(k=1,),(k=2,)                                            #placeholders for the new vertices
		Tri=[(a,b,q1),(b,d,q2),(b,q2,q1),(q1,q2,d),(a,q1,d)]
		q=Best([(c1,c2) for c1 in Candidates(a,b,d) for c2 in Candidates(b,d,a)],Tri,Nrm(X[a],X[b],X[d]))
		q===nothing && continue
		j1=Add!(q[1]); j2=Add!(q[2]); Id(v)=v===q1 ? j1 : v===q2 ? j2 : v
		T[i]=[Id(v) for v in Tri[1]]
		for tr in Tri[2:end]; push!(T,[Id(v) for v in tr]); end
	end
	#2. front edges whose triangles share an inner vertex
	Rim=Dict{Int,Vector{Tuple{Int,Int,Int}}}()          #inner vertex => (triangle, a, b), a→b as in the triangle
	for (i,t) in enumerate(T), k=1:3
		a,b,c=t[k],t[mod1(k+1,3)],t[mod1(k+2,3)]
		haskey(cnt,minmax(a,b)) && Front(a,b) && push!(get!(Rim,c,Tuple{Int,Int,Int}[]),(i,a,b))
	end
	for (c,R) in Rim
		length(R)>1 || continue
		Q(r)=(L=norm(X[r[3]]-X[r[2]]); abs(norm(X[c]-X[r[2]])-L)+abs(norm(X[c]-X[r[3]])-L))
		for (i,a,b) in R[sortperm(Q.(R))][2:end]                     #keep the most regular, split the others
			q1=(k=1,)
			q=Best([(p,) for p in Candidates(a,b,c)],[(a,b,q1),(b,c,q1),(c,a,q1)],Nrm(X[a],X[b],X[c]))
			q===nothing && continue
			j=Add!(q[1])
			T[i]=[a,b,j]; push!(T,[b,c,j]); push!(T,[c,a,j])
		end
	end
	Pout=[collect(1:length(X)) permutedims(reduce(hcat,X))]
	return Pout,permutedims(reduce(hcat,T)),New
end

# Run the CGAL step in a temporary folder: mesh in, mesh out
function GradedRemeshRunCGAL(GradedRemesher,Points,Triangles,h,hf,Fixed;growth=1.3,ProtectFront=true,
                             band=0.0,smooth=10)
	Dir=mktempdir()
	(P,T)=cd(Dir) do
		f=OFFExport(Points,Triangles,size(Triangles,1),size(Points,1),"GradedRemeshIn")
		fx=""
		if !isempty(Fixed)
			fx=joinpath(Dir,"GradedRemeshFixed.txt")
			open(io->foreach(i->println(io,i),Fixed),fx,"w")
		end
		OFFReader(GradedRemesher(f,h,hf;growth=growth,ProtectFront=ProtectFront,FixedVerticesTxt=fx,
		                         band=band,smooth=smooth))
	end
	rm(Dir;recursive=true,force=true)
	return P,T
end

# Front loops: ordered vertex numbers along the mesh border, one vector per loop (no repeat of
# the first point). Each front edge is followed in the direction it has in its own triangle, so
# every loop runs the same way round as the triangles, and a "pinch" vertex where two parts of
# the front touch (four front edges) is handled: each pass through it is a separate visit.
function GradedRemeshBorderLoops(Triangles)
	cnt=Dict{Tuple{Int,Int},Int}()
	for t in eachrow(Triangles), (a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
		k=minmax(a,b); cnt[k]=get(cnt,k,0)+1
	end
	Out=Dict{Int,Vector{Int}}()                     #front edges a→b, as in their triangle
	for t in eachrow(Triangles), (a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
		cnt[minmax(a,b)]==1 && push!(get!(Out,a,Int[]),b)
	end
	Loops=Vector{Vector{Int}}()
	for s in sort(collect(keys(Out)))
		while !isempty(Out[s])
			L=[s]; v=s
			while true
				nxt=get(Out,v,Int[])
				isempty(nxt) && error("GradedRemesh: the mesh border is not a closed loop (are the triangles consistently oriented?)")
				w=popfirst!(nxt)
				w==s && break
				push!(L,w); v=w
			end
			push!(Loops,L)
		end
	end
	return Loops
end

# Points every `per`-th of the way along a closed centripetal Catmull-Rom spline through the
# rows of X (in order). The spline goes through every point of X and does not overshoot or
# loop at uneven spacing. Segment k (X[k]→X[k+1]) is rows (k-1)*per+1 : k*per+1 (wrapping).
function GradedRemeshClosedSpline(X;per=40)
	n=size(X,1); S=Matrix{Float64}(undef,n*per,3); k=0
	P(i)=X[mod1(i,n),:]
	for i=1:n
		p0,p1,p2,p3=P(i-1),P(i),P(i+1),P(i+2)
		t0=0.0; t1=t0+norm(p1-p0)^0.5; t2=t1+norm(p2-p1)^0.5; t3=t2+norm(p3-p2)^0.5
		for t in range(t1,t2,length=per+1)[1:end-1]
			A1=(t1-t)/(t1-t0)*p0+(t-t0)/(t1-t0)*p1
			A2=(t2-t)/(t2-t1)*p1+(t-t1)/(t2-t1)*p2
			A3=(t3-t)/(t3-t2)*p2+(t-t2)/(t3-t2)*p3
			B1=(t2-t)/(t2-t0)*A1+(t-t0)/(t2-t0)*A2
			B2=(t3-t)/(t3-t1)*A2+(t-t1)/(t3-t1)*A3
			k+=1; S[k,:]=(t2-t)/(t2-t1)*B1+(t-t1)/(t2-t1)*B2
		end
	end
	return S
end

# New front points every ~hf along the spline. Each original front edge a→b is cut into
# n=round(arc length/hf) equal pieces, and the triangle that owned a→b becomes a fan onto the new
# points (CGAL tidies the fans up). Returns points and triangles as vectors, and the front loops.
function GradedRemeshResampleFront(X0,Triangles,hf;per=40)
	X=[Vector(r) for r in eachrow(X0)]; T=[Vector(r) for r in eachrow(Triangles)]
	Loops=GradedRemeshBorderLoops(Triangles); NewLoops=Vector{Vector{Int}}()
	for L in Loops
		S=GradedRemeshClosedSpline(X0[L,:];per=per); m=length(L)
		NewL=Int[]
		for k=1:m
			a=L[k]; b=L[mod1(k+1,m)]
			push!(NewL,a)
			seg=S[[((k-1)*per+1):(k*per); mod1(k*per+1,size(S,1))],:]
			s=[0.0; cumsum(norm.(eachrow(diff(seg,dims=1))))]
			n=max(1,round(Int,s[end]/hf))
			n==1 && continue
			ids=Int[]
			for v in range(0,s[end],length=n+1)[2:end-1]
				j=clamp(searchsortedlast(s,v),1,length(s)-1); t=(v-s[j])/(s[j+1]-s[j])
				push!(X,seg[j,:]+t*(seg[j+1,:]-seg[j,:])); push!(ids,length(X))
			end
			append!(NewL,ids)
			#the triangle with edge a→b becomes a fan a,ids...,b onto its third vertex
			ti=findfirst(t->(a in t)&&(b in t),T); t=T[ti]
			c=first(setdiff(t,[a,b]))
			chain=[a;ids;b]
			i=findfirst(==(a),t)
			deleteat!(T,ti)
			#keep the triangle's orientation: its edge runs a→b or b→a
			append!(T, t[mod1(i+1,3)]==b ? [[chain[q],chain[q+1],c] for q=1:length(chain)-1] :
			                               [[chain[q+1],chain[q],c] for q=1:length(chain)-1])
		end
		push!(NewLoops,NewL)
	end
	return X,T,NewLoops
end

# Split the triangle that p lies over into three, p becoming a vertex. Keeps the orientation.
function GradedRemeshInsertPoint!(X,T,p)
	best=Inf; bi=0
	for (k,t) in enumerate(T)
		a,b,c=X[t[1]],X[t[2]],X[t[3]]
		n=cross(b-a,c-a); A2=norm(n); A2<eps() && continue; n/=A2
		q=p-dot(p-a,n)*n
		u=dot(cross(c-b,q-b),n)/A2; v=dot(cross(a-c,q-c),n)/A2; w=1-u-v
		score=abs(dot(p-a,n))+max(0,-min(u,v,w))*sqrt(A2)   #off the plane + outside the triangle
		score<best && (best=score; bi=k)
	end
	push!(X,p); j=length(X); t=T[bi]
	deleteat!(T,bi); append!(T,[[t[1],t[2],j],[t[2],t[3],j],[t[3],t[1],j]])
	return j
end

# Area-weighted vertex normals (follow the triangles' orientation)
function GradedRemeshVertexNormals(X,Triangles)
	N=zeros(size(X))
	for t in eachrow(Triangles)
		fn=cross(X[t[2],:]-X[t[1],:],X[t[3],:]-X[t[1],:])
		for j in t; N[j,:]+=fn; end
	end
	return N./max.(sqrt.(sum(N.^2,dims=2)),eps())
end

function GradedRemeshMeanEdgeLength(X,Triangles)
	mean(norm(X[t[i],:]-X[t[mod1(i+1,3)],:]) for t in eachrow(Triangles) for i=1:3)
end

# Smooth surface normal at p: Gaussian-weighted average of the vertex normals within ~R
function GradedRemeshSurfaceNormal(p,X,N,R)
	d=vec(sqrt.(sum((X.-p').^2,dims=2)))
	w=exp.(-(d./(0.5R)).^2)
	return normalize(vec(sum(N.*w,dims=1)))
end

# Moving least squares: fit a quadric height function z(u,v) to the vertices X near p (in the
# local frame of the surface normal), weighted by distance, and move p onto it.
# A flat set of points gives a plane, so p does not move on a flat crack.
function GradedRemeshMLSProject(p,X,N,R)
	d=vec(sqrt.(sum((X.-p').^2,dims=2)))
	I=findall(d.<R)
	length(I)<8 && (I=partialsortperm(d,1:min(12,length(d))))
	w=exp.(-(d[I]./(0.5R)).^2)
	n=normalize(vec(sum(N[I,:].*w,dims=1)))
	e1=normalize(abs(n[1])<0.9 ? cross(n,[1.0,0,0]) : cross(n,[0,1.0,0])); e2=cross(n,e1)
	Q=X[I,:].-p'
	u=Q*e1; v=Q*e2; z=Q*n
	A=[ones(length(I)) u v u.^2 u.*v v.^2]
	c=(A'*(w.*A))\(A'*(w.*z))
	return p+c[1]*n
end
