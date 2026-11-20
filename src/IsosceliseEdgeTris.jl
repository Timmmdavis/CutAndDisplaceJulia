"""
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)

Makes every rim triangle (a triangle with an edge on the crack front) isosceles about that front edge, with an inner
point of its own: what StressIntensity3D needs (KI is computed from the distance between the triangle's centre and the
midpoint of its front edge, which must be at right angles to the edge). The mesh is only changed by moving inner points
and by splitting triangles:
 - a tip of the front (a triangle with two front edges) is split into five, so each front edge gets its own inner point;
 - a rim triangle whose inner point is used by no other rim triangle and is not on the front: that point is moved onto
   the perpendicular bisector of the front edge: straight across if that leaves the triangles around it no worse,
   else to the position along it (up to the equilateral height) where they are best shaped. A move is only taken if
   no triangle around the point folds or tips over (turns >30°, or leans off the surface more than it already did)
   and the triangles are not clearly (>1°) worse shaped than they were and than a split would make them;
 - otherwise (the inner point is on the front, is shared by several rim triangles, or cannot be moved like that) the
   rim triangle is split by a new point on the bisector, inside it (the most regular of a shared group keeps the
   shared point). If the rim triangle is lopsided, so a split would leave slivers, one of its inner edges is flipped
   first when that gives clearly better triangles (the point across it then becomes the rim triangle's inner point).
Works the same on meshes with holes (an inner front) and on several separate pieces: only the local topology is used.
No triangle is ever deleted and no point on the front is ever moved, so the crack front is unchanged. A rim triangle
that cannot be fixed without folding the mesh is left as it is (with a warning). Only P1, P2 and P3 are used as input;
all the outputs describe the same, new mesh.
"""
function IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)

	(Tri,Pts)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
	X=[Vector{Float64}(r) for r in eachrow(Pts[:,2:4])]; T=[Vector{Int}(r) for r in eachrow(Tri)]
	cnt=Dict{Tuple{Int,Int},Int}()
	for t in T, k=1:3; e=minmax(t[k],t[mod1(k+1,3)]); cnt[e]=get(cnt,e,0)+1; end
	Front(x,y)=get(cnt,minmax(x,y),0)==1                  #the front edges never change, so this stays valid
	OnFront=Set(v for (e,c) in cnt if c==1 for v in e)

	Nrm(p,q,r)=cross(q-p,r-p)
	MinAngle(p,q,r)=minimum(acosd(clamp(dot(u,v)/(norm(u)*norm(v)),-1,1)) for (u,v) in ((q-p,r-p),(p-q,r-q),(p-r,q-r)))
	Isosceles(a,b,c)=abs(norm(X[c]-X[a])-norm(X[c]-X[b]))<=1e-6*norm(X[b]-X[a])
	#points on the bisector of front edge a→b, towards c, from 15% to 100% of Hmax: the equilateral height, or 0.9
	#of the way to the triangle's far side (frac) if that is nearer
	function Bisector(a,b,c;frac=0.9,inside=true)
		e=X[b]-X[a]; L=norm(e); u=e/L; mid=(X[a]+X[b])/2
		w=X[c]-mid; m=normalize(w-dot(w,u)*u)
		cx=dot(w,u); cy=dot(w,m)
		Hmax=sqrt(3)/2*L
		if inside
			sb=cx<L/2 ? (L/2)/(L/2-cx) : Inf; sa=cx>-L/2 ? (L/2)/(cx+L/2) : Inf
			Hmax=min(Hmax,frac*cy*min(sb<=1 ? sb : Inf, sa<=1 ? sa : Inf))
		end
		return [mid+f*Hmax*m for f in range(0.15,1.0,length=12)]
	end
	#best of the candidate layouts: all new triangles keep their orientation (n0: one normal, or one per triangle);
	#largest smallest angle. Returns the layout (nothing if every one folds), or (layout, its smallest angle).
	function Best(cands,tris,n0;withangle=false)
		best=-Inf; arg=nothing
		for q in cands
			pos(v)=v isa Int ? X[v] : q[v.k]
			all(dot(Nrm(pos.(t)...),n0 isa Vector{<:Vector} ? n0[k] : n0)>0 for (k,t) in enumerate(tris)) || continue
			ang=minimum(MinAngle(pos.(t)...) for t in tris)
			ang>best && (best=ang; arg=q)
		end
		return withangle ? (arg,best) : arg
	end
	Add!(q)=(push!(X,q); length(X))
	Failed=0

	#1. Tips of the front: triangle a→b→d with front edges a→b and b→d gets two new inner points, q1 and q2, on the
	#   bisectors of a→b and b→d, and becomes five triangles: each front edge gets its own isosceles rim triangle
	#   (a-b-q1 and b-d-q2). a-b and b-d are the front; b is the tip:
	#
	#                          b                                                   b
	#                        /   \                                               // \\
	#                       /     \                                             //   \\
	#                     /         \                                         /  /   \  \
	#                   /             \                                     /   /     \   \
	#                  /               \                                   /   /       \   \
	#                /                   \                               /     /       \     \
	#              /                       \            ==>            /      /         \      \
	#            /                           \                       /     _-q1‾‾‾‾‾‾‾‾‾‾q2_     \
	#           /                             \                     /   _-‾        ‾--__    ‾-_   \
	#         /                                 \                 /  _-‾                ‾--__  ‾-_  \
	#       /                                     \             / _-‾                        ‾--__‾-_ \
	#      /                                       \           /-‾                                ‾‾-__\
	#    a-------------------------------------------d       a-------------------------------------------d
	#
	for i in 1:length(T)
		t=T[i]; k=findfirst(k->Front(t[k],t[mod1(k+1,3)]) && Front(t[mod1(k+1,3)],t[mod1(k+2,3)]),1:3)
		(k===nothing || Front(t[mod1(k+2,3)],t[k])) && continue      #no tip, or a lone triangle
		a,b,d=t[k],t[mod1(k+1,3)],t[mod1(k+2,3)]
		q1,q2=(k=1,),(k=2,)
		Tri5=[(a,b,q1),(b,d,q2),(b,q2,q1),(q1,q2,d),(a,q1,d)]
		q=Best([(c1,c2) for c1 in Bisector(a,b,d;frac=0.6) for c2 in Bisector(b,d,a;frac=0.6)],Tri5,Nrm(X[a],X[b],X[d]))
		q===nothing && (Failed+=1; continue)
		j1=Add!(q[1]); j2=Add!(q[2]); Id(v)=v===q1 ? j1 : v===q2 ? j2 : v
		T[i]=[Id(v) for v in Tri5[1]]
		for tr in Tri5[2:end]; push!(T,[Id(v) for v in tr]); end
	end

	#2. The other rim triangles, grouped by inner point. The front edge a→b is drawn ===, the mesh carries on below;
	#   α marks the two equal angles of the isosceles rim triangle (60° for an equilateral one).
	#   Each rim triangle a-b-c either has its inner point c moved onto the bisector of a-b (c → c'), taking the
	#   triangles around c with it:
	#
	#    a===================================b               a===================================b
	#    |\                              _-‾ |               | \α                             α/ |
	#    | \                          _-‾    |               |   \                           /   |
	#     | \                      _-‾      |                 |    \                       /    |
	#     |  \                   -‾         |                 |      \                   /      |
	#     |    \              _-‾           |                 |        \               /        |
	#      |    \          _-‾             |                   |         \           /         |
	#      |     \      _-‾                |                   |           \       /           |
	#      |      \  _-‾                   |        ==>        |             \   /             |
	#      |       c‾                      |                   |               c'              |
	#      |      /  ‾-_                   |                   |             /   \             |
	#       |    /      ‾-_               |                     |          /       \          |
	#       |   /          ‾-_            |                     |        /           \        |
	#       |   /             ‾‾-_        |                     |      /               \      |
	#        | /                  ‾-_    |                       |   /                   \   |
	#        |/                      ‾-_ |                       | /                       \ |
	#        .--------------------------‾.                       .---------------------------.
	#
	#   or (c on the front, shared with another rim triangle, or not movable without folding, tilting or worse
	#   triangles) it is split by a new point q on the bisector, inside it:
	#
	#    a===================================b               a===================================b
	#    |                                -‾                 | ‾‾α_                         _α-‾
	#     |                            _-‾                    |    ‾‾-_                 _-‾_-‾
	#     |                         _-‾                       |        ‾‾-_         _-‾‾_-‾
	#     |                       -‾                          |            ‾‾-_q_-‾‾  -‾
	#      |                   _-‾                             |             /     _-‾
	#      |                _-‾                     ==>        |           /    _-‾
	#      |              -‾                                   |          /   -‾
	#       |          _-‾                                      |       /  _-‾
	#       |       _-‾                                         |     / _-‾
	#       |     -‾                                            |    /-‾
	#        | _-‾                                               | _-‾
	#        c‾                                                  c‾
	#
	Rim=Dict{Int,Vector{Tuple{Int,Int,Int}}}()               #inner point => (triangle, a, b), a→b as in the triangle
	for (i,t) in enumerate(T), k=1:3
		a,b,c=t[k],t[mod1(k+1,3)],t[mod1(k+2,3)]
		Front(a,b) && !Front(b,c) && !Front(c,a) && push!(get!(Rim,c,Tuple{Int,Int,Int}[]),(i,a,b))
	end
	Star=Dict{Int,Vector{Int}}()                              #triangles around each point
	for (i,t) in enumerate(T), v in t; push!(get!(Star,v,Int[]),i); end
	Split=Tuple{Int,Int,Int,Int}[]                            #(triangle, a, b, c) to split, after all the moves
	for (c,R) in Rim
		Regularity(r)=(L=norm(X[r[3]]-X[r[2]]); abs(norm(X[c]-X[r[2]])-L)+abs(norm(X[c]-X[r[3]])-L))
		R=R[sortperm(Regularity.(R))]
		keep=0
		if !(c in OnFront)
			(i,a,b)=R[1]
			if Isosceles(a,b,c)
				keep=1
			else
				#move c onto the bisector of a→b (straight across, or along it): every triangle around c must keep its
				#orientation. Taken only if it leaves better shaped triangles than splitting the rim triangle would.
				around=Star[c]; n0=[normalize(Nrm(X[T[s][1]],X[T[s][2]],X[T[s][3]])) for s in around]
				nbar=normalize(sum(Nrm(X[T[s][1]],X[T[s][2]],X[T[s][3]]) for s in around))  #the surface's normal at c
				lim=[min(0.87,dot(n,nbar))-0.02 for n in n0]
				#a triangle around c may turn by <30°, or lean no more than 30° (or than it already did) off the surface's normal
				Upright(n,k)=dot(n,n0[k])>0 && (dot(n,n0[k])>0.87 || dot(n,nbar)>lim[k])
				e=X[b]-X[a]; u=e/norm(e); mid=(X[a]+X[b])/2
				function MoveAngle(q)       #smallest angle around c with c at q; -Inf if a triangle would fold or tilt
					P(v)=v==c ? q : X[v]
					all(Upright(normalize(Nrm(P.(T[s])...)),k) for (k,s) in enumerate(around)) || return -Inf
					return minimum(MinAngle(P.(T[s])...) for s in around)
				end
				q1=(k=1,)
				qs=Best([(p,) for p in Bisector(a,b,c)],[(a,b,q1),(b,c,q1),(c,a,q1)],Nrm(X[a],X[b],X[c]))
				splitang=qs===nothing ? -Inf : minimum(MinAngle(t...) for t in ((X[a],X[b],qs[1]),(X[b],X[c],qs[1]),(X[c],X[a],qs[1])))
				origang=minimum(MinAngle(X[T[s][1]],X[T[s][2]],X[T[s][3]]) for s in around)
				Good(ang)=ang>=min(splitang,origang)-1.0      #not clearly (>1°) worse than before and than a split
				#the smallest move (straight across onto the bisector) if it leaves the triangles no worse than they were,
				#else the best position along the bisector unless that is hardly (<2°) better
				q0=X[c]-dot(X[c]-mid,u)*u
				cands=Bisector(a,b,c;inside=false); angs=MoveAngle.(cands); k=argmax(angs)
				ang0=MoveAngle(q0)
				ang0>=min(origang-1.0,angs[k]-2.0) && (cands[k]=q0; angs[k]=ang0)
				Good(angs[k]) && (X[c]=cands[k]; keep=1)
			end
		end
		for r in R[keep+1:end]; push!(Split,(r[1],r[2],r[3],c)); end
	end
	#Each is split by a point on its bisector, unless flipping one of its inner edges first gives clearly (>2°) better
	#shaped triangles: a lopsided rim triangle (inner point far off the bisector) makes slivers when split, while the
	#point across an inner edge (if it is free: not on the front nor another rim triangle's) may sit much better.
	#
	#   A lopsided rim triangle (c far off the bisector) only splits into slivers, so the inner edge c-a (or b-c) is
	#   flipped first when the point d across it is free: a-b-c + a-c-d become b-c-d + the new rim triangle a-b-d,
	#   which is then split (or kept, if already isosceles):
	#
	#    a=============================b         a=============================b         a=============================b
	#     \  ‾‾--__                    \          \                          / \          \‾α-_                   _α‾/ \
	#      \       ‾‾--__               \          \                       /    \          \   ‾--_           _--‾ /    \
	#       \            ‾‾‾--__        \           \                     /     \           \      ‾-__   __-‾    /     \
	#        \                  ‾‾--__   \           \                  /        \           \         ‾q‾      /        \
	#         \                       ‾‾_c            \               /         _c            \        /      /         _c
	#         \                     _-‾‾              \             /       _-‾‾              \        /    /       _-‾‾
	#          \                __-‾       ==>         \          /     __-‾       ==>         \      /   /     __-‾
	#           \            _-‾                        \       /    _-‾                        \    /  /    _-‾
	#            \       _-‾‾                            \     / _-‾‾                            \  /  / _-‾‾
	#             \  __-‾                                 \  __-‾                                 \ /__-‾
	#              d‾                                      d‾                                      d‾
	#
	Owner=Dict{Tuple{Int,Int},Int}()                          #directed edge => triangle
	for (i,t) in enumerate(T), k=1:3; Owner[(t[k],t[mod1(k+1,3)])]=i; end
	Inner=Set(keys(Rim)); Touched=Set{Int}()
	for (i,a,b,c) in Split
		q1=(k=1,); n0=Nrm(X[a],X[b],X[c])
		(q,ang)=Best([(p,) for p in Bisector(a,b,c)],[(a,b,q1),(b,c,q1),(c,a,q1)],n0;withangle=true)
		plan=(q,[[a,b,q1],[b,c,q1],[c,a,q1]],0,0)            #(new point, new triangles, neighbour, its free point)
		for (x,y) in ((c,a),(b,c))
			j=get(Owner,(y,x),0); (j==0 || j in Touched) && continue
			t=T[j]; d=first(setdiff(t,(x,y)))
			(length(unique(t))<3 || d in OnFront || d in Inner || any(Front(t[k],t[mod1(k+1,3)]) for k=1:3)) && continue
			other=x==c ? [b,c,d] : [a,d,c]
			nj=Nrm(X[t[1]],X[t[2]],X[t[3]]); ns=[n0,n0,n0,nj]
			Isosceles(a,b,d) ? (qf,af)=Best([()],[(a,b,d),Tuple(other)],[n0,nj];withangle=true) :
			                   (qf,af)=Best([(p,) for p in Bisector(a,b,d)],[(a,b,q1),(b,d,q1),(d,a,q1),Tuple(other)],ns;withangle=true)
			qf===nothing && continue
			#the crease between the two flipped triangles may not be sharper than the one it replaces (or 30°)
			crease=dot(normalize(Nrm(X[a],X[b],X[d])),normalize(Nrm(X[other[1]],X[other[2]],X[other[3]])))
			crease<min(0.87,dot(normalize(n0),normalize(nj)))-0.02 && continue
			af>max(ang,-1e9)+2.0 && (ang=af; plan=(qf,Isosceles(a,b,d) ? [[a,b,d],other] : [[a,b,q1],[b,d,q1],[d,a,q1],other],j,d))
		end
		(q,tris,j,d)=plan
		q===nothing && (Failed+=1; continue)
		k=isempty(q) ? 0 : Add!(q[1])
		new=[[v===q1 ? k : v for v in t] for t in tris]
		T[i]=new[1]; push!(Touched,i)
		if j>0
			T[j]=new[end]; push!(Touched,j); push!(Inner,d); new=new[1:end-1]
		end
		for t in new[2:end]; push!(T,t); push!(Touched,length(T)); end
	end
	Failed>0 && @warn "IsosceliseEdgeTris: $Failed rim triangle(s) could not be made isosceles without folding the mesh; left as they were"

	Points=[collect(1:length(X)) permutedims(reduce(hcat,X))]
	Triangles=permutedims(reduce(hcat,T))
	(P1,P2,P3)=CreateP1P2P3(Triangles,Points)
	(FaceNormalVector,MidPoint)=CreateFaceNormalAndMidPoint(Points,Triangles)
	return P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector
end


#Point j of a list of edge points (LeadingPoints, TrailingPoints or InnerPoints from GetSortedEdgesOfMeshList),
#and the triangle it is in. Used by MovingAverageOfStressIntensity, PutOnEdgeOfCirc and CheckIfInsidePreviousBoundary.
function GrabPointNew(PointsIdxList,P1,P2,P3,j)
#Extract the points on the current bit of the edge
Point=[0. 0. 0.]
Indx=0
for k=1:3
    Idx=PointsIdxList[j,k]
    if Idx==0
        continue
    elseif k==1
        Point=P1[Idx,:]
        Indx=Idx
    elseif k==2
        Point=P2[Idx,:]
        Indx=Idx
    elseif k==3
        Point=P3[Idx,:]
        Indx=Idx
    end
end

return Point,Indx
end
