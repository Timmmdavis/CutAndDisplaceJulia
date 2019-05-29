function FindIntersectionOf3DVectors(r1,r2,e1,e2)
#r1 - start location of the first vector [x y z]
#r2 - start location of the 2nd vector   [x y z]
#e1 - direction (cosines) of the first vector [ax ay az]
#e2 - direction (cosines) of the 2nd vector   [ax ay az]

#https://stackoverflow.com/questions/10551555/need-an-algorithm-for-3d-vectors-intersection
        

##For example:
#r1=T.FeMd[idx,:]    #MidPoint of the edge
#r2=Pb[idx,:]        #CornerPoint
#e1=FeM2Ev[idx,:]
#e2=FePb2PcV[idx,:];

u =dot(e1,e2)
t1=dot(r2-r1,e1)
t2=dot(r2-r1,e2)
d1 = (t1-u*t2)/(1.0-u*u)
p1=r1.+d1.*e1
#d2 = (t2-u*t1)/(u*u.-1.0)
#p2=r2.+d2.*e2
if u==1
    println(p1)
    println(p2)
    println("lines are parallel") 
end

#Intersection [x y z]
return p1

end