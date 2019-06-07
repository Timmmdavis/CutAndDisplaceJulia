function CreateMidPoint!(Pa,Pb,Pc,MidPoint)

#Distance formula between two points.   
c=sqrt((Pa[1]-Pb[1])^2+(Pa[2]-Pb[2])^2+(Pa[3]-Pb[3])^2);
a=sqrt((Pb[1]-Pc[1])^2+(Pb[2]-Pc[2])^2+(Pb[3]-Pc[3])^2);
b=sqrt((Pa[1]-Pc[1])^2+(Pa[2]-Pc[2])^2+(Pa[3]-Pc[3])^2);
if a==0 || b==0 || c==0
	error("Slither triangles or triangles with no area exist on your surface")
end
#Calculating midpoint using
#http://mathworld.wolfram.com/Incenter.html
#Gives the same result as MATLABS incenter calc
MidPoint[1]=((Pa[1]*a)+(Pb[1]*b)+(Pc[1]*c))/(a+b+c);
MidPoint[2]=((Pa[2]*a)+(Pb[2]*b)+(Pc[2]*c))/(a+b+c);
MidPoint[3]=((Pa[3]*a)+(Pb[3]*b)+(Pc[3]*c))/(a+b+c);

return MidPoint

end