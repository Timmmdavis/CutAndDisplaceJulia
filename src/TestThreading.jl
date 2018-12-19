function TestThreading()

FillA=zeros(1000,1);
FillB=zeros(1000,1);

Threads.@threads for i=1:100000 #Threads.@threads 

	e12=[0 0 1];
	e13=[0 1 0];

	A=e12'*e13;
	A=acos(A[1]);
	
	B=(e12[1]*e13[1])+(e12[2]*e13[2])+(e12[3]*e13[3]);
	B=acos(B);

	
end	

if any(FillA-FillB.>0)
	error("Two methods not equal")
end

end