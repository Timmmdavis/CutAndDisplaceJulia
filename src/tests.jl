function test()
A=zeros(50,50);
B=view(A,:,1); # Take slice
C=testFunc(B);
end

function testFunc(X)
D=1.0;
@time testFunc2(X,D)
end

function testFunc2(B,D)
h=D
end

