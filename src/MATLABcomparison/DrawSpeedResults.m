FS.NoTris=[2,10,100,2500];
FS.MATLAB=[0.013,0.072,0.907,186.39];
FS.Julia=[0.005,0.013,0.104,2.534];

HS.NoTris=[2,10,100,2500];
HS.MATLAB=[0.37,1.69,17.12,425.94];
HS.Julia=[0.06,0.12,0.87,21.38];

loglog(FS.NoTris,FS.MATLAB)
hold on
loglog(FS.NoTris,FS.Julia)
loglog(HS.NoTris,HS.MATLAB)
loglog(HS.NoTris,HS.Julia)

xlabel('NoTris')
ylabel('Time(s)')
legend MATLAB-FS Julia-FS MATLAB-HS Julia-HS

FS.Relative=FS.MATLAB./FS.Julia
HS.Relative=HS.MATLAB./HS.Julia
figure
hold on
plot(FS.NoTris,FS.Relative)
plot(HS.NoTris,HS.Relative)
xlabel('NoTris')
ylabel('Relative speed')
legend FullSpace HalfSpace