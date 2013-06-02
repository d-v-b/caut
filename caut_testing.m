%%
fxy = 400;
mid = fxy/2;
ft = 500;
f = ones(ft,fxy,fxy);
f(1,mid:mid+3,mid:mid+3) = 2;
f(1,mid-1:mid+1,mid-1:mid+1) = 3;
f(2:end,:,:) = 0;
n = {};
n{1} = ones(1,3,3);
n{2} = ones(1,3,3);
n{2}(1) = 0; n{2}(5) = 0; n{2}(7) = 0;
n{2}(3) = 0; n{2}(9) = 0;
n{3} = n{1};
n{3}(2:2:8) = 0;
g= {};
g{1} = [1,2];
g{2} = [1,2];
g{3} = [1,2];
%%
clear sim
sim = caut(f,n,g);  
sim = sim.runSim;
%%