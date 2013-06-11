%%
fx = 900;
fy = 400;
buffx = 440;
buffy = 1;
ft = 200;
f = ones(ft,fy,fx);
hood = 3;
states = 3;
xrange = [1+buffx:fx-buffx+1];
yrange = [1+buffy:fy];
seed = randi(states,size(f(1,yrange,xrange)));
f(1,yrange,xrange) = seed;

%  for p = 1:5000
f(2:end,:,:) = 0;
n = {};
for i =1:states
    n{i} = randi(2,[1,hood,hood])-1;
    n{i}(1,median(1:hood),median(1:hood)) = 0;
    while isempty(find(n{i}))
        n{i} = randi(2,[1,hood,hood])-1;
        n{i}(1,median(1:hood),median(1:hood)) = 0;
    end
end
g= {};
for i = 1:states
    randgo = randperm(numel(find(n{i})));
    g{i} = randgo(1:randi(numel(find(n{i}))));
end
    
clear sim

sim = caut(f,n,g);  
sim = sim.runSim;
f(1,:,:) = sim.field(end,:,:);
%  end