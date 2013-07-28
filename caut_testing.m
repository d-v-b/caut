%%
% width of the image
fx = 900;
% height of the image
fy = 400;
% width of the buffer on the sides  -- raise this value to have a smaller
% portion of the image seeded
buffx = 440;
buffy = 1;
% number of timesteps to run 
ft = 200;
% initialize the image matrix
f = ones(ft,fy,fx);
% size of the neighborhood around each pixel that influences its color
hood = 3;
% number of colors
states = 3;
% range of x values to use for generating the seed
xrange = [1+buffx:fx-buffx+1];
% range of y values to use for generating the seed
yrange = [1+buffy:fy];
% make the seed, which is usually just a random
seed = randi(states,size(f(1,yrange,xrange)));
% insert the seed
f(1,yrange,xrange) = seed;

% the simulation won't run unless all the images after the first are 0
f(2:end,:,:) = 0;
% here I generate random rules and neighborhoods for the simulation
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
% caut is the object I wrote that does all the work
sim = caut(f,n,g);  
% calling sim.runSim makes the simulation run, and puts an image of it on
% the screen
sim = sim.runSim;
