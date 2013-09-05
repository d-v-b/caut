%%
% width of the image
fx = 600;
% height of the image
fy = 400;
% width of the buffer on the sides  -- raise this value to have a smaller
% portion of the image seeded
x = [1 600];
y = [375 400];

% number of timesteps to run 
ft = 200;
% initialize the image matrix
f = ones(ft,fy,fx);
% size of the neighborhood around each pixel that influences its color
hood = 7;
% number of colors
states = 6;
% range of x values to use for generating the seed
xrange = [x(1):x(2)];
% range of y values to use for generating the seed
yrange = [y(1):y(2)];

for j = 1:2000
% make the seed, which is usually just a random box
seed = randi(states,size(f(1,yrange,xrange)));
% insert the seed
f(1,yrange,xrange) = seed;

% the simulation won't run unless all the images after the first are 0
f(2:end,:,:) = 0;
% here I generate random rules and neighborhoods for the simulation
n = {};
% maximum number of neighbors
nmax = 18;
symm = 0;
for i =1:states
    n{i} = zeros([1,hood,hood]);
    while isempty(find(n{i}))
        % force radial symmetry by generating a (hood-1)/2 X (hood+1)/2 subunit
        % and tiling it around the center
        if symm == 1;
            tile = zeros((hood-1)/2,(hood+1)/2);
            randNb = randperm(numel(tile));
            tile(randNb(1:randi(round(nmax/4)))) = 1;
            n{i}(1,1:(hood-1)/2,1:(hood+1)/2) = tile;
            n{i}(1,1:(hood+1)/2,1+(hood+1)/2:end) = flipdim(tile,1)';
            n{i}(1,1+(hood+1)/2:end,1+(hood-1)/2:end) = flipdim(flipdim(tile,2),1);
            n{i}(1,1+(hood-1)/2:end,1:(hood-1)/2) = flipdim(tile,2)';
            
        else
        randNb = randperm(numel(n{i}));
        randNb(randNb == median(1:(hood^2))) = [];
        n{i}(randNb(1:randi(nmax))) = 1;
        end
    end
end
g= {};
for i = 1:states
    randgo = randperm(numel(find(n{i})));
    g{i} = randgo(1:randi(numel(find(n{i}))));
end
    
% %% make a random colormap
% colmp = zeros(states,3);
% for i = 1:states
%     temp = randi(100,1,3);
%     colmp(i,:) = temp ./ max(temp);
% end
% %%

clear sim
% caut is the object I wrote that does all the work
sim = caut(f,n,g);  
% calling sim.runSim makes the simulation run, and puts an image of it on
% the screen
sim = sim.runSim;
set(gcf,'paperpositionmode','auto')
print(gcf,['E:\cauts\07282013\3n6nm4s' num2str(j)],'-r0','-dpng')
end