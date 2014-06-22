%%
% width of the image
fx = 200;
% height of the image
fy =200;
% width of the buffer on the sides  -- raise this value to have a smaller
% portion of the image seeded
x = [100,100];
y = [100,100];

% number of timesteps to run
ft = 200;
% initialize the image matrix
f = ones(fy,fx,ft);
% Shape of the field
shpe = 'torus';
% size of the neighborhood around each pixel that influences its color
hood = 3;
% number of colors
states = 16;
% intial distribution of colors
seedVals = [1,2];
% no. of nonzero neighborhood spots per time
nrange = [8];
% History sensitivity of the neighborhoods
time = 4;
% permissible transition rules
gorange = [1:time*nrange-1];
%gorange = [1,2,3];
% symmetrical neighborhoods, may break due to change to nrange
symm = 0;

% range of x values to use for generating the seed
xrange = [x(1):x(2)];
% range of y values to use for generating the seed
yrange = [y(1):y(2)];
% whether seed will be tiled
kronMode = 0;
% Scaling factor in y and x
tileSize = [3,3];

% Number of frames between simulated frames for smoothing purposes
smoothing = 1;
% todo(dvb) add magnitude of noise to add to frames of the simulation for display purposes
% todo(dvb) option to kill simulations that are almost certainly boring
% todo(dvb) figure out what makes a simulation interesting
% todo(dvb) various optimization/visualization tricks -- 

% Whether or not to save the simulation
toSave = 0;
% Whether or not to re-use the old simulation as the initialization of the
% next
fieldLoop = 0;

figure(1);
clf

set(gcf,'color','k')
set(gcf,'renderer','opengl')
%colormap(circshift(bone,[0,1]).^1.5)
colormap hot
for j = 1:200    
    if fieldLoop == 1 && j > 1
        f = sim.field(:,:,end-time+1:end);                    
    end
    % make a random seed, which is usually just a random box
    seed = [];
    if kronMode == 1
        % if we use kronmode, then we are taking a random seed and
        % upscaling it to fill the seed area
        for t = 1:time
            seed(:,:,t) = kron(double(randi(seedVals,floor([length(yrange)/tileSize(1),length(xrange)/tileSize(2)]))),ones(tileSize(1),tileSize(2)));
        end
    else
        seed = double(randi(seedVals,size(f(yrange,xrange,1:time))));
    end
    % insert the seed        
        
    f(yrange,xrange,1:time) = seed;    
    f(:,:,time+1:ft) = 0;
   
    % here I generate random rules and neighborhoods for the simulation
    n = {};
    % maximum number of neighbors

    for i =1:states
        n{i} = zeros([hood,hood,time]);
        while isempty(find(n{i}))
                for p = 1:time
                    tmpNb = zeros(hood,hood);
                    randNb = randperm(hood^2);
                    % remove center of future-most neighborhood slice
                    randNb(randNb == median(1:hood^2)) = [];
                    tmpNb(randNb(1:nrange(randi(numel(nrange))))) = 1;
                    n{i}(:,:,p) = tmpNb;
                end            
        end
    end
    g= {};
    for i = 1:states
        randgo = randperm(numel(gorange));
        g{i} = gorange(randgo(1:randi([1,numel(gorange)])));
    end
        
    clear sim
    % caut is the object I wrote that does all the work
    sim = caut(f,n,g);
    sim.smooth = smoothing;
    sim.fieldshape = shpe;
    % calling sim.runSim makes the simulation run, and puts an image of it on
    % the screen
    sim = sim.runSim;
       
    if size(sim.field,3) > 150 && toSave == 1
       % set most frequent color to bottom of colormap
        backg = sim.field(:,:,end);
        backg(backg == mode(backg(1,:))) = 1;
        %imagesc(backg)
       set(gcf,'paperpositionmode','auto')
       set(gcf, 'InvertHardCopy', 'off');
       axis off
       print(gcf,[num2str(states) '_' date '_' num2str(j)],'-r0','-dpng')
    end
end
