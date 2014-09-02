%%
% width of the image
fx = 480;
% height of the image
fy = 320;

% width of the buffer on the sides  -- raise this value to have a smaller
% portion of the image seeded
x = [150,152];
y = [150,152];

% number of timesteps to run
ft = 800;
% initialize the image matrix
f = ones(fy,fx,ft);

% Shape of the field
shpe = 'torus';

% size of the neighborhood around each pixel that influences its color
hood = 3;

% number of colors
states = 8;

% no. of nonzero neighborhood spots per timepoint
nrange = [4];

% History sensitivity of the neighborhoods
time = 4;

% acceptable "go" values, needs to be within [1, number of 1s in
% neighborhood kernel];
goRange(1) = 2;
% to prevent immersed cells from advancing
goRange(2) = nrange*time-1;

% symmetrical neighborhoods, may break due to change to nrange
symm = 0;
% range of x values to use for generating the seed
xrange = [x(1):x(2)];
% range of y values to use for generating the seed``
yrange = [y(1):y(2)];

% Number of frames between simulated frames for smoothing purposes
smoothing = 1;

tosave = 0;
outDir = '/Users/bennettd/Desktop/tmp/';
% how much to scale up the image size when saving  
outImScale = 4;

% todo(dvb) add magnitude of noise to add to frames of the simulation for display purposes
% todo(dvb) option to kill simulations that are almost certainly boring
% todo(dvb) figure out what makes a simulation interesting
figure(1);
clf
set(gcf,'color','k')
set(gcf,'renderer','opengl')
cmap = gray(256).^2;
cmap = cmap(round(linspace(1,length(cmap),states)),:);
colormap(cmap);

for j = 1:4
    % make the seed, which is usually just a random box 
    seed = randi(states,size(f(yrange,xrange,1:time)));
    % insert the seed
    
    f(yrange,xrange,1:time) = seed;
%    f(yrange,xrange,1:time) = a(yrange,xrange,1:time);
    f(:,:,time+1:end) = 0;
   
    % here I generate random rules and neighborhoods for the simulation
    n = {};
    % maximum number of neighbors
    
    for i =1:states
        n{i} = zeros([hood,hood,time]);
        while isempty(find(n{i}))
            % force radial symmetry by generating a (hood-1)/2 X (hood+1)/2 subunit
            % and tiling it around the center
             if symm == 1;
%                 for m = 1:time
%                 tile = zeros((hood-1)/2,(hood+1)/2);
%                 randNb = randperm(numel(tile));
%                 tile(randNb(1:randi(round(max(nrange)/4)))) = 1;
%                 n{i}(m,1:(hood-1)/2,1:(hood+1)/2) = tile;
%                 n{i}(m,1:(hood+1)/2,1+(hood+1)/2:end) = flipdim(tile,1)';
%                 n{i}(m,1+(hood+1)/2:end,1+(hood-1)/2:end) = flipdim(flipdim(tile,2),1);
%                 n{i}(m,1+(hood-1)/2:end,1:(hood-1)/2) = flipdim(tile,2)';
%                 end
             else
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
    end
    g= {};
    for i = 1:states
        randgo = randperm(numel(find(n{i})));
        randgo(randgo < goRange(1)) = [];
        randgo(randgo > goRange(2)) = [];
        g{i} = randgo(1:randi(numel(randgo)));
    end
    
    clear sim
    % caut is the object I wrote that does all the work
    sim = caut(f,n,g);
    sim.smooth = smoothing;
    sim.fieldshape = shpe;
    % calling sim.runSim makes the simulation run, and puts an image of it on
    % the screen
    sim = sim.runSim;
    
    set(gcf,'paperpositionmode','auto')
    axis off
    if tosave == 1;
    if size(sim.field,3) > 150
       % save simulation as tif stack
       outFName = [num2str(states) '_' date '_' num2str(j) '.tif'];
       kronKern = ones(outImScale);
       imwrite(kron(sim.field(:,:,1),kronKern), cmap, [outDir outFName],'tif')   
    for k = 2:size(sim.field,3)
       imwrite(kron(sim.field(:,:,k),kronKern),cmap, [outDir outFName],'tif', 'writemode', 'append');
    end   
%       print(gcf,[num2str(states) '_' date '_' num2str(j)],'-r0','-dpng')
    end
    end
end