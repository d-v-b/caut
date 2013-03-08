%% Charlie Music zone
% Need ~ 20m of simulations that are visually interesting
%
sx = 300;
sy = 300;
hood = 3;
states = 2;
field = [];
hoodlim = 9;
colorsc.cmap = 'gray';
colorsc.caxis = [1 3];
epoch = 20;
field(1:epoch,:,:) = zeros(epoch,sx,sy);
field(1,:,:) =  field(1,:,:) + 1;
seed = ones(3)+1;
seedy = round(sy/2);
seedx = round(sx/2);
field(1,seedx:seedx + size(seed,1)-1,seedy:seedy+size(seed,1)-1) = seed;
hoodseed = randperm((hood^2));
hoodseed(hoodseed == median(1:hood^2)) = [];
hoodseed = hoodseed(1:randi(hoodlim));
% impose a directional bias in nhood
% nhood = randi(2,hood,hood)-1;
nhood = zeros(hood,hood);
nhood(hoodseed) = 1;
%   while sum(nhood(1:median(numel(nhood)))) < 10/dbias * sum(nhood(median(1:numel(nhood))+1:end))
%   nhood = randi(2,hood,hood)-1;
%   end

% while sum(sum(nhood)) < 1
%    nhood = randi(2,hood,hood)-1;
% end
trigs = randperm(sum(sum(nhood == 1)));
go = randi(sum(sum(nhood == 1)),1);
go = trigs(1:go);
stay = [1:(numel(nhood)-1)];
stay = setdiff(stay,go);
show = 1;
sim = caut(field,nhood,states,go,show,colorsc);
sim = sim.runSim;

%%
sim2 = sim.extend(100);

%%
outpath = uigetdir;
outpath = [outpath filesep];

% general simulation area
clear p
clear field
sx = 1600;
sy = 800;
hood = 3;
states = 2;
field = [];
hoodlim = 15;
colorsc.cmap = 'gray';
colorsc.caxis = [1 states];
epoch = 15;

for i = 1
    % field(1,:,:) = randi(3,side);
    %     mid = [round(sidey/2) round(sidex/2)];
    %     field(1,mid(1):mid(1) + hood-1,mid(2):mid(2)+hood-1) = (states-1)+ones(hood,hood);
    % tile field with nstates
    field(1,:,:) = ones(sy,sx);
    % make a matrix full of random ints, then sparsen it
    matscale = [40 80];
    randmat = randi(states,matscale(1),matscale(2));
    fatty = kron(randmat,ones(sy/matscale(1),sx/matscale(2)));
    while size(fatty,1) < sy || size(fatty,2) < sx
        fatty = repmat(fatty,2,2);
    end
    %
    seedstates = states;
    nseeds =20000;
    seedside = [2 15];
    
    for i = 1:nseeds
        tempx = randi(seedside(2) - seedside(1))+seedside(1);
        tempy = randi(seedside(2) - seedside(1))+seedside(1);
        tx = randi(sx - seedside(2));
        ty = randi(sy-seedside(2));
        field(1,ty:ty+tempy-1,tx:tx + tempx -1) = randi(seedstates-1) + ones(tempy,tempx);
    end
    
    % boxdim = [2,2];
    % boxpos = [sy/2,sx/2];
    % randbar = randi(states,boxdim(1),boxdim(2));
    % field(1,boxpos(1):boxpos(1)+boxdim(1)-1,boxpos(2):boxpos(2)+boxdim(2)-1) = randbar;
    
    % field(1,:,:) = fatty(1:sidey,1:sidex);
    
    randpos = randperm(sx-55);
    randx = randpos(1);
    randy = randpos(2);
    randxw = randi(sx/10);
    
    % field(1,1:end,randx:randx+randxw) = randi(states,side,randxw+1);
    
    epoch = 50;
    
    field(2:epoch,:,:) = zeros(epoch-1,size(field,2),size(field,3));
    
    % confine the number of 'on' sites in the neighborhood to 8
    hoodseed = randperm((hood^2));
    hoodseed(hoodseed == median(1:hood^2)) = [];
    hoodseed = hoodseed(1:randi(hoodlim));
    % impose a directional bias in nhood
    % nhood = randi(2,hood,hood)-1;
    nhood = zeros(hood,hood);
    nhood(hoodseed) = 1;
    %   while sum(nhood(1:median(numel(nhood)))) < 10/dbias * sum(nhood(median(1:numel(nhood))+1:end))
    %   nhood = randi(2,hood,hood)-1;
    %   end
    
    % while sum(sum(nhood)) < 1
    %    nhood = randi(2,hood,hood)-1;
    % end
    
    trigs = randperm(sum(sum(nhood == 1)));
    go = randi(sum(sum(nhood == 1)),1);
    go = trigs(1:go);
    
    stay = [1:(numel(nhood)-1)];
    stay = setdiff(stay,go);
    show = 1;
    dump = 0;
    %   field = int8(field);
    %%
    p = caut(field,nhood,states,go,show,colorsc);
    p = p.runSim;
    outfile = [p.simid '.mat'];
    if size(p.field,1) > 10
        disp(['Saving ' outpath outfile])
        save([outpath outfile],'p')
        figfile = [p.simid '_t_' num2str(i)];
        print(gcf,[outpath figfile],'-dpng')
        if size(p.field,1) >= 50
            imagesc(squeeze(p.field(randi(size(p.field,1)),:,:)))
            caxis(p.colorsc.caxis)
            axis off
            axis image
            colormap gray
            set(gcf,'color','k');
            figfile = [p.simid '_t_' num2str(i)];
            print(gcf,[outpath figfile],'-dpng')
        end
    else
        %      disp(['Not saving ' outpath outfile])
    end
    
end
%%
% animations

[infile inpath] = uigetfile('multiselect','on');
if ~iscell(infile)
    infile = {infile};
end
outpath = uigetdir;
%%
for xa = 1:numel(infile)
    load([inpath infile{xa}])
    p.colorsc = 'gray'
    sepind = strfind(infile{xa},'.');
    outname = [outpath filesep 'aut'  infile{xa}(1:sepind-1)]
    if size(p.field,1) > 100
        p.makeGif(outname)
    end
end
%%
figure(1)
filename = [p.simid 'surfwiener2.gif'];
t = surf(wiener2(double(squeeze(p.field(1,:,:)))));
shading interp
set(gca,'cameraposition',[177.9686 81.1290 915.4978])
set(gca,'cameraviewangle',4)
for n = 1:size(p.field,1)
    set(t,'zdata',wiener2(squeeze(double(p.field(n,:,:)))))
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if n == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
