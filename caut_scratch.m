%%
% general simulation area
clear csim
clear field
sx = 300;
sy = 300;
hood = 7;
states = 2;
field = [];
hoodlim = 15;
colorsc.cmap = 'gray';
colorsc.caxis = [1 states];
epoch = 150;
simstats = struct;
simsummary = {};
dump = 1;
if dump == 1;
  outpath = uigetdir;
  outpath = [outpath filesep];
end

%%
for i = 1:5000
  % field(1,:,:) = randi(3,side);
    %     mid = [round(sidey/2) round(sidex/2)];
    %     field(1,mid(1):mid(1) + hood-1,mid(2):mid(2)+hood-1) = (states-1)+ones(hood,hood);
    % tile field with nstates
    field(1,:,:) = ones(sy,sx);
    disp('c1')
    % make a matrix full of random ints, then sparsen it
%     matscale = [40 80];
%     randmat = randi(states,matscale(1),matscale(2));
%     fatty = kron(randmat,ones(sy/matscale(1),sx/matscale(2)));
%     while size(fatty,1) < sy || size(fatty,2) < sx
%         fatty = repmat(fatty,2,2);
%     end
%     %
    seedstates = states;
    nseeds =1;
    seedside = [3 3];
    barsz = [sy 3];
    % barrier blocking simulation growth
    % upper left corner of barrier
    barp = [1 1];
    seed = ones(seedside(1),seedside(2))+(states-1);
    barrier = ones(barsz)-2;
    field(1,sy/2+1:sy/2 + seedside(1),4:seedside(2)+3) = seed;
    field(1,barp(1):barp(1) + barsz(1) -1,barp(2):barp(2)+barsz(2)-1) = barrier;
    disp('c2')
    %     for i = 1:nseeds
%         tempx = randi(seedside(2) - seedside(1))+seedside(1);
%         tempy = randi(seedside(2) - seedside(1))+seedside(1);
%         tx = randi(sx - seedside(2));
%         ty = randi(sy-seedside(2));
%         field(1,ty:ty+tempy-1,tx:tx + tempx -1) = randi(seedstates-1) + ones(tempy,tempx);
%     end
    
    % boxdim = [2,2];
    % boxpos = [sy/2,sx/2];
    % randbar = randi(states,boxdim(1),boxdim(2));
    % field(1,boxpos(1):boxpos(1)+boxdim(1)-1,boxpos(2):boxpos(2)+boxdim(2)-1) = randbar;
    
    % field(1,:,:) = fatty(1:sidey,1:sidex);
    
%     randpos = randperm(sx-55);
%     randx = randpos(1);
%     randy = randpos(2);
%     randxw = randi(sx/10);
%     field(1,1:end,randx:randx+randxw) = randi(states,side,randxw+1);
   
    field(2:epoch,:,:) = zeros(epoch-1,size(field,2),size(field,3));
    disp('c2.1')
    % confine the number of 'on' sites in the neighborhood to 8
    hoodseed = randperm((hood^2));
    disp('c2.2')
    disp('c2.3')
    hoodseed = hoodseed(1:randi(hoodlim));
    % impose a directional bias in nhood
    % nhood = randi(2,hood,hood)-1;
    nhood = zeros(hood,hood);
    disp('c2.4')
    nhood(hoodseed) = 1;
    nhood(median(1:hood^2)) =0;
    %   while sum(nhood(1:median(numel(nhood)))) < 10/dbias * sum(nhood(median(1:numel(nhood))+1:end))
    %   nhood = randi(2,hood,hood)-1;
    %   end
    
    % while sum(sum(nhood)) < 1
    %    nhood = randi(2,hood,hood)-1;
    % end
    disp('c3')
    trigs = randperm(sum(sum(nhood == 1)));
    go = randi(sum(sum(nhood == 1)),1);
    go = trigs(1:go);
    disp('c4')
    stay = [1:(numel(nhood)-1)];
    stay = setdiff(stay,go);
    show = 1;
    dump = 0;
    disp('c5')
    %   field = int8(field);
    csim = caut(field,nhood,states,go,show,colorsc);
    csim = csim.runSim;
    simsummary{i,1} = squeeze(csim.field(end,:,:));
    simsummary{i,2} = csim.nhood;
    simsummary{i,3} = csim.go;
    simsummary{i,4} = i;
    if size(csim.field,1) == epoch & dump ==1
      set(gcf,'position',[100 10 800 800],'color','k')
      axis equal
      subplot('position',[.05 .05 .05 .05]);
      axis off
      imagesc(nhood);
      set(gca,'color','r');
      title(num2str(go),'color','r','fontsize',6)
      outfile = [csim.simid '.mat'];
      disp(['Saving ' outpath outfile])
      save([outpath outfile],'csim')
      figfile = csim.simid;
      print(gcf,[outpath figfile],'-dpng')
      close all
      success = 0;
      while success == 0;
        try
          csim.makeGif
          success = 1;
        catch
          success = 0;
        end
      end
      %     allsims{end+1} = csim;
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
