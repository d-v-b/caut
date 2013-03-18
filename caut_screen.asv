%%
% load/run simulations and stock a database with various parameters of the sims

clear csim
clear field
sx = 300;
sy = 200;
hood = 3;
states = 2;
field = [];
hoodlim = 9;
colorsc.cmap = 'gray';
colorsc.caxis = [1 states];
epoch = 50;
simstats = struct;
simstats.diffs = {};
simstats.go = {};
simstats.nhood = {};
simsummary = {};

count = 1;
%%
while count < 5000;
  % field(1,:,:) = randi(3,side);
  %     mid = [round(sidey/2) round(sidex/2)];
  %     field(1,mid(1):mid(1) + hood-1,mid(2):mid(2)+hood-1) = (states-1)+ones(hood,hood);
  % tile field with nstates
  field = [];
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
  nseeds = 4;
  seedside = [1:4; 1:4]';
  margx = sx/5;
  seedx = round(margx+linspace(0,sx-2*margx,nseeds));
  seedy = round(repmat(sy/2,[1 nseeds]));
  for s = 1:nseeds
    seed = ones(seedside(s),seedside(s))+(states-1);
    field(1,seedy(s)+1:seedy(s) + seedside(s,1),seedx(s)+1:seedx(s) + seedside(s,2)) = seed;
  end  
  
  
  disp('c2')
  
  
  field(2:epoch,:,:) = zeros(epoch-1,size(field,2),size(field,3));
  disp('c2.1')
  
  nhood = zeros(hood,hood);
  while ~any(nhood)
    hoodseed = randperm((hood^2));
    disp('c2.2')
    disp('c2.3')
    hoodseed = hoodseed(1:randi(hoodlim));
    % impose a directional bias in nhood
    % nhood = randi(2,hood,hood)-1;
    disp('c2.4')
    nhood(hoodseed) = 1;
    nhood(median(1:hood^2)) =0;
  end
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
  disp('c5')
  %   field = int8(field);
  csim = caut(field,nhood,states,go,show,colorsc);
  csim = csim.runSim;
  simsummary{i,1} = squeeze(csim.field(end,:,:));
  simsummary{i,2} = csim.nhood;
  simsummary{i,3} = csim.go;
  simsummary{i,4} = i;
  if size(csim.field,1) == epoch & dump ==1
    test = squeeze(sum(sum(permute(csim.field,[3 2 1]) == states)));
    % check if coefficient of variation is decreasing with time and
    % whether the derivative of the sum is constant
    cv = (std(test(1:10))/mean(test(1:10)) > std(test(end-9:end))/mean(test(end-9:end)));
    simstats.diffs{i} = test;
    simstats.nhood{} = csim.nhood;
    simstats.go{i} = csim.go;
    if any(diff(test)) & ~cv
        set(gcf,'position',[100 10 900 600],'color','k')
        subplot(1,2,1)
        imagesc(squeeze(csim.field(end,:,:)));
        axis equal
        axis tight
        subplot(1,2,2)
        imagesc(log(abs(fftshift(fft2(squeeze(csim.field(end,:,:)))))));
        axis equal
        axis tight
        subplot('position',[.05 .05 .05 .05]);
        axis off
        imagesc(nhood);
        axis square
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
    else
      % save image and file as junk
        set(gcf,'position',[100 10 900 600],'color','k')
        subplot(1,2,1)
        imagesc(squeeze(csim.field(end,:,:)));
        axis equal
        axis tight
        subplot(1,2,2)
        imagesc(log(abs(fftshift(fft2(squeeze(csim.field(end,:,:)))))));
        axis equal
        axis tight
        subplot('position',[.05 .05 .05 .05]);
        axis off
        imagesc(nhood);
        axis square
        set(gca,'color','r');
        title(num2str(go),'color','r','fontsize',6)
        outfile = [csim.simid '.mat'];
        disp(['Saving ' outpath outfile])
        save([outpath outfile],'csim')
        figfile = [csim.simid '_loser'];
        print(gcf,[outpath figfile],'-dpng')
        close all
    end
  
  end
    count = count + 1;
  end
