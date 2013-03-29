%%
clear csim
clear field
sx = 500;
sy = 500;
hood = 3;
states = 2;
field = [];
hoodlim = 8;
colorsc.cmap = 'gray';
colorsc.caxis = [1 states];
epoch = 50;
simstats = struct;
simstats.diffs = [];
simstats.go = {};
simstats.go_dec = [];
simstats.nhood_dec = [];
simstats.nhood = [];
show = 1;
dump = 0;


i = 1;

while i < 10;
  close all
  field = [];
  field(1,:,:) = ones(sy,sx);
  field(1, 4*sy/9:5:5*sy/9,4*sx/9:5:5*sx/9) = 2;
  
  % seeds with a box shape
%   seedstates = states;
%   nseeds = 4;
%   seedside = [1:2:7; 1:2:7]';
%   margx = sx/5;
%   
%   seedx = round(margx+linspace(0,sx-2*margx,nseeds));
%   seedy = round(repmat(sy/2,[1 nseeds]));
%   
%   for s = 1:nseeds
%     seed = ones(seedside(s),seedside(s))+(states-1);
%     seed(2:size(seed,1)-1,2:size(seed,2)-1) = 1;
%     field(1,seedy(s)+1:seedy(s) + seedside(s,1),seedx(s)+1:seedx(s) + seedside(s,2)) = seed;
%   end
  
  field(2:epoch,:,:) = 0;
  fresh = 0;
  
  old_nhoods = simstats.nhood_dec;
  old_go = simstats.go_dec;
  
  while fresh == 0
    nhood = zeros(hood,hood);
    while ~any(nhood)
      hoodseed = randperm((hood^2));
      hoodseed = hoodseed(1:randi(hoodlim));
      nhood(hoodseed) = 1;
      nhood(median(1:hood^2)) =0;
    end
    
    trigs = randperm(sum(sum(nhood == 1)));
    go = randi(sum(nhood(:) == 1),1);
    go = trigs(1:go);
    
    temp = zeros(max(go),1);
    temp(go) = 1;
    temp = bin2dec(num2str(temp)');
    
    if ~any(find(old_nhoods == bin2dec(num2str(nhood(:))')))
      if ~any(find(old_go == temp))
        fresh = 1;
    end
    
    end
  end
  % Now check that we haven't run this nhood + go combination already, or a
  % rotation of the former
  
  csim = caut(field,nhood,go,states,show,colorsc);
  csim = csim.runSim;
  if size(csim.field,1) == epoch & dump ==1
    test = squeeze(sum(sum(permute(csim.field,[3 2 1]) == states)));
    % check if coefficient of variation is decreasing with time and
    % whether the derivative of the sum is constant
    cv = (std(test(1:10))/mean(test(1:10)) > std(test(end-9:end))/mean(test(end-9:end)));
    simstats.diffs(i,:) = diff(test);
    simstats.nhood(i,:) = csim.nhood(:);
    simstats.nhood_dec(i) = bin2dec(num2str(csim.nhood(:))');
    temp = zeros(max(csim.go),1);
    temp(csim.go) = 1;
    % Store a binary representation of the transition rule
    simstats.go_dec(i) = bin2dec(num2str(temp)');
    simstats.go(i,:) = temp;
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
  i = i + 1;
end
