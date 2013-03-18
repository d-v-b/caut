%% Charlie Music zone

%%
% playback w/filters -- sims look better a bit anti-aliased, how can we get
% that to haps
% crude-ass AA filter
% generate sinc kernel
outpath = 'c:\users\lol\desktop\plots\charlie3';
kern = [1 4 7 4 1; 4 16 26 16 4; 7 26 41 26 7; 4 16 26 16 4; 1 4 7 4 1]/273;

if sum(size(csim.nhood)) > 2*sqrt(numel(csim.nhood))
  csim.nhood = reshape(csim.nhood,[sqrt(numel(csim.nhood)),sqrt(numel(csim.nhood))]);
end
ex = csim;
%%%%%%%%%%%%%%%%%%%%%%% 
% If we want to re-make the seed with the same rules/nhood
%%%%%%%%%%%%%%%%%%%%%%%
ex.field(1,:,:) = ones(size(ex.field(1,:,:)));
ex.field(2:end,:,:) = ex.field(2:end,:,:) * 0;
sy = size(ex.field,2);
sx = size(ex.field,3);
% start with a 5x5
szy = 1;
szx = 1;
ex.field(1,sy/2:sy/2 + (szy -1),sx/2:sx/2 + (szx -1)) = 1+ones(szy,szy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ex = ex.extend([100,300+(300-size(ex.field,2)),300+(300-size(ex.field,3))]);
% ex.field(2:end,:,:) = ex.field(2:end,:,:)*0;
ex = ex.runSim;
%%
terp.a = interp1(ex.field,0:.1:15);
% terp.b = interp1(ex.field,50.05:0.05:100);
% for i = 1:2*size(ex.field,1) - 1
% terp(i,:
% end

%%

figure;
close all
set(gcf,'renderer','zbuffer')
colormap(lbmap(1256,'bluegray'))
for i = 1:size(terp.a,1)
  offset = [];
  
  four = log(abs(fftshift(fft2(squeeze(terp.a(i,:,:))))));
  four = four(4:end-3,4:end-3);
%   imagesc(four);
%   offset = .5* (size(terp.a,2) - size(four,1));
%   imge = squeeze(terp.a(i,offset+1:end-offset,offset+1:end-offset));
  
%   four(imge > min(min(imge))) = -NaN;
  %    imagesc(log(abs(fftshift(fft2(squeeze(conv2(upsc,kern,'valid')))))));
  %   imagesc(conv2(log(abs(fftshift(fft2(squeeze(squeeze(terp.a(i,:,:))))))),kern,'valid'));
  pcolor(four);
  shading interp
  %  imagesc(log(abs(fftshift(fft2(squeeze(imge(:,:)))))));
  % title(num2str(caxis))
  
  caxis([-4 5.5])
  % colorbar
  pause(.01)
end
clear four
clear imge

%%
figure;
close all
set(gcf,'renderer','zbuffer','menubar','none',...
  'toolbar','none','units','normalized','outerposition',[-.01 -.01 1.3 1.3])
colormap(hot(3000));
% screensize = get(0,'monitorPositions');
% screensize = screensize(3:4);
% set(gcf,'position',[-1 -1 screensize(1) screensize(4) + 50])
for i = 1:size(terp.a,1)
  imagesc(squeeze(terp.a(i,:,:)));
  set(gcf,'color','k')
  
%   set(gcf,'position',[1 1 screensize screensize]) 
  
  %  imagesc(log(abs(fftshift(fft2(squeeze(imge(:,:)))))));
  % title(num2str(caxis))
  
%   caxis([-4 8])
  % colorbar
  pause(.1)
end