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
ex = csim.extend([25,250,150]);
% sy = size(ex.field,2);
% sx = size(ex.field,3);
% ex.field(1,end-1:end,1:2) = ex.field(1,1+(sy/2):2+(sy/2),1+(sy/2):2+(sy/2));
% ex.field(1,1+(sy/2):2+(sy/2),1+(sy/2):2+(sy/2)) = ones(size(ex.field(1,1+(sy/2):2+(sy/2),1+(sy/2):2+(sy/2))));
ex.field(2:end,:,:) = ex.field(2:end,:,:)*0;
ex = ex.runSim;

terp.a = interp1(ex.field,.05:0.05:25);
% terp.b = interp1(ex.field,50.05:0.05:100);
% for i = 1:2*size(ex.field,1) - 1
% terp(i,:
% end

%%
colormap(lbmap(256,'redblue'))
% colormap(gray(256));
for i = 1:size(terp.a,1)
  offset = [];
  
  four = log(abs(fftshift(fft2(squeeze(terp.a(i,:,:))))));
  four = four(4:end-3,4:end-3);
  imagesc(four);
  offset = .5* (size(terp.a,2) - size(four,1));
  imge = squeeze(terp.a(i,offset+1:end-offset,offset+1:end-offset));
  
  four(imge > min(min(imge))) = -NaN;
  %    imagesc(log(abs(fftshift(fft2(squeeze(conv2(upsc,kern,'valid')))))));
  %   imagesc(conv2(log(abs(fftshift(fft2(squeeze(squeeze(terp.a(i,:,:))))))),kern,'valid'));
  imagesc(four);
  %  imagesc(log(abs(fftshift(fft2(squeeze(imge(:,:)))))));
  % title(num2str(caxis))
  
  caxis([-8 5])
  % colorbar
  pause(.1)
  

end
