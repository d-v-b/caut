% code for adding simulation parameters to a database object
% with the following fields:
% neighborhood, represented as a string, which is hard to
% visualize
% path to all data generated from this rule + neighborhood combination
% this object can be used to avoid running redundant simulations, and to
% help guide targeted future screens.

% current need: load a host of .mat files and add their params to a
% database, render all simulations somehow fit for manipulation

% only worried about 2-state systems for now

[infiles inpath] = uigetfile('.mat','Pick simulation files','multiselect','on');
update = 1;
[datfile datpath] = uigetfile('.mat','Pick a database to update, or cancel to make a new one');
if isequal(0,datfile)
  update = 0;
  datpath = [uigetdir(inpath,'set path for new database file') filesep];
  datfname = 'cautDB.mat';
end
%%
clear predat
clear dat
clear tempool
% list of duplicate simulations
dupes = {};
% list of nhood:go combinations loaded
log = [0;0];

tic
disp(['Begin loading ' num2str(numel(infiles)) ' simulation files... ']);
% build aggregate data in a temporary variable to be integrated with the
% main database
for ifile = numel(infiles):-1:1
  % Set this variable to 1 when it's clear the simulation is new
  new = 0;
  try
    load([inpath infiles{ifile}]);
  catch err
    disp(err.message)
    disp('skipping this file')
    continue
  end
  
  % convert binary params to decimals to check for repetitions
  temphood = bin2dec(rot90(num2str(csim.nhood(:))));
  tempgo = zeros(max(csim.go),1);
  tempgo(csim.go) = 1;
  tempgo = bin2dec(num2str(tempgo)');
  
  % Check if these have already been added by checking if any instance of 
  % nhood in the log is on the same row as any instance of go
  if intersect(find(log(1,:) == temphood), find(log(2,:) == tempgo))
    dupes{1,end+1} = ifile;
    dupes{2,end+1} = min(find(log(1,:) == temphood));
  else
    % Add these values to log, and proceed to the next step
    new = 1;
    log = [log [temphood;tempgo]];
  end
  
  if new == 1
  % Store a binary representation of the transition rule
  tempool(ifile).nhood = csim.nhood(:);
  tempool(ifile).nhood_dec = temphood;
  tempool(ifile).go = csim.go;
  tempool(ifile).go_dec = tempgo;

  % the total number of cells in the upstate as a function of time
  tempool(ifile).pop = squeeze(sum(sum(permute(csim.field(:,:,:), [3 2 1]) == 2)));
  % the size of the simulation
  tempool(ifile).size = size(csim.field);
%   tempool(ifile).final = squeeze(csim.field(end,:,:));
  tempool(ifile).path = inpath;
  fname = infiles{ifile}(1:strfind(infiles{ifile},'.')-1);
  % search for any other files with the simid
  others = ls([inpath fname '*']);
  % concatenate each row into a cell array
  for k = 1:size(others,1)
    tempool(ifile).file{k} = deblank(strcat(others(k,:)));
  end
  end
  clear csim
end

disp('Done processing files.')

%%
% Sort simulations by neighborhood, then by go
nbs = {tempool(:).nhood};
% figure out how many distinct neighborhood sizes we have
num_nbs = cellfun(@numel,nbs);
unique_nbs = sort(unique(num_nbs));
unique_nbs(unique_nbs < 1) = [];
sortnbs = {};
sortgo = {};
% initialize the struct to hold the sims w/ different nhoods
predat = tempool(1);
for i = 1:numel(unique_nbs)
  % sort all the neighborhoods of each size and add them to sortnbs
  nbs_ids = find(num_nbs == unique_nbs(i));
  [~, sortnbs{i}] = sort([tempool(nbs_ids).nhood_dec]);
  sortnbs{i} = nbs_ids(sortnbs{i});
  % Sort by go rules converted to decimals
  [~, sortgo{i}] = sort([tempool(nbs_ids).go_dec]);
  sortgo{i} = nbs_ids(sortgo{i});
  % Now put all the sims with same sized neighborhoods in the same column
  % of a struct array
  predat(1:numel(sortnbs{i}),i) = tempool(sortnbs{i});
end
% If applicable update the database, then save
if update == 1
  % Combine tempool with the existing database
  load([datpath datfile])
  % for each cell in the new dataset, check and see if it has a match in
  % the old database. if so, then move on. otherwise, stick it in the
  % appropriate spot, given its position in nhood and go orders
  
  for i = 1:numel(tempool)
    oldhoods = {dat(:,:).nhood_dec};
    % Put the candidate in every position in a cell array of size
    % comparable to dat
    clones = cell(size(oldhoods));
    clones(:) = {tempool(i).nhood_dec};
    % check for nhood identity
    hits = find(cellfun(@isequal,oldhoods,clones));
    % Now check for go identity
    oldgos = [dat(hits).go_dec];
    if ~find(oldgos == tempool(i).go_dec);
      % If none are found, then stick tempool in the main DB
      % Find row position in dat:
      row = find(numel(tempool(i).nhood) == cellfun(@numel,{dat(1,:).nhood}));
      % Now stick it on the end, and then sort
      dat(end+1,row) = tempool(i);
      [~, sortnbs] = sort([dat(:,row).nhood_dec]);
      dat(1:numel(sortnbs),row) = dat(sortnbs,row);
    end
  end
else
  dat = predat;
  save([datpath datfname], 'dat');
  % predat is the new database
end

%%
% % now run some stats
% build matrix of time series of field populations
% find all the monotonically increasing sims
% take the diff of the last half of the simulation, check for lack of negative
% values
hoodsize = 9;
epoch = 50;
row = find(cellfun(@numel,{dat(1,:).nhood}) == hoodsize);
col = find(cellfun(@numel,{dat(:,row).pop}) == epoch);
pops = [dat(col,row).pop];

monocheck = min(diff(pops(end-10:end,:))) > 0;
compl = find(~monocheck);
mono = find(monocheck);
%%
% look at all the populations
close all
figure
plot(pops);


