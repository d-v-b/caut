
%%%%%%%%%%%%%%%%%%%%%%%
% Cellular automaton class
% Properties:
% Methods:

% --Davis Bennett--
%--  03/19/2013  --

%%%%%%%%%%%%%%%%%%%%%%%%

% TODO (DVB): increase dimensionality of nhood along time and state axis --
% different states can have different nhoods.
% TODO (DVB): increase dimensionality of go rules along state axis --
% different states can have different rules for advancing

classdef caut
    
    properties
        states % int: states for a cell
        nhood % array: neighbors around a cell. a matrix with a NaN at the cell
        % to be updated; values at all other positions determine the weight of
        % that position in the updating step
        go % cell: conditions for a cell to advance; if size(go,1) > 1,
        % then there are multiple go rules, at most one per state.
        field % array: the simulation itself
        show % struct: determines how simulation is displayed
        colorsc % colorscale for visualizing system
        fieldshape % geometry of the field -- torus or bounded plane
        %qualParams % struct containing parameters for the quality checks
        %behavior % struct containing boolean-valued fields that describe the behavior of the sim
    end
    
    properties (Dependent = true)
        simid % the identification code for this simulation
    end
    
    methods
        function id = get.simid(obj)
            % generate a unique ID for each simulation. The id must contain
            % information about the size of the field, the length of the
            % epoch, the total states, the go conditions, and the behavior
            % of the neighborhood
            stat = ['s' num2str(max(obj.states))];
            go = '-g';
            for i = 1:length(obj.go)
                if i == 1
                    go = strcat(go, [num2str(dec2hex(obj.go(i)))]);
                else
                    go = strcat(go, ['-' num2str(dec2hex(obj.go(i)))]);
                end
            end
            
            hood = '-n';
            k = num2str(find(obj.nhood(:)));
            for i = 1:size(k,1)
                hood = [hood '-' num2str(dec2hex(str2num(k(i,:))))];
            end
            
            fld = '-f';
            for i = 1:3
                if i == 1
                    fld = strcat(fld, [num2str(size(obj.field,i))]);
                else
                    fld = strcat(fld, ['-' num2str(size(obj.field,i))]);
                end
            end
            id = [stat go hood fld];
        end
        
        function obj = caut(fld_matrix, nbr_matrix, go_matrix, varargin)
            % Constructs the intial condition of the simulation
            par = inputParser;
            defcol.cmap = 'gray';
            defcol.caxis = [];
            par.addRequired('field')
            par.addRequired('nhood')
            par.addRequired('go')
            par.addParamValue('colorsc','gray')
            par.addParamValue('nstates',numel(unique(fld_matrix(1,:,:))))
            par.addParamValue('show',1)
            par.addParamValue('fieldshape','plane')
            par.parse(fld_matrix, nbr_matrix, go_matrix,varargin{:})
            
            obj.field = par.Results.field;
            obj.nhood = par.Results.nhood;
            obj.go = par.Results.go;
            obj.states = par.Results.nstates;
            
            % If the input rules weren't cell arrays, make them so.
            if ~iscell(obj.go)
                obj.go = {obj.go};
                obj.go(2:obj.states) = obj.go(1);
            end
            
            if ~iscell(obj.nhood)
                obj.nhood = {obj.nhood};
                obj.nhood(2:obj.states) = obj.nhood(1);
            end
            
            % At most there is a 1:1 rule:state pairing, where a rule is a
            % combination of go and nhood matrices. Store nhood and go as
            % cell arrays with length == nstates, with each entry
            % corresponding to the go rule/nhood for that state.
            
            if ~isequal(numel(obj.go), obj.states)
                if numel(obj.go) == 1
                    obj.go(2:nstates) = obj.go(1);
                else
                    error('Number of go conditions must equal number of states or 1')
                end
            end
            
            if ~isequal(numel(obj.nhood), obj.states)
                if numel(obj.nhood) == 1
                    obj.nhood(2:nstates) = obj.nhood(1);
                else
                    error('Number of neighborhoods must equal number of states or 1')
                end
            end
            
            obj.show = par.Results.show;
            obj.colorsc = par.Results.colorsc;
            obj.fieldshape = par.Results.fieldshape;
        end
        
        % method to run the simulation
        function obj =  runSim(obj)
            tic
            nstates = obj.states;
            % 1 gives a toroidal geometry, 0 gives a plane
            shape = strcmp(obj.fieldshape,'torus');
            % Set params
            % First dimension is time
            epoch = size(obj.field,1);
            sy = size(obj.field,2);
            sx = size(obj.field,3);
            
            hood = obj.nhood;
            % Flip all dimensions of nhood in prep for convolution
            for n = 1:nstates
                for d = 1:numel(size(hood{n}));
                   hood{n} = flipdim(hood{n},d);
                end
            end
            
            % Determine the size of the pad necessary for edges by finding the
            % max non-time dimension of the neighborhoods
            % Also determine maximum time dimension of nhoods
            
            dists = [];
            mems = [];
            for p = 1 : nstates
                temp = size(hood{p});
                mems(end+1) = temp(1);
                dists(end+1) = max(temp(2:end));
            end
            sensdist = (max(dists)-1)/2;
            memdist = max(mems);
                        
            % generate all the subscript coordinates
            ncells = numel(obj.field(1,:,:));
            [allsubs(:,1), allsubs(:,2)] = ind2sub([sy sx],1:ncells);

            % Translate subscript indices from field into their corresponding linear
            % indices in subfield
            allinds = sub2ind([sy + 2*sensdist sx + 2*sensdist], allsubs(:,1) + sensdist, allsubs(:,2) + sensdist);
            
            % Figure out where to start based on where the first 0 is in field
            init = min(1+floor(find(permute(obj.field, [3 2 1]) == 0)/ncells));
            if isempty(init)
                disp('No init value (e.g., 0) found in field. Cannot proceed.')
                return
            end
            
%             % Prepare figure
%             if obj.show == 1;
%                 fig = figure('color','k','position',[10 10 700 700]);
%                 imge = imagesc(squeeze(obj.field(1,:,:)));
%                 set(gca,'units','normalized')
%                 set(gca,'position',[0 0 1 1])
%                 axis image
%                 axis off
%             end
%             
            % Make neighborhood indices
            % thanks to mattj on matlab help for this
            % Dimensions of the padded field.
            %             sz = [sy sx] + 2 * sensdist;
            %             % Distance from central cell to edge of nhood
            %             w = sensdist + 1;
            %             % Width of nhood
            %             ww = w + sensdist;
            %             % Form an (x,y) grid with the size of nhood
            %             [ii ,jj] = ndgrid(1:ww);
            % Convert the grid to linear indices that can be added to an
            % index to find the linear indices of all the cells around it
            % TODO (DVB) assign the neighborhood at this point by paring
            % down jumps
            %             jumps = sub2ind(sz, ii, jj) - sub2ind(sz, w, w);
            %             jumps = jumps(:);
            %             nbr_inds = bsxfun(@plus,inds,jumps(:).');
            
            % TODO (DVB) What is the best way to handle neighborhoods if
            % the neighborhood can be 3-dimensional (past-dependent) and
            % state-dependent? Adding dimensions makes ultimately linearizing it
            % increasingly difficult. Let's forget about linearizing it for
            % now, since that doesn't scale well at all.
            nextstates = [2:nstates 1];            
                     
            % divide the image into chunks and operate on those in
            % parallel
            divs = [];
            parInds = [];
            nWorkers = matlabpool('SIZE');
            if nWorkers > 1
                divs = floor(sy/nWorkers);
                for i = 1:nWorkers
                    parInds(i,:) = [1, divs] + ((i-1) * divs);
                end
                parInds(end,end) = sy;
                
            end
            
            for i = init:(epoch)
                % Pad the field by giving it a toroidal geometry or
                % embedding the entire array in 0s, where the width of the
                % padding is determined by the largest neighborhood matrix.
                nextfield = squeeze(obj.field(i-1,:,:));
                if shape == 1
                    subfield = maketorus(squeeze(obj.field(i-memdist:i,:,:)),sensdist);
                else
                    subfield = zeros(memdist,sy+2*sensdist,sx+2*sensdist);
                    subfield(:,sensdist+1:end-sensdist,sensdist+1:end-sensdist) = obj.field(i-memdist:i-1,:,:);
                end
                                
                matches = [];
                
                for n = 1:nWorkers
                    locInds = parInds(n,:);
                    parIm{n} = squeeze(subfield(:,sensdist+(locInds(1):locInds(2)),:));
                    
                end
                
                parfor n = 1:nWorkers
                    locInds = parInds(n,:);
                    locIm = squeeze(obj.field(i-1,locInds(1):locInds(2),:));
                    locSubIm = subfield(:,sensdist+(locInds(1):locInds(2)));
                    for k = 1:nstates
                        loc_hood = hood{k};
                        loc_go = obj.go{k};
                        convolved = (locIm == k).*squeeze(convn(locSubIm == nextstates(k),loc_hood,'valid'));
                        matches{k} = find(sum(bsxfun(@eq, convolved(:), loc_go),2));
                        nextfield(matches) = nextstates(k);
                    end
                end
                
%                      for k = 1:nstates
%                     loc_hood = hood{k};
%                     loc_go = obj.go{k};
%                     convolved = (squeeze(obj.field(i-1,:,:)) == k).*squeeze(convn(subfield == nextstates(k),loc_hood,'valid'));
%                     
%                     matches = find(sum(bsxfun(@eq, convolved(:), loc_go),2));
%                     nextfield(matches) = nextstates(k);
%                 end
%                 
                % Put data in obj.field
                obj.field(i,:,:) = nextfield;
                % Plot it
                if obj.show == 1;
%                     colormap(obj.colorsc.cmap);
                    imagesc(squeeze(obj.field(i,:,:)));
                    drawnow
                end
                
                if i > obj.states && isequal(obj.field(i-obj.states,:,:),obj.field(i,:,:)) || isequal(obj.field(i-1,:,:),obj.field(i,:,:))
                    % Wipe the redundant future times
                    obj.field(i+1:end,:,:) = [];
                    obj.behavior.froze = 1;
                    disp('Abort: Simulation oscillating or crystallized')
                    break
                end
               
                %
                %                 if i > 2
                %                     % if i == 2 then there's nothing to check before updating.
                %                     % otherwise, we can subtract the previous two fields and
                %                     % only operate on the set of all cells within sendist of a
                %                     % cell that changed on the last update.
                %                     % grab the indices of all the cells that changed
                %                     deltas = find(squeeze(subfield(end,:,:)) - maketorus(squeeze(obj.field(i-2,:,:)),sensdist));
                %                     % add indices of their neighbors
                %                     delta_nbrs = bsxfun(@plus,deltas,jumps(:).');
                %                     delta_nbrs(delta_nbrs < 1 | delta_nbrs > numel(subfield)) = [];
                %                     % we only want the unique values falling within allinds
                %                     % convert to subscripts
                %                     % TODO (DVB) check whether this calculation is
                %                     % redundant
                %                     deltas = unique(vertcat(deltas(:),  delta_nbrs(:)));
                %                     [dy, dx] = ind2sub(size(subfield),deltas);
                %                     % toss all cells outside the field
                %                     deltas(dy' <= sensdist | dy' > (sy + sensdist) | dx' <= sensdist | dx' > (sx + sensdist)) = [];
                %                     inds = deltas;
                %                     % add to this list all the cells in the neighborhood of
                %                     % these cells
                %                 else
                %                     inds = allinds;
                %                 end
                
                % pre-calculate weighted neighborhoods
                % thanks to mattj on matlab help for these
                % TODO (DVB) implement looping through different go
                % conditions
                
                %                 NeighborTable = nbr_inds;
                %                 for n = 1:size(nhood,1);
                %                     NeighborTable(n,:) = bsxfun(@times,nbr_inds, squeeze(nhood(n,:)).');
                %                     NeighborTable(n,:) = sum(bsxfun(@eq, NeighborTable(n,:), nextstates(inds)),2);
                %                 end
                %
                %                 for k = 1:numel(obj.go)
                %                     matches = find(sum(bsxfun(@eq, NeighborTable, obj.go{k}),2));
                %                     subfield(inds(matches)) = nextstates(inds(matches));
                %                 end
                % TODO (DVB) test efficacy of using convn instead of
                % manually applying the neighborhood filter. this requires
                % first flipping the neighborhood across all dimensions but
                % convn should be really really fast.
                %                 NeighborTable = nbr_inds;

            end
            
            if obj.show == 1
                disp(['Simulation completed in ' num2str(toc)  's']);
            end
            
            function torus = maketorus(mat,distance)
                % given an array mat with size [z,y,x], return an array of
                % size = [z,y + 2 * distance, x + 2 * distance] where the last [distance] rows
                % & columns have been copied to the opposite side of the array,
                % respectively. So given mat = [1 2 3; 4 5 6; 7 8 9] and
                % distance = 1, maketorus(mat,distance) returns
                % [9 7 8 9 7; 3 1 2 3 1; 6 4 5 6 4; 9 7 8 9 7; 3 1 2 3 1];
                torus = zeros(size(mat,1), size(mat,2)+ 2*distance, size(mat,3) + 2*distance);
                % Bottom right
                torus(:,(2*distance+1):end,(2*distance+1):end) = circshift(mat,[0,-distance,-distance]);
                % Top left
                torus(:,1:end-(2*distance),1:end-(2*distance)) = circshift(mat,[0,distance,distance]);
                % Bottom left
                torus(:,(2*distance+1):end,1:end-(2*distance)) = circshift(mat,[0,-distance,distance]);
                % Top right
                torus(:,1:end-(2*distance),(2*distance+1):end) = circshift(mat,[0,distance,-distance]);
            end
        end
        
        % increase the length and size of a simulation
        function obj = extend(obj, ext)
            if length(ext) > 3
                error('Extension must be a vector with three or fewer elements')
            elseif length(ext) < 3
                % stuff ext with zeros if it has fewer than 3 elements
                ext(3) = 0;
            elseif ext(1) < 1
                error('Must add positive, non-zero amount of time')
            elseif ext(2) < 2 || ext(3) < 2
                ext(ext < 2) = 0;
                warning('Lengths < 2 are ignored')
            end
            % extension is a vector with three values
            % extension(1) = time to add
            % extension([2 3]) = x and y coords to add
            % for now extending area is padded with ones. if
            % diff(extension(2),size(field,2)) is odd, then extension(2) =
            % extension(2) - 1
            ep = size(obj.field,1);
            sy = size(obj.field,2);
            sx = size(obj.field,3);
            
            newfield = ones(ep,sy+ext(2),sx + ext(3));
            newfield(:,1 + ext(2) / 2 : sy + ext(2) / 2,...
                1 + ext(3) / 2 : sx + ext(3) / 2) = obj.field;
            newfield(end+1:end+ext(1),:,:) = zeros(ext(1),...
                sy + ext(2),sx + ext(3));
            obj.field = newfield;
        end
        
        function playBack(obj,varargin)
            % play a movie of the simulation by calling imagesc repeatedly
            % on obj.field
            % args in: a caut object you wish to see animated
            % rate: time between frames. default is 24 fps, i.e. ~.05s
            % refresh
            % colors: a colormap for the image. default is obj.colorsc
            % reps: number of times to loop the video. default is 1
            % fig: a handle to a figure for plotting. default behavior
            % generates a new figure for each playBack call
            par = inputParser;
            par.addRequired('obj')
            par.addOptional('fig',1,@ishandle)
            par.addOptional('rate',24,@(x) x > 0)
            par.addOptional('colors',obj.colorsc.cmap,@(x) ischar(x) || (ismatrix(x) && size(x,2) == 3))
            par.addOptional('reps',1,@(x) x > 0)
            par.parse(obj,varargin{:});
            
            fig = par.Results.fig;
            colors = par.Results.colors;
            rate = 1/par.Results.rate;
            reps = par.Results.reps;
            
            if ishandle(fig)
                set(0,'currentfigure',fig);
            else
                figure(fig)
            end
            imge = imagesc(squeeze(obj.field(1,:,:)));
            colormap(colors)
            axis off
            axis image
            axis manual
            pause(rate)
            for k = 1:reps
                for i = 1:size(obj.field,1)
                    set(imge,'cdata',squeeze(obj.field(i,:,:)))
                    colormap(colors)
                    pause(rate)
                end
            end
        end
        
        function vectout = qualCheck(obj)
        % apply various measures to a completed simulation and return a
        % vector of booleans specifying which checks failed
        field = obj.field;
        % first check the distribution of states in the final frame
        
        
        
        end
        
        function out = makeGif(obj, outfname)
            % save sim as a .gif
            % check outpath for validity
            if exist('outfname','var')
                sepind = strfind(outfname,filesep);
                pathname = outfname(1:sepind(end));
                filename = outfname(sepind(end)+1:end);
            else
                filename = [obj.simid '.gif'];
                pathname = [pwd filesep];
            end
            toplot = uint8(obj.field);
            
            screensize = get(0,'ScreenSize');
            figpos(1) = 10;
            figpos(2) = 10;
            figpos(3) = 2.5*size(toplot,3);
            figpos(4) = 2.5*size(toplot,2);
            
            if ~isdir(pathname)
                disp([outfname ' does not contain a valid path. Using pwd.'])
                pathname = [];
            end
            fig = figure('color','k','position',figpos);
            set(0,'currentfigure',fig);
            imge = imagesc(squeeze(toplot(end,:,:)));
            axis off
            axis image
            drawnow
            fframes = 100;
            colormap(obj.colorsc.cmap)
            if ~isempty(obj.colorsc.caxis)
                caxis(obj.colorsc.caxis);
            end
            
            for n = [ -1 * ones(1,fframes) 1:size(toplot,1)]
                % throw up the last frame so the thumbnail works
                if  n == -1
                    % check if figure still exists
                    if ishandle(fig)
                        set(imge,'cdata',squeeze(toplot(end,:,:)));
                        drawnow
                    else
                        fig = figure('color','k','position',figpos);
                        imge = imagesc(squeeze(toplot(end,:,:)));
                        axis off
                        axis image
                        drawnow
                    end
                else
                    % check if figure still exists
                    if ishandle(fig)
                        set(imge,'cdata',squeeze(toplot(n,:,:)));
                        drawnow
                    else
                        fig = figure('color','k','position',figpos);
                        imge = imagesc(squeeze(toplot(n,:,:)));
                        axis off
                        axis image
                        drawnow
                    end
                end
                frame = getframe(1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if n == -1;
                    imwrite(imind,cm,[pathname filename],'gif', 'Loopcount',inf,'delaytime',0);
                else
                    imwrite(imind,cm,[pathname filename],'gif','WriteMode','append','delaytime',0);
                end
            end
            close
            out = 1;
        end
        
        % snip out a section of the simulation. returns a (rnge,sidey,sidex) array.
        function out = snap(obj,rnge)
            par = inputParser;
            par.addRequired('obj')
            par.addOptional('rnge',size(obj.field,1),@(x) true)
            par.parse(obj,rnge);
            rnge = par.Results.rnge;
            out = double(squeeze(obj.field(rnge,:,:)));
        end
    end
end

