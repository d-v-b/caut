classdef caut
    % Cellular automaton, cyclic for now
    properties
        states % int: states for a cell
        nhood % array: neighbors around a cell. a matrix with a NaN at the cell
        % to be updated; values at all other positions determine the weight of
        % that position in the updating step
        go % array: conditions for a cell to advance
        stay % array: conditions for a cell to stay put
        field % array: the simulation itself
        show % struct: determines how simulation is displayed
        dump % struct: params for saving simulation
        colorsc % colorscale for visualizing system
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
            
            field = '-f';
            for i = 1:3
                if i == 1
                    field = strcat(field, [num2str(size(obj.field,i))]);
                else
                    field = strcat(field, ['-' num2str(size(obj.field,i))]);
                end
            end
            id = [stat go hood field];
        end
        
        function obj = caut(field, nbrs, nstates, ingo,toshow,colors)
            % Constructs the intial condition of the simulation
            par = inputParser;
            defcol.cmap = 'gray';
            defcol.caxis = [];
            par.addRequired('field',@(x) numel(size(x)) > 2)
            par.addRequired('nbrs',@(x) size(x,1) == size(x,2))
            par.addRequired('ingo')
            par.addOptional('nstates',max(unique(field(1,:,:))),@(x) x > 1)
            par.addOptional('toshow',1,@(x) x == 1 || x == 0)
            par.addOptional('colors',defcol,@(x) isstruct(x))
            
            par.parse(field,nbrs,ingo,nstates,toshow,colors)
            
            obj.field = par.Results.field;
            obj.nhood = par.Results.nbrs;
            obj.go = par.Results.ingo;
            obj.states = par.Results.nstates;
            obj.show = par.Results.toshow;
            obj.colorsc = par.Results.colors;
        end
        
        % function to run the simulation
        function obj =  runSim(obj)
            tic
            nstates = obj.states;
            % determine the size of the pad necessary for toroidal geometry
            sensdist = (size(obj.nhood,1)-1)/2;
            % set params
            epoch = size(obj.field,1);
            sy = size(obj.field,2);
            sx = size(obj.field,3);
            [allsubs(:,1) allsubs(:,2)] = ind2sub([sy sx],1:numel(obj.field(1,:,:)));
            ncells = numel(obj.field(1,:,:));
            trans = circshift(1:nstates,[0 -1]);
            allinds = sub2ind([sy + 2*sensdist sx + 2*sensdist], allsubs(:,1) + sensdist, allsubs(:,2) + sensdist);
            % Figure out where to start based on where the first 0 is in field
            init = min(1+floor(find(permute(obj.field, [3 2 1]) == 0)/ncells));
            if isempty(init)
              disp('No init value found in field.')
              return
            end
            % make neighborhood indices
            sz=[sy sx] + 2*sensdist;
            w=sensdist+1;
            ww=w+sensdist;
            [ii,jj]=ndgrid(1:ww);
            jumps= sub2ind(sz,ii,jj) - sub2ind(sz, w, w) ;
            jumps=jumps(:);
            nhood=obj.nhood(:);
            
            if obj.show == 1;
                fig = figure('color','k');
                imge = imagesc(squeeze(obj.field(1,:,:)));
                if ~isempty(obj.colorsc.caxis)
                    caxis(obj.colorsc.caxis);
                end
                axis image
                axis off
            end
            
            for i = init:(epoch)
                % make a torus
                subfield = zeros(size(squeeze(obj.field(1,:,:)))+2*sensdist);
                % Bottom right
                subfield((2*sensdist+1):end,(2*sensdist+1):end) = circshift(squeeze(obj.field(i-1,:,:)),[-sensdist,-sensdist]);
                % Top left
                subfield(1:end-(2*sensdist),1:end-(2*sensdist)) = circshift(squeeze(obj.field(i-1,:,:)),[sensdist,sensdist]);
                % Bottom left
                subfield((2*sensdist+1):end,1:end-(2*sensdist)) = circshift(squeeze(obj.field(i-1,:,:)),[-sensdist,sensdist]);
                % Top right
                subfield(1:end-(2*sensdist),(2*sensdist+1):end) = circshift(squeeze(obj.field(i-1,:,:)),[sensdist,-sensdist]);
                % to do: check which area(s) need to be updated
                temp = squeeze(obj.field(i-1,:,:));
                % pre-calculate all transitions
                nextstates = temp(:) + 1;
                % values above nstates are ignored
                nextstates(nextstates ==  nstates+1) = 1;
                % pre-calculate weighted neighborhoods
                % thanks to mattj on matlab help for these
                NeighborTable=subfield(bsxfun(@plus,allinds,jumps(:).'));
                NeighborTable=bsxfun(@times,NeighborTable, nhood(:).');
                % loop through cells
                
                for k = 1:ncells
                    neighbors = NeighborTable(k,:);
                    quorum = sum(neighbors == nextstates(k));
                    if ~isempty(find(obj.go == quorum,1))
                        temp(k) = nextstates(k);
                    end
                end
                
                obj.field(i,:,:) = temp;
                if obj.show == 1;
                    colormap(obj.colorsc.cmap);
                    set(imge,'cdata',squeeze(obj.field(i,:,:)));
                    drawnow
                end
                
                if i > obj.states && isequal(obj.field(i-obj.states,:,:),obj.field(i,:,:)) || isequal(obj.field(i-1,:,:),obj.field(i,:,:))
                    % Wipe the redundant future times
                    obj.field(i+1:end,:,:) = [];
                    disp('Abort: Simulation oscillating')
                    break
                end
            end
            disp(['Simulation completed in ' num2str(toc)  's'])
            function z = zipintersect(x,y)
                % quickly find the intersection of two positive integer sets
                if ~isempty(x)&&~isempty(y)
                    P = zeros(1, max(max(x),max(y)) ) ;
                    P(x) = 1;
                    z = y(logical(P(y)));
                else
                    z = [];
                end
            end
            
        end
        
        % increase the length of a simulation
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
        
        % save sim as a .gif
        function out = makeGif(obj, outfname)
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

