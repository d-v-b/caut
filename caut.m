
%%%%%%%%%%%%%%%%%%%%%%%
% Cellular automaton class
% Properties:
% Methods:

% --Davis Bennett--
%--  03/19/2013  --

%%%%%%%%%%%%%%%%%%%%%%%%


classdef caut

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
        
        function obj = caut(field, nbrs, ingo, nstates ,toshow, colors)
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
            % generate all the subscript coordinates
            [allsubs(:,1) allsubs(:,2)] = ind2sub([sy sx],1:numel(obj.field(1,:,:)));
            ncells = numel(obj.field(1,:,:));
            % translate subscript indices from field into their corresponding linear 
            % indices in subfield
            allinds = sub2ind([sy + 2*sensdist sx + 2*sensdist], allsubs(:,1) + sensdist, allsubs(:,2) + sensdist);
            % Figure out where to start based on where the first 0 is in field
            init = min(1+floor(find(permute(obj.field, [3 2 1]) == 0)/ncells));
            if isempty(init)
              disp('No init value (e.g., 0) found in field. Cannot proceed.')
              return
            end
            
            % make neighborhood indices
            % thanks to mattj on matlab help for this
            sz=[sy sx] + 2*sensdist;
            w=sensdist+1;
            ww=w+sensdist;
            [ii,jj]=ndgrid(1:ww);
            jumps= sub2ind(sz,ii,jj) - sub2ind(sz, w, w) ;
            jumps=jumps(:);
            nhood=obj.nhood(:);
            
            if obj.show == 1;
                fig = figure('color','k','position',[10 10 700 700]);
                imge = imagesc(squeeze(obj.field(1,:,:)));
                set(gca,'units','normalized')
                set(gca,'position',[0 0 1 1])
                if ~isempty(obj.colorsc.caxis)
                    caxis(obj.colorsc.caxis);
                end
                axis image
                axis off
            end
             
            for i = init:(epoch)
                % To do: accelerate simulation by tracking which cells
                % chnaged and only checking the smallest area containing those 
                % changed cells, plus sensdist.
                % we can check by storing which cells changed using linear
                % coordinates, then checking which indices are closest to
                % each axis and making a box with those edges
                
                % make a torus
                subfield = maketorus(squeeze(obj.field(i-1,:,:)),sensdist);
                % pre-calculate all transitions
                nextstates = subfield(:) + 1;
                % values above nstates are ignored
                nextstates(nextstates ==  nstates+1) = 1;
                
                if i > 2
                    % if i == 2 then there's nothing to check before updating.
                    % otherwise, we can subtract the previous two fields and 
                    % only operate on the set of all cells within sendist of a
                    % cell that changed on the last update. 
                    % grab the indices of all the cells that changed
                    deltas = find(subfield - maketorus(squeeze(obj.field(i-2,:,:)),sensdist));
                    % add indices of their neighbors
                    delta_nbrs = bsxfun(@plus,deltas,jumps(:).');
                    delta_nbrs(delta_nbrs < 1 | delta_nbrs > numel(subfield)) = [];
                    % we only want the unique values falling within allinds
                    % convert to subscripts (possibly inline, if this is
                    % slow)
                    deltas = unique(vertcat(deltas(:),  delta_nbrs(:)));
                    [dy dx] = ind2sub(size(subfield),deltas);
                    % toss all cells outside the field
                    deltas(dy' <= sensdist | dy' > (sy + sensdist) | dx' <= sensdist | dx' > (sx + sensdist)) = [];
                    inds = deltas;
                    % add to this list all the cells in the neighborhood of
                    % these cells
                else
                    inds = allinds;
                end
                
                % pre-calculate weighted neighborhoods
                % thanks to mattj on matlab help for these
                NeighborTable=subfield(bsxfun(@plus,inds,jumps(:).'));
                NeighborTable=bsxfun(@times,NeighborTable, nhood(:).');
                % see how vectorizing speeds things up
                NeighborTable = sum(bsxfun(@eq, NeighborTable, nextstates(inds)),2);
                NeighborTable = find(sum(bsxfun(@eq, NeighborTable, obj.go),2));
                subfield(inds(NeighborTable)) = nextstates(inds(NeighborTable));

                obj.field(i,:,:) = subfield(1+sensdist:(end-sensdist),1+sensdist:(end-sensdist));
                if obj.show == 1;
                    colormap(obj.colorsc.cmap);
                    set(imge,'cdata',squeeze(obj.field(i,:,:)));
                    drawnow
                end
                
                if i > obj.states && isequal(obj.field(i-obj.states,:,:),obj.field(i,:,:)) || isequal(obj.field(i-1,:,:),obj.field(i,:,:))
                    % Wipe the redundant future times
                    obj.field(i+1:end,:,:) = [];
                    disp('Abort: Simulation oscillating or crystallized')
                    break
                end
            end
            
            if obj.show == 1
                disp(['Simulation completed in ' num2str(toc)  's']);
            end
            function torus = maketorus(mat,distance)
            % given an array of size = [y,x], return an array of 
            % size = [y,x] + 2 * distance where the last [distance] rows
            % columns have been copied to the opposite side of the array,
            % respectively. So given mat = [1 2 3; 4 5 6; 7 8 9] and
            % distance = 1, maketorus(mat,distance) returns
            % [9 7 8 9 7; 3 1 2 3 1; 6 4 5 6 4; 9 7 8 9 7; 3 1 2 3 1];
            torus = zeros(size(squeeze(mat))+2*distance);
            % Bottom right
            torus((2*distance+1):end,(2*distance+1):end) = circshift(squeeze(mat),[-distance,-distance]);
            % Top left
            torus(1:end-(2*distance),1:end-(2*distance)) = circshift(squeeze(mat),[distance,distance]);
            % Bottom left
            torus((2*distance+1):end,1:end-(2*distance)) = circshift(squeeze(mat),[-distance,distance]);
            % Top right
            torus(1:end-(2*distance),(2*distance+1):end) = circshift(squeeze(mat),[distance,-distance]);
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

