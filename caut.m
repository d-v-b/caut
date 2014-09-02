
%%%%%%%%%%%%%%%%%%%%%%%
% Cellular automaton class
% --Davis Bennett--
%--  03/19/2013  --
%%%%%%%%%%%%%%%%%%%%%%%%
% todo(dvb) replace convolution with fft2, or at least try it
% todo(dvb) make field dimensions [space,space,time] instead of
% [time,space,space] in order to reduce the number of squeezes necessary
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
        qualParams % struct containing parameters for the quality checks
        behavior % struct containing boolean-valued fields that describe the behavior of the sim
        smooth % integer representing the number of interpolation frames 
        % to insert between simulation frames
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
            % Initialize the caut object
            par = inputParser;
            par.addRequired('field')
            par.addRequired('nhood')
            par.addRequired('go')
            par.addParamValue('nstates',-1+numel(unique(fld_matrix(:,:,:))))
            par.addParamValue('show',1)
            par.addParamValue('fieldshape','plane')
            par.addParamValue('smooth',0)
            par.parse(fld_matrix, nbr_matrix, go_matrix,varargin{:})
            
            obj.field = par.Results.field;
            obj.nhood = par.Results.nhood;
            obj.go = par.Results.go;
            obj.states = par.Results.nstates;
            obj.smooth = par.Results.smooth;
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
                end
            end
            
            if ~isequal(numel(obj.nhood), obj.states)
                if numel(obj.nhood) == 1
                    obj.nhood(2:nstates) = obj.nhood(1);
                else
                    obj.states = numel(obj.nhood);
                end
            end
            
            obj.show = par.Results.show;
            obj.fieldshape = par.Results.fieldshape;
        end
        
        % function to run the simulation
        function obj =  runSim(obj)
            runtime = tic;
            tic;
            nstates = obj.states;
            % 1 gives a toroidal geometry, 0 gives a plane
            shape = strcmp(obj.fieldshape,'torus');
            % Set params
            % First dimension is time
            epoch = size(obj.field,3);
            sy = size(obj.field,1);
            sx = size(obj.field,2);
            % Need to track the time between frames for proper smoothing
            dt = [];
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
            mems = 1;            
            for p = 1 : nstates
                temp = size(hood{p});
                if numel(temp) == 3 
                    mems(end+1) = temp(3);
                else
                    mems(end+1) = 1;
                end
                dists(end+1) = max(temp(1:2));
            end
            sensdist = (max(dists)-1)/2;
            memdist = max(mems);

            % Figure out where to start based on where the first 0 is in field           
            ncells = numel(obj.field(:,:,1));
            init = min(1+floor(find(obj.field == 0)/ncells));
            if isempty(init)
                disp('No init value (e.g., 0) found in field. Cannot proceed.')
                return
            end
                                             
            for i = init:(epoch)
                
                % Pad the field by giving it a toroidal geometry or
                % embedding the entire array in 0s, where the width of the
                % padding is determined by the largest neighborhood matrix.
                
                dt(1) = toc;
                
                nextfield = squeeze(obj.field(:,:,i-1));
                if shape == 1
                    subfield = padarray(squeeze(obj.field(:,:,i-memdist:i-1)),[sensdist,sensdist,0],'circular');
                else
                    subfield = zeros(sy+2*sensdist,sx+2*sensdist,memdist);
                    subfield(sensdist+1:end-sensdist,sensdist+1:end-sensdist,:) = obj.field(:,:,i-memdist:i-1);
                end
                
                active = unique(subfield);
                
                for m = 1:numel(active)                    
                    k = active(m);
                    %convolved = (squeeze(obj.field(:,:,i-1)) == k).*squeeze(convn(subfield == 1+mod(k,nstates),hood{k},'valid'));
                    convolved = (squeeze(obj.field(:,:,i-1)) == k).*squeeze(convn(subfield == k,hood{k},'valid'));
                    matches = find(sum(bsxfun(@eq, convolved(:), obj.go{k}),2));
                    nextfield(matches) = 1+mod(k,nstates);
                end
                
                % Put data in obj.field
                obj.field(:,:,i) = nextfield;                               
                
                % todo(dvb) optimize rejection criteria for speed
                if i > obj.states && isequal(obj.field(:,:,i-obj.states),obj.field(:,:,i)) || isequal(obj.field(:,:,i),obj.field(:,:,i-1))
                    % Wipe the redundant future times
                    obj.field(:,:,i+1:end) = [];
                    disp('Abort: Simulation oscillating or crystallized')
                    break
                end
                
                % Store the elapsed time for this frame in dt(2)
                dt(2) = abs(toc - dt(1));
                
                % plot it
                if obj.show == 1
                    
                    if obj.smooth > 1
                    
                   	smoothMat = obj.field(:,:,i) - obj.field(:,:,i-1);
                    smoothMat = smoothMat./(obj.smooth - 1);
                    for s = 1:(obj.smooth - 1)
                        imagesc(obj.field(:,:,i-1) + smoothMat * s);
                        drawnow                        
                        pause(dt(2))
                    end
                    
                    else
                        imagesc(obj.field(:,:,i))
                        drawnow                        
                    end
                end
                
            end
            
            disp(['Simulation completed in ' num2str(toc(runtime))  's']);    
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
            imge = imagesc(squeeze(obj.field(:,:,1)));
            colormap(colors)
            axis off
            axis image
            axis manual
            pause(rate)
            for k = 1:reps
                for i = 1:size(obj.field,1)
                    set(imge,'cdata',obj.field(:,:,i))
                    colormap(colors)
                    pause(rate)
                end
            end
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
            imge = imagesc(squeeze(toplot(:,:,end)));
            axis off
            axis image
            drawnow
            fframes = 100;
            colormap(obj.colorsc.cmap)
            if ~isempty(obj.colorsc.caxis)
                caxis(obj.colorsc.caxis);
            end
            
            for n = [ -1 * ones(1,fframes) 1:size(toplot,3)]
                % throw up the last frame so the thumbnail works
                if  n == -1
                    % check if figure still exists
                    if ishandle(fig)
                        set(imge,'cdata',squeeze(toplot(:,:,end)));
                        drawnow
                    else
                        fig = figure('color','k','position',figpos);
                        imge = imagesc(squeeze(toplot(:,:,end)));
                        axis off
                        axis image
                        drawnow
                    end
                else
                    % check if figure still exists
                    if ishandle(fig)
                        set(imge,'cdata',squeeze(toplot(:,:,n)));
                        drawnow
                    else
                        fig = figure('color','k','position',figpos);
                        imge = imagesc(squeeze(toplot(:,:,n)));
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
            par.addOptional('rnge',size(obj.field,3),@(x) true)
            par.parse(obj,rnge);
            rnge = par.Results.rnge;
            out = double(squeeze(obj.field(:,:,rnge)));
        end
    end
end

