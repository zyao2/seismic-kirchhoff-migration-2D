close all;
clear;
clc;

%% Read in velocity model data

fname='inputParameters.txt';
nx=readParam(fname,'nx');
nz=readParam(fname,'nz');
gridSize=readParam(fname,'gridSize');
dx=gridSize;dz=gridSize;
% Define frequency parameter for ricker wavelet
f = readParam(fname,'f');
shotWatchInterval=readParam(fname,'shotWatchInterval');
shotInterval=readParam(fname,'shotInterval');
shotGatherFileDir=readParam(fname,'shotGatherFileDir');
if exist(shotGatherFileDir, 'dir')
    rmdir(shotGatherFileDir,'s');
    mkdir(shotGatherFileDir);
else
    mkdir(shotGatherFileDir);
end

 velocityModel=vmodel(nz,nx);
nBoundary = 20;
x = (0:nx) * dx;
z = (0:nz) * dz;
% grids and positions of shot array
xShotGrid=1:shotInterval:nx;
nShots = length(xShotGrid);
xShot = (xShotGrid-1) * dx;

xRecGrid=1:nx;
xRec = (xRecGrid-1) * dx;


    figure;
    set(gcf,'Position',[100 100 1200 400])

    subplot(1,3,1);
    imagesc(x, z, velocityModel)
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Velocity Model');
    hold on;
    hShotPos = plot(xShot(1), z(1), 'w*');
    hold off;
    colormap(seismic);

    %% Create shot gathers
vmin = min(velocityModel(:));
vmax = max(velocityModel(:));
dt = 0.75*(dz/vmax/sqrt(2));
nt = round(sqrt((dx*nx)^2 + (dz*nz)^2)*2/vmin/dt + 1);
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions
V = [repmat(velocityModel(:, 1), 1, nBoundary),velocityModel, repmat(velocityModel(:, end), 1, nBoundary)];
% number of approximation order for differentiator operator
nDiffOrder = 3;
% generate shot signal
Ssignal=ricker(f, nt, dt);
for ixs = 1:nShots    
    xs = xShotGrid(ixs); % shot position on x    
    % generate shot signal
    source = zeros([size(V), nt]);
    % source(1, xs, 1) = 1; % impulse input
    source(1, xs + nBoundary, :) = reshape(Ssignal, 1, 1, nt);  
    % plot shot position
    set(hShotPos, 'XData', x(xs), 'YData', z(1));  
    % generate shot record
    tic;
    %[data, snapshot] = fwd2d(V, source, nDiffOrder, nBoundary, dz, dx, dt);
    [data, snapshot] = fwd2d_mex(V, source, nDiffOrder, nBoundary, dz, dx, dt);
   % data=data(KK:end,:);
    timeForward = toc;  
    data = data(nBoundary+1:end-nBoundary,:)'; 
    shotGatherFilename = [shotGatherFileDir, sprintf('/shotGather%d.mat', ixs)];
    save(shotGatherFilename, 'data', '-v7.3');
    
   
        if(ixs==floor(ixs/shotWatchInterval)*shotWatchInterval+1)
            start_t = 1;
        else
            start_t = nt;
        end

        for it = start_t:nt        
            % plot shot record evolution (true)
            ds = zeros(nt, nx);
            ds(1:it, :) = data(1:it, :);
            subplot(1,3,2)
            imagesc(x, t, ds)
            xlabel('Distance (m)'), ylabel('Time (s)')
            title(sprintf('Shot No. %d at x = %dm', ixs, x(xs)));
            caxis([-0.1 0.1])       
            % plot wave propagation (true)
            subplot(1,3,3)
            imagesc(x, z, snapshot(1:end-nBoundary, nBoundary+1:end-nBoundary, it))
            xlabel('Distance (m)'), ylabel('Depth (m)')
            title(sprintf('Wave Propagation (True) t = %.3f', t(it)));
            caxis([-0.14 1])
            colormap(seismic);
            drawnow;
        end
    
end


