%% CLEANUP

clear all
close all
clc

%% FOLDERS

if exist('D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'D:\GoogleDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
elseif exist('G:\My Drive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'G:\My Drive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
elseif exist('D:\GDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database','dir')
    data_folder = 'D:\GDrive\LCL\Ludovico Cademartiri\Work\projects\ARBD\database';
end
if exist('D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
elseif exist('C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='C:\Users\lcade\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
elseif exist('D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr','dir')
    output_folder='D:\OneDrive - Università degli Studi di Parma\Manuel Dedola\outlines\spherical boundary conditions for brownian dynamics simulations\submissions\JChemPhys\revision\code\freeghost_nocorr';
end
    
toolbox_folder = '..\ARBD_toolbox';
utilities_folder =  fullfile('..', '..', 'Utilities');
utilities_folder = genpath(utilities_folder);
addpath(utilities_folder)
addpath(data_folder)
addpath(toolbox_folder)
addpath(output_folder)
cmaplibrary=loadColormaps();
filenameseries='2D_%d_%d.mat';

%% FLAGS

moviemaking=true;
allpos=true;
analyze=true;
plotting=true;

%% SIMULATION PARAMETERS
S.rp=1;
S.dp=2*S.rp;

sigma=2*S.rp;
S.rc=3*sigma;
maxsteps=1e5;
reps=10;
S.stdx=S.rp/10;
nodisps=1e6;
S.timestep=1;
S.esdiff=S.stdx^2/(2*S.timestep);
S.kbT=1;

%% VARIABLES

Ns=round(logspace(log10(1e2),log10(1e3),2)',0);
phis=logspace(log10(0.5),log10(0.5),1)';
epsilons=logspace(log10(0.1),log10(1),3)';
pots=(1:2)';
bcs=(1:2)';

%% UTILITIES
uhex=[1,0;0.5,sqrt(3)/2];
[~, cmdout] = system('wmic cpu get L2CacheSize, L3CacheSize /value');
tokens = regexp(cmdout, '\d+', 'match');
S.cacheSizeMB = (max(str2double(tokens))/1024)/feature('numCores');
eps=1e-6;
COLORMAPS=loadColormaps();

%% CONDITIONS
q=1;
for i0=1:numel(Ns)
	N=Ns(i0);
	for i1=1:numel(phis)
		phi=phis(i1);
		for i2=1:numel(pots)
			pot=pots(i2);
			for i3=1:numel(bcs)
                bc=bcs(i3);
                for i4=1:numel(epsilons)
                    epsilon=epsilons(i4);
				    c(q,:)=[q,N,phi,pot,bc,epsilon];
				    q=q+1;
                end
			end
		end
	end
end
c(c(:,6)~=1 & c(:,4)~=1,:)=[];
c(:,1)=(1:size(c,1))';
%% LOOP

for ic=9:size(c,1)
    for irep=1:reps

        posfilename=sprintf(filenameseries,ic,irep);
        if ~exist(posfilename,'file')
	        S.N=c(ic,2);
	        S.phi=c(ic,3);
	        S.pot=c(ic,4);
	        S.bc=c(ic,5);
            S.epsilon=c(ic,6);
	        S.A=(S.N*pi*S.rp^2)/S.phi;
        
            % DEFINE DISPLACEMENT LIBRARY
	        BD=single(build_noise_library(S.stdx,nodisps));
	        BD(:,3)=[];
	        qd=1;
        
            % INITIALIZE POS MATRIX
            POS=single(zeros(S.N,2,maxsteps));
        
            % DEFINE FORCE AND CLAMP
            if S.pot~=0
	            H=pot_force(S.pot,S.rc,30000,sigma,S.epsilon);
                H_interpolant = griddedInterpolant(H(:,1), H(:,2), 'linear', 'nearest');
                S.pot_clamp=mcdClamp2D(nodisps, BD, S.esdiff, 1, H, 1);
            end	
        
	        % DEFINE STARTING CONFIGURATION
            if S.bc==1
                S.R=sqrt(S.A/pi);
                S.L=2*S.R;
                S.br=S.R;
                maxcoeff=ceil(S.L/S.dp);
		        p=single([]);
		        temp=(-maxcoeff:maxcoeff)';
		        tempn=numel(temp);
		        for i0=-maxcoeff:maxcoeff
			        p=[p;temp,i0.*ones(tempn,1)];
		        end
		        p=p(:,1).*uhex(1,:).*S.dp+p(:,2).*uhex(2,:).*S.dp;
                mask=vecnorm(p,2,2)>S.R-S.rp;
                p(mask,:)=[];
                if size(p,1)<S.N
                    disp('too dense')
                    pause
                end
		        p=p(randsample(size(p,1),S.N),:);
		        pnorms=vecnorm(p,2,2);
		        idxgp=pnorms>(S.br-S.rc);
		        pg=p(idxgp,:)-(2*S.br).*(p(idxgp,:)./pnorms(idxgp,:));
	        elseif S.bc==2
		        S.L=sqrt(S.A);
                S.br=0.5*S.L;
                maxcoeff=ceil(sqrt(2)*S.L/S.dp);
		        p=single([]);
		        temp=(-maxcoeff:maxcoeff)';
		        tempn=numel(temp);
		        for i0=-maxcoeff:maxcoeff
			        p=[p;temp,i0.*ones(tempn,1)];
		        end
		        p=p(:,1).*uhex(1,:).*S.dp+p(:,2).*uhex(2,:).*S.dp;
                mask=sum(p>S.L/2-S.rp | p<-S.L/2+S.rp,2)>0;
                p(mask,:)=[];
                if size(p,1)<S.N
                    disp('too dense')
                    pause
                end
		        p=p(randsample(size(p,1),S.N),:);
	        end
	        
	        % SIMULATION	
	        for it=1:maxsteps
	        
		        % collect brownian displacements
		        bdisp=BD(qd:qd+S.N-1,:);
		        qd=qd+S.N;
		        if qd+S.N>nodisps
			        BD=BD(randperm(nodisps),:);
			        qd=1;
		        end
		        % collect potential displacements
		        if S.pot~=0
			        if S.bc==1
				        potdisp = potential_displacements_2D_v13([p;pg], S, H, H_interpolant, 0, S.cacheSizeMB);
				        potdisp=potdisp(1:S.N,:);
			        elseif S.bc==2
				        potdisp = potential_displacements_2D_v13(p, S, H, H_interpolant, 0, S.cacheSizeMB);
			        end
                    bdisp=bdisp+potdisp;
		        end
		        % apply displacements
		        p=p+bdisp;
                
                % oob
                if S.bc==1
			        % 1. Check for particles that left the sphere
                    d = vecnorm(p, 2, 2);
                    idx_out = d > S.br;
                    
                    % 2. "Wrap" them to the antipodal point
                    % Instead of sorting, we simply teleport particle 'i' 
                    % from R+eps to -R+eps (reflection through origin).
                    % This preserves the index 'i' in the array.
                    if any(idx_out)
                         p(idx_out,:) = p(idx_out,:) - (2*S.br) .* (p(idx_out,:) ./ d(idx_out));
                    end
                    
                    % 3. Regenerate Ghosts (for forces in next step)
                    % Recalculate norms after the wrap
                    d = vecnorm(p, 2, 2);
                    idx_near = d > (S.br - S.rc);
                    % Generate ghosts at antipodal points
                    pg = p(idx_near,:) - (2*S.br) .* (p(idx_near,:) ./ d(idx_near));
                elseif S.bc==2
                    p=mic_wrap_positions(p, S);
                end
        
                % counter
                if mod(it,1e3)==0
                    disp([ic,log10(it)])
                end
		        
                % movie
                if mod(it,maxsteps/10)==0 & moviemaking
                    scatter(p(:,1),p(:,2),30,'filled')
                    axis equal
                    if S.bc==1
                        xlim([-S.R S.R])
                        ylim([-S.R S.R])
                    elseif S.bc==2
                        xlim([-S.br S.br])
                        ylim([-S.br S.br])
                    end            
                    drawnow
                end
        
                % position storage
                if allpos
                    POS(:,:,it)=p;
                end
            end
            save([output_folder,'\',posfilename], 'POS', 'S', 'c', '-v7.3');
        else
            load(posfilename)
            c(:,1)=(1:size(c,1))';
        end
    
        %% ANALYSIS
        if analyze
            p=POS;    
            clear POS
            dt = 1;
            p(:,:,1:floor(maxsteps/4))=[];
            [N, dim, T_steps] = size(p);

            % --- 1. COM CORRECTION & UNWRAPPING ---
            COMw = mean(p,1);   % wrapped COM
            pt=p-COMw;
            if S.bc==2
                pu = zeros(size(p));
                pu(:,:,1) = p(:,:,1);
                for t = 2:T_steps
                    dp = p(:,:,t) - p(:,:,t-1);
                    dp(dp >  S.br) = dp(dp >  S.br) - 2*S.br;
                    dp(dp < -S.br) = dp(dp < -S.br) + 2*S.br;
                    pu(:,:,t) = pu(:,:,t-1) + dp;
                end
                COMu = (mean(pu,1));
                put=pu-COMu;
            end

			% K-WINDOW
			k_fundamental = 2*pi/(2*S.br);
			k_max=pi/S.rp;
            k_mags=(k_fundamental:k_fundamental:k_max)';
            k_mags=sort([k_mags;k_fundamental.*[0.5;sqrt(2);sqrt(3);pi]]);
            % k_mags = k_fundamental.*[1/2,1,sqrt(2),sqrt(3),2,3,pi,4,5,6]'; 
            nK = length(k_mags);
			
			% THETA-WINDOW
            thetas = (0:10:170)';
            thetas_rad = deg2rad(thetas); 
            nTheta = length(thetas_rad);
            
			% OUTPUT ARRAY INITIALIZATION
            S_std = zeros(nTheta, nK);  
            S_mask = zeros(nTheta, nK); 
            Deff_std = zeros(nTheta, nK);
            Deff_mask = zeros(nTheta, nK);
            
			% DYNAMICS
			max_lag = min(floor(T_steps/3), maxsteps/10); 
            
            fprintf('Calculating Structure and Dynamics...\n');
    
            % --- 2. DEFINE MASK ---
            % Mask is defined on the COM-corrected 'p'
            if S.bc == 2 % PBC
				R_mask = S.L/2;
				dist_sq = sum(p.^2, 2); 
				mask = squeeze(dist_sq < R_mask^2); 
				N_eff = mean(sum(mask, 1));				
            else % SBC
                mask = true(N, 1, T_steps);
                N_eff = N;
            end
    
            % --- 3. LOOP K-VECTORS ---
            for th_i = 1:nTheta
				
				% initialize angles and utilities
                theta = thetas_rad(th_i);
                cos_th = cos(theta); sin_th = sin(theta);
                
                for k_i = 1:nK
				
					% initialize k-vectors and utilities
                    k_val = k_mags(k_i);
                    qx = k_val * cos_th; qy = k_val * sin_th;
                    
                    % A. S(k)
					px = squeeze(pt(:,1,:));
                    py = squeeze(pt(:,2,:));

                    % Standard
                    rho_static = sum(exp(1i * -(qx * px + qy * py)), 1);
                    S_std(th_i, k_i) = mean(abs(rho_static).^2) / S.N;
    
                    % Masked
                    if S.bc==2
                        rho_mask = sum(mask .* exp(1i * (qx*px + qy*py)), 1);
                        S_mask(th_i, k_i) = mean(abs(rho_mask).^2) / N_eff;
                    end
                    
                    % B. DYNAMICS (Use Unwrapped)

                    
                    if S.bc==1
                        px = squeeze(pt(:,1,:)); py = squeeze(pt(:,2,:));
                        rho_dyn = sum(exp(1i * -(qx * px + qy * py)), 1);
                    elseif S.bc==2
                        pux = squeeze(put(:,1,:)); puy = squeeze(put(:,2,:));
                        rho_dyn = sum(exp(1i * -(qx * pux + qy * puy)), 1);
                        rho_dyn_mask = sum(mask .* exp(1i * (qx*pux + qy*puy)), 1);
                    end

                    Deff_std(th_i, k_i) = get_deff(rho_dyn, max_lag, dt, k_val);
                    if S.bc==2
					    Deff_mask(th_i, k_i) = get_deff(rho_dyn_mask, max_lag, dt, k_val);
                    end
                    disp([ic,irep,th_i/nTheta,k_i/nK])
                end
            end
        end
        if analyze
            SSTD(:,:,irep)=S_std;
            DEFFSTD(:,:,irep)=Deff_std;
            if S.bc==2
                SMASK(:,:,irep)=S_mask;
                DEFFMASK(:,:,irep)=Deff_mask;
            end
        end
    end
    %% PLOTTING
    if plotting
        S_std=mean(SSTD,3);
        Deff_std=mean(DEFFSTD,3);
        if S.bc==2
            S_mask=mean(SMASK,3);
            Deff_mask=mean(DEFFMASK,3);
        end
        
        figure('Color','k'); tiledlayout(2,2);
        if S.bc == 2, tit='PBC'; else, tit='SBC'; end
        X = rad2deg(thetas_rad);
        
        nexttile; imagesc(X, k_mags, S_std'); axis xy; colorbar;
        title([tit ' Standard S(k)']); xlabel('Theta'); ylabel('k');

        if S.bc==2           
            nexttile; imagesc(X, k_mags, S_mask'); axis xy; colorbar;
            title([tit ' Masked S(k)']); xlabel('Theta'); ylabel('k');
        end
    
        nexttile; imagesc(X, k_mags, Deff_std'); axis xy; colorbar;
        title([tit ' Standard Deff']); xlabel('Theta'); ylabel('k');
        
        if S.bc==2
            nexttile; imagesc(X, k_mags, Deff_mask'); axis xy; colorbar;
            title([tit ' Masked Deff']); xlabel('Theta'); ylabel('k');
        end

        figure
        plot(X,S_std(:,1))
        hold
        plot(X,S_std(:,2))
        plot(X,S_std(:,3))
        plot(X,S_std(:,4))
        plot(X,S_std(:,5))
        plot(X,S_std(:,6))
        plot(X,S_std(:,7))
        plot(X,S_std(:,8))
        plot(X,S_std(:,9))
        title([tit ' S(k,theta) for raw box']);
        if S.bc==2
            figure
            plot(X,S_mask(:,1))
            hold
            plot(X,S_mask(:,2))
            plot(X,S_mask(:,3))
            plot(X,S_mask(:,4))
            plot(X,S_mask(:,5))
            plot(X,S_mask(:,6))
            plot(X,S_mask(:,7))
            plot(X,S_mask(:,8))
            plot(X,S_mask(:,9))
            title([tit ' S(k,theta) for circular inscribed domain']);
        end
        figure
        plot(X,Deff_std(:,1))
        hold
        plot(X,Deff_std(:,2))
        plot(X,Deff_std(:,3))
        plot(X,Deff_std(:,4))
        plot(X,Deff_std(:,5))
        plot(X,Deff_std(:,6))
        plot(X,Deff_std(:,7))
        plot(X,Deff_std(:,8))
        plot(X,Deff_std(:,9))
        title([tit ' Deff(k,theta) for raw box']);
        if S.bc==2
            figure
            plot(X,Deff_mask(:,1))
            hold
            plot(X,Deff_mask(:,2))
            plot(X,Deff_mask(:,3))
            plot(X,Deff_mask(:,4))
            plot(X,Deff_mask(:,5))
            plot(X,Deff_mask(:,6))
            plot(X,Deff_mask(:,7))
            plot(X,Deff_mask(:,8))
            plot(X,Deff_mask(:,9))
            title([tit ' Deff(k,theta) for circular inscribed domain']);
        end
        figure
        scatter(p(:,1,end),p(:,2,end),50,'filled')
        axis equal
        title([tit ' final positions of particles']);
    end
end

% --- HELPER: PERIODIC CENTER OF MASS ---
function COM = get_pbc_com(pos, br)
    % Map positions to angles [-pi, pi]
    L = 2*br;
    theta = (pos + br) / L * 2 * pi; 
    
    % Average the unit vectors (circular mean)
    xi = mean(cos(theta), 1);
    zeta = mean(sin(theta), 1);
    
    % Get mean angle
    mean_theta = atan2(zeta, xi);
    
    % Map back to Cartesian [-br, br]
    COM = mean_theta / (2*pi) * L - br;
end

% --- HELPER: DEFF FIT ---
function D = get_deff(rho_t, max_lag, dt, k_val)
    if any(isnan(rho_t)), D=NaN; return; end
    rho_c = rho_t - mean(rho_t);
    acf = xcorr(rho_c, max_lag, 'biased');
    acf = acf(max_lag+1:end);
    F_norm = abs(acf) / abs(acf(1));
    time = (0:max_lag) * dt;
    idx = F_norm < 0.9 & F_norm > 0.2;
    if sum(idx) < 5
        D = NaN; 
    else
        p = polyfit(time(idx), log(F_norm(idx)), 1);
        D = -p(1) / k_val^2;
    end
end

%%
plot_polar_surface(thetas_rad, k_mags, S_std,COLORMAPS)
%%

function plot_polar_surface(thetas, ks, Deff,COLORMAPS)
    % PLOT_POLAR_SURFACE
    % Inputs:
    %   thetas: 1D array of angles in radians
    %   ks:     1D array of k magnitudes
    %   Deff:   2D matrix (nTheta x nK)
    
    % 1. SYMMETRIZE TO 360 (Optional but recommended for visual impact)
    % If your data is 0-180, we mirror it to 0-360 to show the full "Flower"
    if max(thetas) <= pi + 0.1
        fprintf('Mirroring data from 180 to 360 for visualization...\n');
        thetas_full = [thetas; thetas + pi;2*pi];
        Deff_full = [Deff; Deff;Deff(1,:)]; % Duplicate data for the second half
    else
        thetas_full = thetas;
        Deff_full = Deff;
    end
    % 2. Create the Grid (Theta on Rows, K on Cols)
    % Note: Check dimensions. 
    % If z_export is (nTheta x nK), we need grids to match.
    [K_GRID, TH_GRID] = meshgrid(ks, thetas_full);
    
    % 3. Flatten to Columns
    % Column 1: Angle in Degrees
    Col_Theta = rad2deg(TH_GRID(:));
    
    % Column 2: Radius (k)
    Col_Radius = K_GRID(:);
    
    % Column 3: Intensity (Deff)
    Col_Z = Deff_full(:);

    OriginData=[Col_Theta,Col_Radius,Col_Z];


    % 2. CREATE POLAR GRID
    % Create 2D grids. 
    % Note: meshgrid(x,y) puts x on columns and y on rows.
    [TH, R] = meshgrid(thetas_full, ks);
    
    % 3. CONVERT TO CARTESIAN (X, Y)
    % These are the coordinates for the floor of the plot
    X = R .* cos(TH);
    Y = R .* sin(TH);
    
    % 4. PREPARE Z-DATA
    % Deff_full is usually (nTheta x nK). 
    % meshgrid produced (nK x nTheta). We must TRANSPOSE Deff.
    Z = Deff_full'; 
    
    % 5. PLOT SURFACE
    figure('Color','w', 'Position', [100 100 1000 800]);
    
    % surf(X, Y, Z) creates the 3D topology
    h = surf(X, Y, Z);
    
    % 6. STYLING (The "Publication" Look)
    shading interp;          % Smooth out the grid lines
    colormap(flipud(COLORMAPS.magma));         % High contrast colormap
    colorbar;
    
    % Lighting to show 3D texture
    camlight left; 
    lighting gouraud;
    material dull; 
    
    % Axis cleanup
    axis equal;              % Make circles look circular
    axis off;                % Hide the square box axes
    view(0, 90);             % Start with Top-Down view (Heatmap style)
    
    % Add a title
    title('D_{eff}(k, \theta) Polar Topography', 'FontSize', 16);
    
    % 7. ADD REFERENCE RINGS (The "Target")
    % Draw circles at specific k values to visualize the shells
    hold on;
    for k_val = ks
        theta_ring = 0:0.01:2*pi;
        plot3(k_val*cos(theta_ring), k_val*sin(theta_ring), ...
              ones(size(theta_ring))*max(Z(:)), ...
              'k-', 'Color', [1 1 1 0.3]); % Faint white rings
    end
    
    % 8. INSTRUCTIONS
    disp('Rotate the plot using the 3D Rotate tool to see the "Mountains"!');
end