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
addpath(data_folder)
addpath(toolbox_folder)
addpath(output_folder)
filenameseries='2D_%d_%d.mat';

%% FLAGS

moviemaking=true;
allpos=true;
analyze=false;
plotting=false;

%% SIMULATION PARAMETERS
S.rp=1;
S.dp=2*S.rp;

sigma=2*S.rp;
S.rc=3*sigma;
maxsteps=1e5;
reps=1;
S.stdx=S.rp/10;
nodisps=1e6;
S.timestep=1;
S.esdiff=S.stdx^2/(2*S.timestep);
S.kbT=1;

%% VARIABLES

Ns=round(logspace(log10(1e2),log10(1e4),3)',0);
phis=logspace(log10(0.1),log10(0.5),5)';
epsilons=logspace(log10(0.1),log10(1),3)';
pots=(1:2)';
bcs=(1:2)';

%% UTILITIES
uhex=[1,0;0.5,sqrt(3)/2];
[~, cmdout] = system('wmic cpu get L2CacheSize, L3CacheSize /value');
tokens = regexp(cmdout, '\d+', 'match');
S.cacheSizeMB = (max(str2double(tokens))/1024)/feature('numCores');
eps=1e-6;

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

%% LOOP

for ic=1:size(c,1)
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
			        ptemp=[p;pg];
                    pnorms=vecnorm(ptemp,2,2);
                    [~,idxs]=sort(pnorms);
                    ptemp=ptemp(idxs,:);
                    p=ptemp(1:S.N,:);
			        pnorms=vecnorm(p,2,2);
			        idxgp=pnorms>(S.br-S.rc);
			        pg=p(idxgp,:)-(2*S.br).*(p(idxgp,:)./pnorms(idxgp,:));
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
        end
    
        %% ANALYSIS
        if analyze
            p=POS;    
            clear POS
            dt = 1;
            p(:,:,1:maxsteps/2)=[];
            [N, dim, T_steps] = size(p);
            
            % Define K-vectors and Angles
            % Probe fundamental modes based on geometry
            if S.bc == 2 % PBC
                L_scale = 2*S.br;
            else % SBC
                L_scale = 2*S.br; 
            end
            
            k_fundamental = 2*pi/L_scale;
            k_mags = k_fundamental * [1]; % Low k, Med k, High k
            
            thetas = (0:10:360)';
            thetas_rad = deg2rad(thetas(1:end-1)); % Exclude 360 duplicate
            
            nK = length(k_mags);
            nTheta = length(thetas_rad);
            
            % Initialize outputs
            S_std = zeros(nTheta, nK);  Deff_std = zeros(nTheta, nK);
            S_mask = zeros(nTheta, nK); Deff_mask = zeros(nTheta, nK);
            
            max_lag = min(floor(T_steps/3), maxsteps/10); 
            time_lags = (0:max_lag) * dt;
        
            fprintf('Calculating Structure and Dynamics...\n');
        
            % --- 1. UNWRAP & CENTER OF MASS ---
            % Unwrap is needed for PBC Standard calculation
            dp = diff(p, 1, 3);
            
            % Always Apply PBC logic to dp for unwrapping purposes (even for SBC COM check)
            if S.bc == 2 
                dp(dp > S.L/2) = dp(dp > S.L/2) - S.L;
                dp(dp < -S.L/2) = dp(dp < -S.L/2) + S.L;
            end
            
            pu = cat(3, p(:,:,1), zeros(size(dp)));
            pu(:,:,2:end) = dp;
            pu = cumsum(pu, 3);
            
            % COM Subtraction (Wrapped and Unwrapped)
            COM = mean(p, 1);       p_cent = p - COM;
            COMu = mean(pu, 1);     pu_cent = pu - COMu;
        
            % --- 2. DEFINE MASK ---
            % Circular mask (Inscribed for PBC, Full for SBC)
            if S.bc == 1, R_mask = S.br; else, R_mask = S.L/2; end
            
            dist_sq = sum(p_cent.^2, 2); 
            mask = dist_sq < R_mask^2; % Logical mask (N x 1 x T)
            N_eff = mean(sum(mask, 1));
        
            % --- 3. LOOP K-VECTORS ---
            for th_i = 1:nTheta
                theta = thetas_rad(th_i);
                cos_th = cos(theta); sin_th = sin(theta);
                
                for k_i = 1:nK
                    k_val = k_mags(k_i);
                    qx = k_val * cos_th; qy = k_val * sin_th;
                    
                    % --- CALC DENSITY MODES ---
                    
                    % A. STANDARD CALCULATION
                    if S.bc == 2 
                        % PBC: Must use Unwrapped trajectories
                        % squeeze() ensures 2D matrix (N x T)
                        phase = -(qx * squeeze(pu_cent(:,1,:)) + qy * squeeze(pu_cent(:,2,:)));
                        rho_std_t = sum(exp(1i * phase), 1);
                    else
                        % SBC: Must use Wrapped trajectories (Unwrapping is invalid)
                        phase = -(qx * squeeze(p_cent(:,1,:)) + qy * squeeze(p_cent(:,2,:)));
                        rho_std_t = sum(exp(1i * phase), 1);
                    end
                    
                    % B. MASKED CALCULATION (Circular Crop)
                    % Always use Wrapped + Mask
                    phase_m = -(qx * squeeze(p_cent(:,1,:)) + qy * squeeze(p_cent(:,2,:)));
                    % Multiply by mask (logical 0/1)
                    rho_mask_t = sum(exp(1i * phase_m) .* squeeze(mask), 1);
        
                    % --- STORE S(k) ---
                    S_std(th_i, k_i) = mean(abs(rho_std_t).^2) / S.N;
                    S_mask(th_i, k_i) = mean(abs(rho_mask_t).^2) / N_eff;
                    
                    % --- FIT Deff ---
                    Deff_std(th_i, k_i) = get_deff(rho_std_t, max_lag, dt, k_val);
                    Deff_mask(th_i, k_i) = get_deff(rho_mask_t, max_lag, dt, k_val);
                end
            end
        end
    
        %% PLOTTING
        if plotting
            figure('Color','k'); tiledlayout(2,2);
            
            if S.bc == 2, tit='PBC'; else, tit='SBC'; end
            X = rad2deg(thetas_rad);
            
            nexttile; imagesc(X, k_mags, S_std'); axis xy; colorbar;
            title([tit ' Standard S(k)']); xlabel('Theta'); ylabel('k');
            
            nexttile; imagesc(X, k_mags, S_mask'); axis xy; colorbar;
            title([tit ' Masked S(k)']); xlabel('Theta'); ylabel('k');
        
            nexttile; imagesc(X, k_mags, Deff_std'); axis xy; colorbar;
            title([tit ' Standard Deff']); xlabel('Theta'); ylabel('k');
        
            nexttile; imagesc(X, k_mags, Deff_mask'); axis xy; colorbar;
            title([tit ' Masked Deff']); xlabel('Theta'); ylabel('k');
        end
    end
end

% --- HELPER FUNCTION FOR FIT ---
function D = get_deff(rho_t, max_lag, dt, k_val)
    if any(isnan(rho_t)), D=NaN; return; end
    
    rho_c = rho_t - mean(rho_t);
    acf = xcorr(rho_c, max_lag, 'biased');
    acf = acf(max_lag+1:end);
    
    F_norm = abs(acf) / abs(acf(1));
    time = (0:max_lag) * dt;
    
    % Fit window: F between 0.9 and 0.2
    idx = F_norm < 0.9 & F_norm > 0.2;
    
    if sum(idx) < 5
        D = NaN; 
    else
        y = log(F_norm(idx));
        x = time(idx);
        p = polyfit(x, y, 1);
        D = -p(1) / k_val^2;
    end
end
	