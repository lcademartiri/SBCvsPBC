%% CLEANUP

clear all
close all
clc
addpath('..\ARBD_toolbox');

%% FLAGS

moviemaking=false;
allpos=true;

%% SIMULATION PARAMETERS
S.rp=1;
dp=2*S.rp;
epsilon=1;
sigma=2*S.rp;
S.rc=3*sigma;
maxsteps=1e6;
reps=1;
S.stdx=S.rp/50;
nodisps=1e6;
S.timestep=1;
S.esdiff=S.stdx^2/(2*S.timestep);
S.kbT=1;

%% VARIABLES

Ns=round(logspace(log10(1e3),log10(1e3),1)',0);
phis=logspace(log10(1e-1),log10(0.5),5)';
pots=(1:2)';
bcs=(1:2)';

%% UTILITIES
uhex=[1,0;0.5,sqrt(3)/2];
[~, cmdout] = system('wmic cpu get L2CacheSize, L3CacheSize /value');
tokens = regexp(cmdout, '\d+', 'match');
S.cacheSizeMB = (max(str2double(tokens))/1024)/feature('numCores');

%% CONDITIONS
q=1;
for i0=1:numel(Ns)
	N=Ns(i0);
	for i1=1:numel(phis)
		phi=phis(i1);
		for i2=1:numel(pots)
			pot=pots(i2);
			for i3=1:numel(bcs)
				c(q,:)=[q,N,phi,pot,bcs(i3)];
				q=q+1;
			end
		end
	end
end

%% LOOP

for ic=18
	S.N=c(ic,2);
	S.phi=c(ic,3);
	S.pot=c(ic,4);
	S.bc=c(ic,5);
	S.A=(S.N*pi*S.rp^2)/S.phi;

    % DEFINE DISPLACEMENT LIBRARY
	BD=single(build_noise_library(S.stdx,nodisps));
	BD(:,3)=[];
	qd=1;

    % INITIALIZE POS MATRIX
    POS=single(zeros(S.N,2,maxsteps));

    % DEFINE FORCE AND CLAMP
    if S.pot~=0
	    H=pot_force(S.pot,S.rc,30000,sigma,epsilon);
        H_interpolant = griddedInterpolant(H(:,1), H(:,2), 'linear', 'nearest');
        S.pot_clamp=mcdClamp2D(nodisps, BD, S.esdiff, 1, H, 1);
    end	

	% DEFINE STARTING CONFIGURATION
    if S.bc==1
        S.R=sqrt(S.A/pi);
        S.L=2*S.R;
        S.br=S.R;
        maxcoeff=ceil(S.L/dp);
		p=single([]);
		temp=(-maxcoeff:maxcoeff)';
		tempn=numel(temp);
		for i0=-maxcoeff:maxcoeff
			p=[p;temp,i0.*ones(tempn,1)];
		end
		p=p(:,1).*uhex(1,:).*dp+p(:,2).*uhex(2,:).*dp;
        mask=vecnorm(p,2,2)>S.R-S.rp;
        p(mask,:)=[];
        if size(p,1)<S.N
            disp('too dense')
            pause
        end
		p=p(randsample(size(p,1),S.N),:);
    elseif S.bc==2
		S.L=sqrt(S.A);
        S.br=0.5*S.L;
        maxcoeff=ceil(sqrt(2)*S.L/dp);
		p=single([]);
		temp=(-maxcoeff:maxcoeff)';
		tempn=numel(temp);
		for i0=-maxcoeff:maxcoeff
			p=[p;temp,i0.*ones(tempn,1)];
		end
		p=p(:,1).*uhex(1,:).*dp+p(:,2).*uhex(2,:).*dp;
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
            potdisp = potential_displacements_2D_v13(p, S, H, H_interpolant, 0, S.cacheSizeMB);
            bdisp=bdisp+potdisp;
		end
		% apply displacements
		p=p+bdisp;
        
        % oob
        if S.bc==1
            idxswap=vecnorm(p,2,2)>S.br;
            if sum(idxswap)>0
                p(idxswap,:)=p(idxswap,:)-(2*S.br).*(p(idxswap,:)./vecnorm(p(idxswap,:),2,2));
            end
        elseif S.bc==2
            p=mic_wrap_positions(p, S);
        end

        % counter
        if mod(it,1e3)==0
            disp([ic,log10(it)])
        end
		
        % movie
        if mod(it,10)==0 & moviemaking
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
	
	% Inputs:
	%   p        : (N x 2 x T) matrix. 
	%              N = particles, 2 = (x,y), T = time steps.
	%              If T=1, only S is calculated; F and Deff will be NaN.
	%   k_mags   : Vector of wavenumber magnitudes to probe (e.g., 0.1:0.1:10)
	%   thetas   : Vector of angles (in radians) to probe (e.g., 0:pi/4:pi)
	%   dt       : Time step between frames in p (for scaling Deff and time axis)
	%
	% Outputs:
	%   S_map    : (nTheta x nK) Static Structure Factor S(theta, |k|)
	%   F_map    : (nTheta x nK x nLags) Intermediate Scattering Function
	%   Deff_map : (nTheta x nK) Effective diffusion coefficient
	%   time_lags: (1 x nLags) Time vector corresponding to the 3rd dim of F_map
	%
	% Definitions:
	%   rho(k,t) = sum_j exp(-i * k . r_j(t))
	%   S(k)     = < |rho(k,t)|^2 > / N
	%   F(k,tau) = < rho(k, t+tau) * rho(-k, t) > / (N * S(k))
	%   Deff(k)  = Extracted via fit: F(k,tau) ~ exp(-k^2 * Deff * tau)
    p=POS;
    clear POS
    dt=1;
	kmin=2*pi/S.L;
	kmax=2*pi/S.rp;
	k_mags=kmin;
	thetas=(0:5:360)';
	thetas=(2*pi).*(thetas/360);
	thetas(end,:)=[];
	[N, dim, T_steps] = size(p);
	
	if dim ~= 2
		error('Input p must be Nx2xT');
	end

	nK = length(k_mags);
	nTheta = length(thetas);
	
	% Initialize outputs
	S_map = zeros(nTheta, nK);
	Deff_map = zeros(nTheta, nK);
	
	% Max lag time for ISF (usually don't need full correlation length)
	% Using 50% of trajectory or max 1000 frames to save memory/time
	max_lag = min(floor(T_steps/2), 1000); 
	if T_steps == 1, max_lag = 0; end
	
	time_lags = (0:max_lag) * dt;
	F_map = zeros(nTheta, nK, max_lag + 1);

	% Pre-compute exponentials for all k-vectors to avoid loops
	% We process one (theta, k) pair at a time to manage memory, 
	% or vectorize if N is small. Here we loop over k-space for clarity/memory safety.
	
	fprintf('Calculating Structure and Dynamics...\n');
	
	for th_i = 1:nTheta
		theta = thetas(th_i);
		cos_th = cos(theta);
		sin_th = sin(theta);
		
		for k_i = 1:nK
			k_val = k_mags(k_i);
			
			% Define q vector components
			qx = k_val * cos_th;
			qy = k_val * sin_th;
			
			% --- 1. Calculate Density Mode rho(k, t) for all t ---
			% rho_t will be 1 x T_steps complex vector
			rho_t = zeros(1, T_steps);
            rhou_t = zeros(1, T_steps);
            rhoc_t = zeros(1, T_steps);
            rhocu_t = zeros(1, T_steps);
            % unwrapped coordinates
            % unwrap along each coordinate
            pu = zeros(size(p));
            pc = cell(size(p),1);
			
			% Vectorized over particles for each frame
			% Note: If T is huge, might be better to loop frames.
			% Here assuming T fits in memory.
            COM = squeeze(mean(p,1)); % 2 x T			
			
			% 1. VECTORIZED UNWRAPPING
			% Calculate differences between all consecutive frames at once
			% dp is (N x 2 x T-1)
			dp = diff(p, 1, 3); 

			% Apply Minimum Image Convention to all dp at once
			dp(dp > S.L/2)  = dp(dp > S.L/2) - S.L;
			dp(dp < -S.L/2) = dp(dp < -S.L/2) + S.L;

			% Reconstruct unwrapped paths using cumulative sum
			% We prepend the first frame so dimensions match
			pu = cat(3, p(:,:,1), zeros(size(dp))); % Initialize
			pu(:,:,2:end) = dp;                     % Fill deltas
			pu = cumsum(pu, 3);                     % Integrate: [p1, p1+dp1, p1+dp1+dp2...]

			% 2. CENTER OF MASS (Vectorized)
			% Calculate COM for Unwrapped (2 x T)
			% mean(pu, 1) gives 1 x 2 x T -> squeeze -> 2 x T
			COMu = squeeze(mean(pu, 1))'; 

			% Prepare COM arrays for broadcasting subtraction (1 x 2 x T)
			COM_wrapped_reshaped = reshape(COM, [1, 2, T_steps]);
			COM_unwrapped_reshaped = reshape(COMu, [1, 2, T_steps]);

			% Subtract COMs from all particles at all times simultaneously
			% MATLAB R2016b+ supports implicit expansion (p is Nx2xT, COM is 1x2xT)
			p_centered = p - COM_wrapped_reshaped;
			pu_centered = pu - COM_unwrapped_reshaped;

			% 3. FOURIER COEFFICIENTS (Vectorized)
			% A. Wrapped
			% Calculate phase for all N particles and T steps: (N x T) matrix
			% We squeeze p(:,1,:) to get N x T matrices
			phase_wrapped = -(qx * squeeze(p_centered(:,1,:)) + qy * squeeze(p_centered(:,2,:)));
			rho_t = sum(exp(1i * phase_wrapped), 1); % Sum down columns (particles)

			% B. Unwrapped
			phase_unwrapped = -(qx * squeeze(pu_centered(:,1,:)) + qy * squeeze(pu_centered(:,2,:)));
			rhou_t = sum(exp(1i * phase_unwrapped), 1);

			% C. Spherical (Masking approach)
			if S.bc == 2
				% Calculate distance squared for all particles/times
				% p_centered is N x 2 x T
				dist_sq = sum(p_centered.^2, 2); % Result is N x 1 x T
				
				% Create a logical mask (N x T)
				% 1 if inside sphere, 0 if outside
				mask = squeeze(dist_sq) < S.br^2; 
				
				% Multiply the complex exponential by the mask (zeros out external particles)
				% Then sum.
				rhoc_t = sum(exp(1i * phase_wrapped) .* mask, 1);
			end


			
			% --- 2. Static Structure Factor S(theta, |k|) ---
			% S(k) = <|rho(k)|^2> / N
			S_val = mean(abs(rho_t).^2) / N;
			S_map(th_i, k_i) = S_val;
			
			% If only one frame, skip dynamics
			if T_steps > 1
				% --- 3. Intermediate Scattering Function F(k, tau) ---
				% Standard autocorrelation calculation
				% Vectorized correlation using xcorr is technically possible but 
				% manual loop handles normalization better for "F=1 at t=0".
				
				acf = zeros(1, max_lag + 1);
				count = zeros(1, max_lag + 1);
				
				% Compute autocorrelation
				% F(tau) = < rho(t+tau) * conj(rho(t)) >
				for tau = 0:max_lag
					% Elements we can multiply (t and t+tau)
					val = rhou_t(1+tau : end) .* conj(rhou_t(1 : end-tau));
					acf(tau+1) = mean(val);
				end
				
				% Normalize F usually so F(0) = 1 (or S(k))
				% Strictly F(k,t) / S(k) goes 1->0. 
				% The definition requested implies F(|k|,theta,t), usually normalized by N.
				F_curve = real(acf) / N; 
				
				% Store normalized F (Usually F(k,0) = S(k)).
				F_map(th_i, k_i, :) = F_curve;
				
				% --- 4. Effective Diffusion Deff ---
				% Model: F(k, t) = S(k) * exp(-k^2 * Deff * t)
				% ln(F(k,t)/S(k)) = -k^2 * Deff * t
				
				% We normalize by S(k) for fitting
				F_norm = F_curve / S_val;
				
				% Fit range: Fit the initial decay (e.g., down to 0.1 or first 10-20 points)
				% Find index where F drops below noise or 1/e
				cutoff_idx = find(F_norm < 0.1, 1);
				if isempty(cutoff_idx), cutoff_idx = max_lag+1; end
				
				% Don't fit just 1 point
				fit_len = max(5, floor(cutoff_idx/2));
				fit_len = min(fit_len, max_lag+1);
				
				y_fit = log(max(F_norm(1:fit_len), 1e-6)); % avoid log(0)
				x_fit = time_lags(1:fit_len)';
				
				% Linear regression: y = slope * x
				% slope = -k^2 * Deff
				if k_val > 1e-6
					b = x_fit \ y_fit.'; % Simple least squares
					slope = b;
					Deff_calc = -slope / (k_val^2);
				else
					Deff_calc = 0;
				end
				
				Deff_map(th_i, k_i) = max(0, Deff_calc); % Ensure positive
			else
				F_map(th_i, k_i, :) = NaN;
				Deff_map(th_i, k_i) = NaN;
            end
            disp([ic,th_i,k_i])
		end
    end
end

