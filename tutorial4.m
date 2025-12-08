%% CLEANUP

clear
close all
clc
rng('shuffle')

%% flags
io_enabled=true;
plotting=false;


%% parameters

N=100;
rp=1;
sigma=rp*2;
rc=3*sigma;
epsilon=2;
phi=0.85;
ba=(N*pi*rp^2)/phi;
br=sqrt(ba)/2;
bl=2*br;
stdx=rp/100;
D0=1;
T=0.1;
timestep=stdx^2/(2*D0);
steps=1e5;
storageinterval=1;
movieinterval=20;
eps=1e-12;

%% initial positions

% generate force interpolant and displacement library
DISP=build_noise_library(stdx,1e6);
H=pot_force(1,rc,30000,sigma,epsilon);
H_interpolant = griddedInterpolant(H(:,1), H(:,2), 'linear', 'nearest');
DISP=single(DISP(:,1:2));

% initial positions
u1=[1,0];
u2=[0.5,sqrt(3)/2];
v1=u1*2*(rp+eps);
v2=u2*2*(rp+eps);
n=2*ceil(bl/(2*rp));
p=zeros((2*n+1)^2,2);
q=1;
for ix=-n:n
    for iy=-n:n
        p(q,:)=ix*v1+iy*v2;
        q=q+1;
    end
end
p=p+rp/2+2*eps;
idxoob=sum((p>(bl-rp/2.5-eps) | p<(0+rp/2.5+eps)),2)>0;
p(idxoob,:)=[];
p=p(randperm(size(p,1)),:);
p=p(1:N,:);
p_unwrapped=p;
qd=1;
POS=single(zeros(N,2,steps/storageinterval));
if plotting
    figure
    xlim([0 bl])
    ylim([0 bl])
    axis equal
    % --- Circle Parameters ---
    N_points = 10; % Resolution
    
    % 1. Generate angle vector
    theta = linspace(0, 2*pi, N_points);
    
    % 2. Calculate coordinates
    xc = rp * cos(theta);
    yc = rp * sin(theta);
end
for t = 1:steps
    % collect brownian motions
    bdisp=DISP(qd:qd+N-1,:);
    qd=qd+N;
    if qd+N>1e6
        DISP=DISP(randperm(1e6),:);
        qd=1;
    end
    % squared distances (MIC corrected)
    dx = p(:,1)-p(:,1)'; dx = dx-bl*round(dx/bl);
    dy = p(:,2)-p(:,2)'; dy = dy-bl*round(dy/bl);
    r2 = dx.^2+dy.^2; r2(1:N+1:end)=Inf;
    % apply cutoff
    mask = r2 < rc^2;
    r2_c = r2(mask); 
    r=sqrt(r2_c);
    invr=1./r;
    % calculate LJ force and displacement
    f_scalar=H_interpolant(r);
    f_div_r = f_scalar .* invr; 
    fx = zeros(N,N); fy = zeros(N,N);
    fx(mask) = f_div_r .* dx(mask);
    fy(mask) = f_div_r .* dy(mask);
    potdisp=(D0/T)*[sum(fx,2), sum(fy,2)]*timestep;
    overshoot=vecnorm(potdisp,2,2)>2*stdx;
    potdisp(overshoot,:)=2*stdx.*(potdisp(overshoot,:)./vecnorm(potdisp(overshoot,:),2,2));
    % applying displacement
    bdisp=potdisp + bdisp;
    p = p + bdisp;
    p_unwrapped = p_unwrapped + bdisp;
    p = mod(p, bl);
    % storing positions
    if mod(t,storageinterval)==0
        POS(:,:,t/storageinterval)=p_unwrapped;
    end
    % plotting
    if plotting & mod(t,movieinterval)==0
        for imovie=1:N 
            tempx=xc+p(imovie,1);
            tempy=yc+p(imovie,2);
            if imovie==1
                plot(tempx, tempy, 'w-', 'LineWidth', 2);
                hold on
            else
                plot(tempx, tempy, 'w-', 'LineWidth', 2);
            end
        end
        hold off      
        xlim([0 bl])
        ylim([0 bl])
        axis equal
        drawnow
    end
    if mod(t,1000)==0
        disp(100*t/steps)
    end
end
%%

DISP=diff(POS,1,3);
DISP=reshape(permute(DISP, [1, 3, 2]), [], 2);
rhovsaz=[];
[az,rho]=cart2pol(DISP(:,1),DISP(:,2));
azedges=linspace(-pi,pi,36)';
azedges=[azedges(1:end-1,:),azedges(2:end,:)];
for ia=1:size(azedges,1)
    idx=find(az>=azedges(ia,1) & az<azedges(ia,2));
    rhovsaz(ia)=mean(rho(idx,1));
end
rhovsaz=rhovsaz./mean(rho);
figure
hold
plot(azedges(:,1)/pi*180,rhovsaz')

%% timelapse10
POS10=POS(:,:,(1:10:size(POS,3)));
DISP10=diff(POS10,1,3);
DISP10=reshape(permute(DISP10, [1, 3, 2]), [], 2);
rhovsaz10=[];
[az10,rho10]=cart2pol(DISP10(:,1),DISP10(:,2));
for ia=1:size(azedges,1)
    idx=find(az10>=azedges(ia,1) & az10<azedges(ia,2));
    rhovsaz10(ia)=mean(rho10(idx,1));
end
rhovsaz10=rhovsaz10./mean(rho10);

plot(azedges(:,1)/pi*180,rhovsaz10')
%% timelapse100
POS100=POS(:,:,(1:100:size(POS,3)));
DISP100=diff(POS100,1,3);
DISP100=reshape(permute(DISP100, [1, 3, 2]), [], 2);
rhovsaz100=[];
[az100,rho100]=cart2pol(DISP100(:,1),DISP100(:,2));
for ia=1:size(azedges,1)
    idx=find(az100>=azedges(ia,1) & az100<azedges(ia,2));
    rhovsaz100(ia)=mean(rho100(idx,1));
end
rhovsaz100=rhovsaz100./mean(rho100);

plot(azedges(:,1)/pi*180,rhovsaz100')

%% timelapse1000
POS1000=POS(:,:,(1:1000:size(POS,3)));
DISP1000=diff(POS1000,1,3);
DISP1000=reshape(permute(DISP1000, [1, 3, 2]), [], 2);
rhovsaz1000=[];
[az1000,rho1000]=cart2pol(DISP1000(:,1),DISP1000(:,2));
for ia=1:size(azedges,1)
    idx=find(az1000>=azedges(ia,1) & az1000<azedges(ia,2));
    rhovsaz1000(ia)=mean(rho1000(idx,1));
end
rhovsaz1000=rhovsaz1000./mean(rho1000);

plot(azedges(:,1)/pi*180,rhovsaz1000')