clear all
close all
rng('shuffle')

S.br=1e-6;
S.rp=1e-8;
S.stdx=1e-9;
S.N=1000;
eps=1e-12;

init_r = S.br * rand(S.N, 1).^(1/3);
init_az = 2*pi * rand(S.N, 1);
init_el = asin(2*rand(S.N, 1) - 1);

[x, y, z] = sph2cart(init_az, init_el, init_r);
p0 = [x, y, z];
DISP=build_noise_library(S.stdx,1e7);
qd=1;
debugrhoedges=linspace(0,S.br,1001)';
debugrhoedgesdisp=linspace(0,10*S.stdx,1001)';
debugazedges=linspace(-pi,pi,360)';
debugeledges=linspace(-pi/2,pi/2,180)';

DEBUG{1,1}=zeros(numel(debugazedges)-1,1);
qs=0;
history=[0,0];
while true
    qs=qs+1;
    draw=randi(1e7,[S.N 1]);
    p1 = p0+DISP(draw, :);
    [az,~,~]=cart2sph(p1(:,1),p1(:,2),p1(:,3));
    [temp,debugazedges]=histcounts(az,debugazedges);                   
    DEBUG{1,1}=DEBUG{1,1}+temp';
    rho=vecnorm(p1,2,2);
    idxswap=rho>S.br;
    p1(idxswap,:)=p1(idxswap,:)-(2*S.br).*(p1(idxswap,:)./rho(idxswap,:));
    p0=p1;
    if mod(qs,1e3)==0
        history(end+1,:)=[qs,(max(DEBUG{1,1})-min(DEBUG{1,1}))/mean(DEBUG{1,1})];
        disp([log10(qs),(max(DEBUG{1,1})-min(DEBUG{1,1}))/mean(DEBUG{1,1})])
    end
end
