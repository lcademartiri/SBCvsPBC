clear all
close all
rng('shuffle')

load('conditionsfordummy.mat');
V.c(:,[5,14])=[];
V=table2array(V.c);
V=V(:,[1,2,3,5,10,12,16]);
V(:,end+1)=(1.38e-23*298)./(6*pi.*V(:,7).*V(:,4));
V(:,end+1)=sqrt(2*0.7.*V(:,8).*V(:,5));
V(V(:,2)~=1,:)=[];
V(:,[1,2,4,5,7,8])=[];
V=unique(V,'rows');

debugazedges=linspace(-pi,pi,360)';
for ic=1:size(V,1)
    S.br=V(ic,2);
    S.stdx=V(ic,3);
    S.N=V(ic,1);
    
    init_r = S.br * rand(S.N, 1).^(1/3);
    init_az = 2*pi * rand(S.N, 1);
    init_el = asin(2*rand(S.N, 1) - 1);
    
    [x, y, z] = sph2cart(init_az, init_el, init_r);
    p0 = [x, y, z];
    DISP=build_noise_library(S.stdx,1e7);
    qd=1;
    
    DEBUG{ic,1}=zeros(numel(debugazedges)-1,1);
    qs=0;
    history{ic,1}=[0,0];
    while qs<=1e7
        qs=qs+1;
        draw=randi(1e7,[S.N 1]);
        p1 = p0+DISP(draw, :);
        [az,~,~]=cart2sph(p1(:,1),p1(:,2),p1(:,3));
        [temp,debugazedges]=histcounts(az,debugazedges);                   
        DEBUG{ic,1}=DEBUG{ic,1}+temp';
        rho=vecnorm(p1,2,2);
        idxswap=rho>S.br;
        p1(idxswap,:)=p1(idxswap,:)-(2*S.br).*(p1(idxswap,:)./rho(idxswap,:));
        p0=p1;
        if mod(qs,1e3)==0
            convergence=(max(DEBUG{ic,1})-min(DEBUG{ic,1}))/mean(DEBUG{ic,1});
            history{ic,1}(end+1,:)=[qs,convergence];
            disp([ic,log10(qs),convergence])
        end
    end
end