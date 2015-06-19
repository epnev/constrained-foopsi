clear;
T = 10000;
gr1 = [0.85,0.9];
g1 = poly(gr1);
gr2 = [0.7,0.98];
g2 = poly(gr2);
ps = 0.01;
sps = rand(1,T)<ps;

r = gamrnd(repmat([2,2],T,1),1,T,2);
r = r ./ repmat(sum(r,2),1,2);
ct1 = filter(1,g1,sps.*r(:,1)')';
ct2 = filter(1,g2,sps.*r(:,2)')';
c = ct1 + ct2;
sn = .125;
spt = [sps.*r(:,1)';sps.*r(:,2)']';
y = c + sn*randn(size(c));

figure;plot(1:T,y);

%%

grin1 = [0.3,0.95];
grin2 = [0.5,0.85];
gin1 = poly(grin1);
gin2 = poly(grin2);

gin{1} = -gin1(2:3)';
gin{2} = -gin2(2:3)';

[cm,b,c1,sp] = cvx_mixed_foopsi(y(:),0,[0;0],sn,0,gin,ones(T,1),true(T,1));
%% iterate a few times
ITER = 15;
gr_ = [grin1(:),grin2(:)];
c_ = cm;
for iter = 1:ITER
    rp = randperm(2);
    for i = 1:2
        [tau_new, c_new, accept] = propose_new_time_constant(y(:)-c_(:,[1:rp(i)-1,rp(i)+1:2]),sp(:,i),-1./log(gr_(:,i)),0,200,[0.2,0.5],sn,0);
        c_(:,rp(i)) = c_new;
        gr_(:,rp(i)) = exp(-1./tau_new(:));
        gn{rp(i)} = [sum(gr_(:,rp(i)));-prod(gr_(:,rp(i)))];
        disp(accept)
    end
    [c_,b,c1,sp] = cvx_mixed_foopsi(y(:),0,[0;0],sn,0,gn,ones(T,1),true(T,1));
    disp(corr(spt,sp))
end
    

