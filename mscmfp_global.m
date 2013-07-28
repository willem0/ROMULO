%mscmfp_global;
%inputs: PY, PG, PGr

%should PY be normalized all the time?
%presumed coherent
t0 = tic;

PGr = reshape(permute(PG,[1 3 2]),M*KF,lrs*lzs);
PGr_save = PGr;
PY_save = PY;

PGr = PGr ./ (ones(M*KF,1) * normByCols(PGr));
PY = PY(:)/norm(PY(:));

PGr0 = PGr;
PY0 = PY;
gam = 0.5;
p = 2/(1+abs(gam));
th = 2*((1-abs(gam))/2)^(p/2);

YPPG = abs(PY0'*PGr0);
[y,nx] = sort((YPPG).^p,'descend');
nxe = max(find(y(1)+y>th));

% for l=1:lzs*lrs     %global search over all locations
reg = ones(lzs*lrs,1);
for ll=1:nxe     %global search over all locations
    l = nx(ll);
    PGl = PGr0(:,l);
    nxle = max(find(y(ll)+y>th));
    PGrl = PGr0(:, nx(1:nxle) );
    PGrl = PGrl - PGl*(PGl'*PGrl);
    PYl   = PY0  - PGl*(PGl'*PY );
    
    ht = abs(PYl'*PGrl).^2 ./ (sum(abs(PGrl).^2) + 2^16*eps);
    
    [rek1,cmnx] = max(ht);
    nxg(l) = nx(cmnx);
    A = [PGl PGrl(:,cmnx)];
    reg(l) = 1-norm(A*pinv(A) * PYl)/norm(PYl);
    if mod(l,1000)==1, toc(t0), pause(.01); end;
end;

[y,nxgg] = min(reg);
reg = y;
rhat = [grid_pts(:,nxgg), grid_pts(:,nxg(nxgg))];


return;
%30 minutes per script without speedup
