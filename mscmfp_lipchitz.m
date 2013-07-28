%mscmfp_lipchitz
maxt = 0:.01:6;
maxt2 = maxt.^2;
maxm = 0*maxt;
minm = 0*maxt+1;

% hold on;
for k=1:lzs*lrs
% for k=1:lzs*lrs
    Gk = G(:,k);
    [ptk,nx] = sort(normByCols(DE*(grid_pts-grid_pts(:,k)*ones(1,lrs*lzs))));
    nxm = min(find(ptk>6));
    ptk = ptk(2:nxm);
    nx = nx(2:nxm);
    GGk = abs(Gk'*G(:,nx));
    dGGk = 1-GGk;
    nxx = ceil(100*ptk); 
    nxx(nxx>601) = 601;
    nx2 = find(dGGk>maxm(nxx));
    maxm(nxx(nx2)) = dGGk(nx2);
    nx2 = find(dGGk<minm(nxx));
    minm(nxx(nx2)) = dGGk(nx2);
    pause(.01);
    if mod(k,100)==0
        nx3 = find(diff(maxm)<0);
        while ~isempty(nx3)
            maxm(nx3+1) = maxm(nx3);
            nx3 = find(diff(maxm)<0);
        end;
        nx3 = find(diff(minm)<0);
        while ~isempty(nx3)
            minm(nx3) = minm(nx3+1);
            nx3 = find(diff(minm)<0);
        end;
        disp(num2str([k,toc/60])); pause(.01); 
    end;
end;
