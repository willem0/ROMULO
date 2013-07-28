function [maxx,nx] = minimax_assignment(xx,yy)
% function [maxx,nx] = minimax_assignment(xx,yy)
% 
% Finds the one-to-one assignment between 
% the column vectors in xx and the column vectors in yy
% that minimizes the maximum distance between the resulting pairs:
% maxx := max(sqrt(  sum((yy(:,nx) - xx)).^2 ))
% 
% (such approaches generally work well under uniform or gaussian assumptions
% on errors between the two sets of vectors)
if nargin==1
    D = xx;
else
    if norm(size(xx)-size(yy))~=0
        error('xx and yy must be the same size');
    end;
    [M,N] = size(xx);
    
    %first column of DD corresponds to yy(:,1)
    D = sqrt((ones(N,1)*(sum(xx.^2,1)))' + ones(N,1)*sum(yy.^2,1) - 2*xx'*yy);
end;

[nx,maxx] = munkres(D.^2);
if maxx==0
    return;
end;
Dk = D.^2/maxx;
l2 = sqrt(maxx/length(nx));
for k=1:10
    [nx,y] = munkres(Dk);
    Dk = (Dk/y).^10;
end;
maxx = max(diag(D(:,nx)));

maxx / l2;
% [y,dx] = max(min(D));
% maxx = y;
% nx = dx;
