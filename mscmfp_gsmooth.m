function out=gsmooth(in,sigma)
% function out=gsmooth(in,sigma)
if nargin<2
    sigma=5;
end;
L = round(10*sigma);
xx=-L:L;
gg = exp(-xx.^2/(2*sigma^2))';

%new version
in1 = 0*in+1;
%commented this and the next 3 lines out to do a 1-D gsmooth
% for k=1:2
    [M,N] = size(in);
    UpS = max(2*M,5*L);
    if M>1
        in (UpS,N) = 0;
        in1(UpS,N) = 0;
        ggk = gg; ggk(UpS)=0; fggk = fft(ggk);
        in  = ifft( fft(in ).*(fggk*ones(1,N)) );
        in1 = ifft( fft(in1).*(fggk*ones(1,N)) );
        in  = in ((L+1):(M+L),:);
        in1 = in1((L+1):(M+L),:);
    end;
%     in = in';
%     in1 = in1';
% end;
out = in./in1;

% [M,N] = size(in);
% if 4*M*N>1e7    out = 0;    return;    end;
% if 4*M*N>1e4
%     in1 = in*0+1;
%     gmask(M*2,N*2) = 0;
%       in1(M*2,N*2) = 0;
%        in(M*2,N*2) = 0;
%     out =        ifft2(fft2(gmask).*fft2(in));
%     out = out ./ ifft2(fft2(gmask).*fft2(in1));
%     out = out((L+1):(L+M),(L+1):(L+N));
% end;

% %old version
% gmask = gg*gg';
% out = conv2(in,gmask)./conv2(in*0+1,gmask);
% out = out(L+1:end-L,L+1:end-L);
