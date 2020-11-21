%% Further segment nuclei from tag-segmentation
function [NUC, Bn, Ln]= Nuclear2(D,varargin)



% Segment DAPI nuclei using watershed strategy
Dl=bwdist(~D);
if ~isempty(varargin)
    s=varargin{1};
h=fspecial('disk',s); %start with 16
Dl=imfilter(Dl,h);
end
Dl=-Dl;
Dl(~D)=-Inf; % -Inf
%NUC=watershed(Dl,8);
NUC=watershed(Dl);
NUC(~D)=0;
NUC=logical(NUC);





[Bn, Ln, Nn]=bwboundaries(NUC,4,'noholes');












