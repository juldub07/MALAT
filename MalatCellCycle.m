function [MLnuc, Mmit, SLnuc, Smit]=MalatCellCycle(Dpp)

% THis function aims at identifying mitotic vs interphasic nuclei based on
% DAPI signal intensity and shape from max projection
%
% Input argument:
%
% Dpp  maximum projection of the DAPI channel
%
% Output arguments:
%
% MLnuc binary mask of interphasic nuclei
% Mmit  binary mask of mitotic nuclei
% SLnuc properties of interphasic nuclei
% Smit  properties of mitotic nuclei


% tic
% IS=bfopen(dvfile);
% nZ=size(IS{1,1},1)/3;
% 
% Dp=cat(3,IS{1,1}{nZ*2+1:end,1});
% Dpp=max(Dp,[],3);
[a,b]=size(Dpp);
Dppt=imtophat(Dpp,strel('disk',100));
Dpptg=imgaussfilt(Dppt,5);
Dr=imresize(Dpptg,0.25);

% Intensity-based identification of mitotic nuclei
DWmit=Dr>9000;

% Segmentation of all nuclei
[~,Lnuc]=WatershedDAPI(Dr,10,4,100,1);

% calculate region properties of both populations
Smit=regionprops(DWmit,Dr,'Area','MeanIntensity','PixelIdxList','PixelValues','Eccentricity','Solidity');
SLnuc=regionprops(Lnuc,Dr,'Area','PixelIdxList','MeanIntensity','PixelValues');



M=[SLnuc.MeanIntensity];
An=[SLnuc.Area];
tf=M<1000 | An<100;
SLnuc(tf)=[];
for i=1:numel(SLnuc)
    SLnuc(i).Var=std(double(SLnuc(i).PixelValues))/SLnuc(i).MeanIntensity;
end

STD=stdfilt(Dr);
RG=rangefilt(Dr);
ENT=entropyfilt(Dr);

for j=1:numel(Smit)
    Smit(j).Std=mean(STD(Smit(j).PixelIdxList));
    Smit(j).RG=mean(RG(Smit(j).PixelIdxList));
    Smit(j).Ent=mean(ENT(Smit(j).PixelIdxList));
end

% filter mitotic nuclei based on texture
A=[Smit.Area];
E=[Smit.Ent];
tf=A<100 | E<4;
Smit(tf)=[];

nmit=numel(Smit);
nLnuc=numel(SLnuc);

% get interphasic nuclei only
k=false(nLnuc,1);
for i=1:nmit
    A=Smit(i).PixelIdxList;
    for j=1:nLnuc
        if ~isempty(intersect(A,SLnuc(j).PixelIdxList))
            k(j)=true;
        end
    end
end

SLnuc(k)=[];

% create interphasic mask
MLnuc=false(size(Lnuc));
Pnuc=vertcat(SLnuc.PixelIdxList);
MLnuc(Pnuc)=true;
MLnuc=imresize(MLnuc,[a b]);
MLnuc=bwmorph(MLnuc,'thicken',5);
%SMLnuc=regionprops(MLnuc,'Centroid');


% create mitotic mask
Mmit=false(size(Lnuc));
Pmit=vertcat(Smit.PixelIdxList);
Mmit(Pmit)=true;
Mmit=imresize(Mmit,[a,b]);
Mmit=bwmorph(Mmit,'thicken',10);



%SMmit=regionprops(Mmit,'Centroid');
% 
% %Dps=sum(Dp,3);
% Dps=mean(Dp,3);
% 
% Data=regionprops(Lnuc,Dps,'Area','Centroid','PixelValues','MeanIntensity','Eccentricity','MajorAxisLength','MinorAxisLength');
% for i=1:numel(Data)
%     Data(i).IntInt=sum(Data(i).PixelValues);
%     Data(i).ARatio=Data(i).MajorAxisLength/Data(i).MinorAxisLength;
% end
% 
% Data=rmfield(Data,{'PixelValues'});
% tf=[Data.MeanIntensity]<2000;
% Data(tf)=[];
% tfg1=[Data.IntInt]<=8e7;
% DataG1=Data(tfg1);
% DataSG2M=Data(~tfg1);

% BLnuc=bwboundaries(MLnuc);
% Bmit=bwboundaries(Mmit);
% 
% figure, imshow(Dpp,[]);hold on;
% 
% for i=1:length(BLnuc)
%     bb=BLnuc{i};
%     plot(gca,bb(:,2),bb(:,1),'-g');
% end
% 
% CC=[SMLnuc.Centroid];
% 
% for i=1:2:numel(CC)
%     text(CC(i),CC(i+1),int2str((i+1)/2),'Color','y');
% end
% 
% 
% for i=1:length(Bmit)
%     bb=Bmit{i};
%     plot(gca,bb(:,2),bb(:,1),'-r');
% end
% 
% CC=[SMmit.Centroid];
% 
% for i=1:2:numel(CC)
%     text(CC(i),CC(i+1),int2str((i+1)/2),'Color','m');
% end

%toc


