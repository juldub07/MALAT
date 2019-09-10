function  T=malat2(dvfile)

% This script gathers metrics related to number and size of MALAT and NEAT
% lncRNA aggregates in the nucleus and the cytoplasm of cultured cells.
% This script has been optimized for only one imaging setting using a
% DeltaVision deconvolution microscope (GE).
%
% Input argument:
%
% dvfile: name of the deconvolved dv file (*D3D.dv) to be analyzed. This
% file must be a z-stack with 3 channels: DAPI, MALAT FISH, NEAT FISH
%
% Output argument:
%
% T: a table containing single cell (one cell per row) metrics:
% NucArea: size of the nucleus (in pixels)
% NucIntensity: Mean DAPI intensity
% NucVar: variance of nuclear signal
% NucEcc: Eccentricity of nuclear mask
% datanucMALAT: number of malat aggreagates in the nucleus
% datacytoMALAT: number of malat aggregates in the cytoplasm
% datanucNEAT: number of NEAT aggregates in the nucleus
% datacytoNEAT: number of NEAT aggregates in the cytoplasm
% size_cellcompartment_lncRNA_stats_: mean or median size (in pixel) of a
% NEAT or MALAT aggregates in the nucleus or cytoplasm

tic

%% load file and arrange channels in distinct stacks

Istack=bfopen(dvfile);
ni=size(Istack{1,1},1);
nic=ni/3;
Rstack=cat(3,Istack{1,1}{1:nic,1});
Gstack=cat(3,Istack{1,1}{nic+1:2*nic,1});
Dstack=cat(3,Istack{1,1}{2*nic+1:3*nic,1});

Dproj=max(Dstack,[],3); %maximum projection of DAPI stack

%% Identify mitotic vs interphasic nuclei

[MLnuc, Mmit]=MalatCellCycle(Dproj);

%% Create nuclear mask

Dg=imgaussfilt(Dproj,4);
DW=imbinarize(Dg);
DWf=imfill(DW,'holes');
DWfc=imclose(DWf,strel('disk',10));
%Dproj=adapthisteq(Dproj);
%Gproj=max(Gstack,[],3);
%T=graythresh(Gproj);
%DW=imbinarize(Dproj);
%DW=imclearborder(DW);

%% Segment and quantify MALAT and NEAT granules in each cellular compartment for a given nuclear population (mitotic-Mmit or interphasic- MLnuc)

[a,b,c]=size(Gstack);
GW=false(a,b,c);
RW=false(a,b,c);

% Dr=imresize(Dproj,0.25);
% [~, Lnuc]=WatershedDAPI(Dr,8,2,50);
%LnucR=Lnuc>0;
%LnucR=Lnuc;
%LnucR=imresize(LnucR,[a b]);

%%%comment/uncomment LnucR=MLnuc or LnucR=Mmit to analyze either
%%%interphasic or mitotic nuclei

%LnucR=MLnuc;
LnucR=Mmit;

%%%

if any(LnucR(:))
    
    LnucRL=bwlabel(LnucR);
    
    statsNUC=regionprops(LnucRL,Dproj,'Area','PixelIdxList','PixelValues','MajorAxisLength','MinorAxisLength','Eccentricity','MeanIntensity');
    %A=[statsNUC.Area];
    %tf=A<1000;
    %statsNUC(tf)=[];
    % LnucR=false(a,b);
    % Pix=vertcat(statsNUC.PixelIdxList);
    % LnucR(Pix)=true;
    LcytoR=bwmorph(LnucR,'thicken',50);
    statsCYTO=regionprops(LcytoR,'PixelIdxList');
    LcytoR=bwlabel(LcytoR);
    
    
    for i=1:numel(statsNUC)
        Pn=statsNUC(i).PixelIdxList;
        for j=1:numel(statsCYTO)
            Pc=statsCYTO(j).PixelIdxList;
            if ~isempty(intersect(Pn,Pc))
                LcytoR(Pc)=i;
            end
        end
    end
    
    
    datanucMALAT=zeros(numel(statsNUC),1);
    datacytoMALAT=zeros(numel(statsNUC),1);
    datanucNEAT=zeros(numel(statsNUC),1);
    datacytoNEAT=zeros(numel(statsNUC),1);
    
    for i=1:c
        Gstack(:,:,i)=imtophat(Gstack(:,:,i),strel('disk',15));
        Rstack(:,:,i)=imtophat(Rstack(:,:,i),strel('disk',15));
    end
    
    TG=graythresh(max(Gstack,[],3));
    TR=graythresh(max(Rstack,[],3));
    %T=graythresh(Gstack(:,:,13));
    %T=0.0549;
    
    
    for j=1:c
        GW(:,:,j)=imbinarize(Gstack(:,:,j),TG);
        RW(:,:,j)=imbinarize(Rstack(:,:,j),TR);
    end
    
    
    for k=1:max(LcytoR(:))
        
        MASKN=LnucRL==k;
        MASKC= (LcytoR==k) & ~DWfc & ~MASKN;
        GWn=GW & MASKN;
        GWc=GW & MASKC;
        RWn=RW & MASKN;
        RWc=RW & MASKC;
        
        
        
        
        
        
        
        
        
        
        
        
        
        % %rd=randi(numel(statsCYTO));
        % for k=1:numel(statsCYTO)
        %     %Box=statsCYTO(k).BoundingBox;
        %     MASKN=false(size(DWfc));
        %     MASKC=false(size(DWfc));
        %     Pix=vertcat(statsNUC(k).PixelIdxList);
        %     MASKN(Pix)=true;
        %     %MASKN=imcrop(MASKN,Box);
        %     %MASKNN=imcrop(LnucR,Box);
        %     %MASKG=imcrop(Gproj,Box);
        %     %figure, imshow(MASKG,[]);
        %     %MASKN=repmat(MASKN,1,1,c);
        %     Pixc=vertcat(statsCYTO(k).PixelIdxList);
        %     MASKC(Pixc)=true;
        %     %MASKC=imcrop(MASKC,Box);
        %     %DWc=imcrop(DWfc,Box);
        %     MASKCC=MASKC & ~DWc & ~MASKN;
        %
        %     %MASKC=repmat(MASKC,1,1,c);
        %     GWcn=false([size(MASKC), c]);
        %     GWcc=false([size(MASKC),c]);
        %     RWcn=false([size(MASKC),c]);
        %     RWcc=false([size(MASKC),c]);
        %     for p=1:c
        %     GWcn(:,:,p)=GW(:,:,p) & MASKN;
        %     GWcc(:,:,p)=GW(:,:,p) & MASKCC;
        %     RWcn(:,:,p)=RW(:,:,p) & MASKN;
        %     RWcc(:,:,p)=RW(:,:,p) & MASKCC;
        %     end
        
        
        
        
        %GWcn=bwareaopen(GWcn,10,6);GWcc=bwareaopen(GWcc,10,6);
        %RWcn=bwareaopen(RWcn,10,6);RWcc=bwareaopen(RWcc,10,6);
        
        CGn=bwconncomp(GWn); CGc=bwconncomp(GWc);
        CRn=bwconncomp(RWn); CRc=bwconncomp(RWc);
        datanucMALAT(k)=CGn.NumObjects;datacytoMALAT(k)=CGc.NumObjects;
        datanucNEAT(k)=CRn.NumObjects;datacytoNEAT(k)=CRc.NumObjects;
        
        sizenucMALATmedian(k)=median(cellfun(@numel,CGn.PixelIdxList));sizenucMALATmean(k)=mean(cellfun(@numel,CGn.PixelIdxList));
        sizecytoMALATmedian(k)=median(cellfun(@numel,CGc.PixelIdxList));sizecytoMALATmean(k)=mean(cellfun(@numel,CGc.PixelIdxList));
        
        sizenucNEATmedian(k)=median(cellfun(@numel,CRn.PixelIdxList));sizenucNEATmean(k)=mean(cellfun(@numel,CRn.PixelIdxList));
        sizecytoNEATmedian(k)=median(cellfun(@numel,CRc.PixelIdxList));sizecytoNEATmean(k)=mean(cellfun(@numel,CRc.PixelIdxList));
 
        %%% uncomment section below for visual quality control of random
        %%% cells
        
        %  if k==rd
        % %     figure, imshow(MASKN);
        % %     figure, imshow(MASKC);
        %      figure, imshow(MASKCC);
        % figure, imshow(max(GWcc,[],3));
        % figure, imshow(max(RWcc,[],3));
        %     [d,e,f]=size(GWcn);
        %     [x, y, z]=meshgrid(1:e,1:d,1:f);
        %     figure, isosurface(x,y,z,GWcn,0.5);
        %     figure, isosurface(x,y,z,RWcn,0.5);pause
        %  end
        
        %%%
    end
    
    
    NucArea=[statsNUC.Area]';
    NucIntensity=[statsNUC.MeanIntensity]';
    NucCell=struct2cell(statsNUC);
    NucVar=cellfun(@double,NucCell(6,:),'UniformOutput',0)';
    NucVar=cellfun(@var,NucVar)./NucIntensity;
    NucEcc=[statsNUC.Eccentricity]';
    sizenucMALATmedian=sizenucMALATmedian';
    sizenucMALATmean=sizenucMALATmean';
    sizecytoMALATmedian=sizecytoMALATmedian';
    sizecytoMALATmean=sizecytoMALATmean';
    sizenucNEATmedian=sizenucNEATmedian';
    sizenucNEATmean=sizenucNEATmean';
    sizecytoNEATmedian=sizecytoNEATmedian';
    sizecytoNEATmean=sizecytoNEATmean';
    T=table(NucArea,NucIntensity,NucVar,NucEcc,datanucMALAT,datacytoMALAT,...
        datanucNEAT,datacytoNEAT,sizenucMALATmedian,sizenucMALATmean,...
        sizecytoMALATmedian,sizecytoMALATmean,sizenucNEATmedian,...
        sizenucNEATmean,sizecytoNEATmedian,sizecytoNEATmean);
    
else
    T=table;
end
%Ttmp=[datanucMALAT,datacytoMALAT,datanucNEAT,datacytoNEAT];
%tf0=sum(Ttmp,2)==0;
%T(tf0,:)=[];
%statsNUC=rmfield(statsNUC,{'BoundingBox','PixelIdxList'});

toc

