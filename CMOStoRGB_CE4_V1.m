%% This program composites RGB images from the hyperspectral data of Chinese Chang'e-4, Yutu-2 
% The CE-4 VNIS (hyperspectral) data can be downloaded from the official website: https://moon.bao.ac.cn/ce5web/searchOrder_dataSearchData.search
% Specific the "outpu_path" and "file_path" before using
% Specific the "COR_mod" before using
% If using the white balance coefficients calculated from the calibration panel (i.e. use the Sun as standard white), the color temperatue is about 5500K
% Version：V1, the initial version
% Date: Feb.2023
% Author: Zhenxing Zhao, China, Beijing, NSSC, Chinese Academy of Sciecnes
% (2020-2025)
% If you have any questions, please contact cooperzhaozx@gmail.com 
% My website (Chinese): cooperzzx.com
restoredefaultpath
clear
close all

%% initialization
output_path='*\CMOS';
file_path='*CE4-data\VNIS'; %the file path that store the CE-4 VNIS 2B data
RGB_modtb=readtable('CMOS2RGBCorMod.xlsx'); %read the coefficients of different color mode
RGB_mod=table2array(RGB_modtb);
Max_normal='N';%Y:normalize all the images using the same constant，N:normalize each the images according to its maxminum
COR_mod='CIE';% "3bands":450,540,565nm, "Stock":Stockman trichromatic vision, "CIE":CIE1931 trichromatic vision
Wb_man=[0.808,0.914,1.498];%default white balance value
CMOS_WL=450:5:945; % the wavelength range of the CE4 CMOS hyperspectral data
if strcmpi(COR_mod,'CIE')
    Cor_ind=[RGB_mod(:,1),RGB_mod(:,2:4)];
    COR_mod_name='CIE';
elseif strcmpi(COR_mod,'Stock')
    Cor_ind=[RGB_mod(:,1),RGB_mod(:,5:7)];
    COR_mod_name='Stock';
elseif strcmpi(COR_mod,'3bands')
    Cor_ind=[RGB_mod(:,1),RGB_mod(:,8:10)];
    COR_mod_name='3bands';
end
%% start
tic
cd(file_path);
dir_VC=dir(fullfile(file_path,'*VC_SCI*.2B'));%the files of CMOS in calibration mode (calibration panel spectral data)
file_VC={dir_VC.name}';
file_VC=cell2mat(file_VC);   % convert cell to matrix.
fileNum_VC = size(file_VC,1); % count the total number of files.
Wb_VC_all=zeros(fileNum_VC,3);
dir_VD=dir(fullfile(file_path,'*VD_SCI*.2B'));%the files of CMOS in detection mode (lunar soil spectral data)
file_VD={dir_VD.name}';
file_VD=cell2mat(file_VD);   % convert cell to matrix.
fileNum_VD = size(file_VD,1); % count the total number of files.
Wb_VD_all=zeros(fileNum_VD,3);
cd(output_path)
figure();
for i=1:fileNum_VC
    basename=file_VC(i,:);
    sgtitle_name=['N',basename(end-7:end-5),'-D',num2str(day_file_n),'-VC-',COR_mod_name];
    filename_VC=[file_path,'\',basename];
    fID_VC = fopen(filename_VC);
    VC_sor0 = fread(fID_VC,'float')';%read the radiance data
    fclose(fID_VC);
    VC_sor0=reshape(VC_sor0,[65536,100]);%
    [Wb_VC,VC_RGB,VC_RGB_Wb]=CMOS2RGB(CMOS_WL,VC_sor0,Cor_ind,'self');%"self": white balance according the data itself
    VC_sor=reshape(VC_sor0,[256,256,100]);
    VC_sor=permute(VC_sor,[2 1 3]);
    Rad_1band=VC_sor(:,:,23);%23th:560nm
    Wb_VC_all(i,:)=Wb_VC;
    M_VC=mean(VC_sor0(:,23));
    STD_VC=std(VC_sor0(:,23));
    if M_VC<0.05 %if there is no sunlight on the calibarion panel, abandon this file. 0.05 threshold is based on experience
        disp([filename_VC,' is an dark VC measurement']);
        Wb_VC_all(i,:)=[];
    end
    if STD_VC>0.05 %if there is a lot of shadow on the calibarion panel, abandon this file. 0.05 threshold is based on experience
        disp([filename_VC,' is an shadowy VC measurement']);
        Wb_VC_all(i,:)=[];
    end
    VC_RGB=VC_RGB/max(max(VC_RGB));
    VC_RGB=reshape(VC_RGB,[256,256,3]);
    VC_RGB=permute(VC_RGB,[2 1 3]);
    VC_RGB_Wb=VC_RGB_Wb/max(max(VC_RGB_Wb));
    VC_RGB_Wb=reshape(VC_RGB_Wb,[256,256,3]);
    VC_RGB_Wb=permute(VC_RGB_Wb,[2 1 3]);
    subplot(1,3,1)
    imagesc(Rad_1band);
    title(['Rad ',num2str(CMOS_WL(23)),' nm']);
    colormap gray
    set(gca, 'XTick',[],'YTick',[]);
    subplot(1,3,2)
    imagesc(VC_RGB);
    title('NO WB');
    set(gca, 'XTick',[],'YTick',[]);
    subplot(1,3,3)
    imagesc(VC_RGB_Wb);
    title(['WB (alpha: ',num2str(Wb_VC(1),'%1.2f'),', ',num2str(Wb_VC(2),'%1.2f'),', ',num2str(Wb_VC(3),'%1.2f'),')']);
    set(gca, 'XTick',[],'YTick',[]);
    sgtitle(sgtitle_name,'Interpreter','none');
    set(gcf,'Color','w');
    set(gcf,'Position',[500,50,1410,400]);
    img_file_name=[sgtitle_name,'.jpg'];
    exportgraphics(gcf,img_file_name,'Resolution',300)
end
%%
if fileNum_VC==0
    disp([file_path,' no VC file, Wb_VD will be manually assigned']);
    Wb_VD=[1.01327568727212,1,0.924242434049424];
elseif isempty(Wb_VC_all) %is there is no calibration data
    disp([file_path,' invalid VC file, Wb_VD will be manually assigned']);
    Wb_VD=[1.01327568727212,1,0.924242434049424];
else
    Wb_VD=mean(Wb_VC_all,1); % white balance coefficients according the data of calibration panel
end
figure();
for i=1:fileNum_VD
    basename=file_VD(i,:);
    sgtitle_name=['N',basename(end-7:end-5),'-D',num2str(day_file_n),'-VD-',COR_mod_name];
    filename_VD=[file_path,'\',basename];
    fID_VD = fopen(filename_VD);
    VD_sor0 = fread(fID_VD,'float')';
    fclose(fID_VD);
    VD_sor0=reshape(VD_sor0,[65536,100]);%
    [Wb_VD_out,VD_RGB,VD_Wb_VC]=CMOS2RGB(CMOS_WL,VD_sor0,Cor_ind,'mannual',Wb_VD);
    VD_sor=reshape(VD_sor0,[256,256,100]);
    VD_sor=permute(VD_sor,[2 1 3]);
    VD_RGB=VD_RGB/max(max(VD_RGB));
    VD_RGB=reshape(VD_RGB,[256,256,3]);
    VD_RGB=permute(VD_RGB,[2 1 3]);
    max_VD_Wb_VC(i)=max(max(VD_Wb_VC));
    VD_Wb_VC=reshape(VD_Wb_VC,[256,256,3]);
    VD_Wb_VC=permute(VD_Wb_VC,[2 1 3]);
    subplot(2,2,1)
    Band_i=VD_sor(:,:,23);
    imagesc(Band_i);
    title(['Rad ',num2str(CMOS_WL(23)),' nm']);
    colormap gray
    set(gca, 'XTick',[],'YTick',[]);
    subplot(2,2,2)
    imagesc(VD_RGB);
    title('NO WB');
    set(gca, 'XTick',[],'YTick',[]);
    subplot(2,2,3)
    imagesc(VD_Wb_VC./max_VD_Wb_VC(i));
    title(['VC WB (alpha: ',num2str(Wb_VD(1),'%1.2f'),', ',num2str(Wb_VD(2),'%1.2f'),', ',num2str(Wb_VD(3),'%1.2f'),')']);
    set(gca, 'XTick',[],'YTick',[]);
    
    %self white balance
    subplot(2,2,4)
    [Wb_VD_auto,VD_RGB,VD_Wb_self]=CMOS2RGB(CMOS_WL,VD_sor0,Cor_ind,'self');
    VD_Wb_self=VD_Wb_self/max(max(VD_Wb_self));
    VD_Wb_self=reshape(VD_Wb_self,[256,256,3]);
    VD_Wb_self=permute(VD_Wb_self,[2 1 3]);
    imagesc(VD_Wb_self);
    title(['self WB (alpha: ',num2str(Wb_VD_auto(1),'%1.2f'),', ',num2str(Wb_VD_auto(2),'%1.2f'),', ',num2str(Wb_VD_auto(3),'%1.2f'),')']);
    set(gca, 'XTick',[],'YTick',[]);
    sgtitle(sgtitle_name,'Interpreter','none');
    set(gcf,'Color','w');
    set(gcf,'Position',[300,50,1000,990]);
    img_file_name=[sgtitle_name,'.jpg'];
    exportgraphics(gcf,img_file_name,'Resolution',500)
    
    % brightness stretching
    max_VD_RGB(i)=max(max(VD_RGB));
    if strcmpi(Max_normal,'Y')
        VD_none=VD_RGB./0.8;
        VD_Wb_VC_export=VD_Wb_VC./0.8;
        Band_i=Band_i./0.04;
    else
        VD_none=VD_RGB./max_VD_RGB(i);
        VD_Wb_VC_export=VD_Wb_VC./max_VD_Wb_VC(i);
        Band_i=Band_i./max(max(Band_i));
    end
    
    %save the images
    VD_u16_1=uint16(VD_Wb_VC_export*65535);
    imwrite(VD_u16_1,[sgtitle_name,'-Wb_VC.tiff']) %after white balance
    VD_none=reshape(VD_none,[256,256,3]);
    VD_none=permute(VD_none,[2 1 3]);
    VD_u16_2=uint16(VD_none*65535);
    imwrite(VD_u16_2,[sgtitle_name,'-Wb_none.tiff']) %none white balance
    VD_u16_3=uint16(Band_i*65535);
    imwrite(VD_u16_3,[sgtitle_name,'-',num2str(CMOS_WL(23)),'nm','.tiff']) %save lossless photo
end
toc
%%
%input：
%WL: wavelength of CMOS hyperspectral data，CE4:450-2395@5nm
%CMOS_dat: CMOS data,65536×100
%Cor_coef: color coefficient of 3-color vision
%WB_mod='none/self/mannual'; none: no white balance; self: white balance
%according to itself，mannual: white balance manually (e.g. coefficients from the calibration panel)
%Wb_in=[Alpha_R,Alpha_G,Alpha_B] white balance coefficients of RGB channel
%output：
%Wb: white balance coefficients [R,G,B], normalize according G channel
%C_RGB: result without white balance
%C_RGB_Wb: result after white balance
function [Wb,C_RGB,C_RGB_Wb]=CMOS2RGB(CMOS_WL,CMOS_dat,Cor_coef,WB_mod,Wb_in) %
Cor_WL=Cor_coef(:,1);
CMOS_SG=smoothdata(CMOS_dat','sgolay',15);%smooth the data because of the high noise level
CMOS_inte=interp1(CMOS_WL(7:100),CMOS_SG(7:100,:),Cor_WL,'linear','extrap');%extrapolate to the short wavelength range; there are artifacts of CE-4 CMOS 450-475nm
CMOS_inte=CMOS_inte';
C_RGB(:,1)=CMOS_inte*Cor_coef(:,2);%Red
C_RGB(:,2)=CMOS_inte*Cor_coef(:,3);%Green
if strcmpi(WB_mod,'none') %no white balance
    C_RGB_Wb=C_RGB;
elseif strcmpi(WB_mod,'self') %auto white balance
    R_avg=mean(C_RGB(:,1));
    G_avg=mean(C_RGB(:,2));
    B_avg=mean(C_RGB(:,3));
    K=(R_avg+G_avg+B_avg)/3;
    Wb=[K/R_avg,K/G_avg,K/B_avg];
    Wb=Wb./Wb(2); %normalize according G channel
    C_RGB_Wb=C_RGB.*Wb;
elseif strcmpi(WB_mod,'mannual') %manual coefficients
    Wb=Wb_in;
    Wb=Wb./Wb(2); %normalize according G channel
    C_RGB_Wb=C_RGB.*Wb;
end
end
