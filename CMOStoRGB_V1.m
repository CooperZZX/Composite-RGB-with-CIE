%% 本程序用于将嫦娥4的CMOS高光谱进行RGB图像合成,对每一个像素进行滤波，选择多样
%请在使用前将月昼数据的文件夹名称改为‘DayXX’的格式，以便程序自动提取月昼参数
%使用前需要更改三个路径参数
%roless滤波开销太大，SG就要快不少
%使用定标板文件（VC）计算白平衡参数，即以太阳为白色，大约等效为5500K色温
%使用方法：
%修改COR_mod以选择采用那种色彩模式合成RGB图像
%版本说明：V1
%最初版本
restoredefaultpath
clear
close all
addpath 'G:\code\GA-PLSR-CE4' %本程序以及所需函数所在路径
addpath G:\code\Other-data
% addpath 'G:\code\CE-data'
output_path='J:\CE4-VNIS-exported\CMOS';

% load CE4_VNIS_par.mat %载入CMOS中SWIR区域的id，VNIS波长，对太阳光谱的响应等参数

RGB_modtb=readtable('CMOS2RGBCorMod.xlsx'); %读取颜色模式
RGB_mod=table2array(RGB_modtb);

% Day_name=char(regexp(filePath,'Day+\d+\','match'));%利用正则表达式匹配文件夹名称中的月昼编号，\d表示任意数字
% Day_name=Day_name(1:end-1);
% Day_num=str2double(Day_name(4:end));
% if Day_num<7
%     warning('第七月昼之前的数据需要乘以校正系数!!!');
% end
% if Day_num<7
%     disp('------Day num less than 7, calibration matrix was used!!!!')
%     Ref_CMOS_0=Ref_CMOS_0.*CMOS_cali';
%     Ref_SWIR_0=Ref_SWIR_0.*SWIR_cali';
% end
%%
%可以控制day_0和day_end，选择需要反演的光谱Day的范围，也可以设置为其中的某一天
filePath='J:\CE4-data\DayXX\VNIS'; %初始化数据文件所在路径
day_0=10;
day_end=10;
Day_normal='N';%Y代表对一个光谱进行亮度归一化，N代表每一张CMOS照片的亮度以自己的最大值进行拉伸
COR_mod='CIE';% 3bands 代表450,540,565nm，Stock 代表Stockman3色视觉，CIE 代表CIE1931
% WB_mod='none/self/mannual'; %none代表无白平衡，self代表以自身自动白平衡，mannual代表手动值（如来自每个月昼的定标板的白平衡数值）
Wb_man=[0.808,0.914,1.498];%默认RGB的白平衡值
CMOS_WL=450:5:945;
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
tic
for day_file_n=day_0:day_end
    if day_file_n<10
        filePath(17)=num2str(day_file_n); %数据文件所在路径
        filePath(16)='0';
    else
        filePath(16:17)=num2str(day_file_n); %数据文件所在路径
    end
cd(filePath);

dir_VC=dir(fullfile(filePath,'*VC_SCI*.2B'));%定标板光谱文件
file_VC={dir_VC.name}';
file_VC=cell2mat(file_VC);   % convert cell to matrix.
fileNum_VC = size(file_VC,1); % count the total number of files.
Wb_VC_all=zeros(fileNum_VC,3);
dir_VD=dir(fullfile(filePath,'*VD_SCI*.2B'));%月壤光谱文件
file_VD={dir_VD.name}';
file_VD=cell2mat(file_VD);   % convert cell to matrix.
fileNum_VD = size(file_VD,1); % count the total number of files.
Wb_VD_all=zeros(fileNum_VD,3);
cd(output_path)
figure();
for i=1:fileNum_VC
    basename=file_VC(i,:);
    sgtitle_name=['N',basename(end-7:end-5),'-D',num2str(day_file_n),'-VC-',COR_mod_name];
    filename_VC=[filePath,'\',basename];
    fID_VC = fopen(filename_VC);
    VC_sor0 = fread(fID_VC,'float')';%读取辐照度数据
    fclose(fID_VC);
    VC_sor0=reshape(VC_sor0,[65536,100]);%
    [Wb_VC,VC_RGB,VC_RGB_Wb]=CMOS2RGB(CMOS_WL,VC_sor0,Cor_ind,'self');%self代表以自身自动白平衡
    VC_sor=reshape(VC_sor0,[256,256,100]);
    VC_sor=permute(VC_sor,[2 1 3]);
    Rad_1band=VC_sor(:,:,23);
    Wb_VC_all(i,:)=Wb_VC;
    M_VC=mean(VC_sor0(:,23));
    STD_VC=std(VC_sor0(:,23));
    if M_VC<0.05%如果定标板上无阳光，则抛弃这一个定标板文件，0.05来自经验参数
        disp([filename_VC,' is an dark VC measurement']);
        Wb_VC_all(i,:)=[];
    end
    if STD_VC>0.05%如果定标板上部分面积未被照射到，则抛弃这一个定标板文件，0.05来自经验参数
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
    imagesc(Rad_1band);%第23波段是560nm
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
%figure();
% colormap gray ;
% for k=1:100
%     imagesc(VC_sor(:,:,k));
%     pause(0.2);
% end %播放每一个波段的画面 可以看到全局的坏点和800nm之后的左右不均匀性
end
%%
if fileNum_VC==0
    disp([filePath,' no VC file, Wb_VD will be manually assigned']);
    Wb_VD=[1.01327568727212,1,0.924242434049424];
elseif isempty(Wb_VC_all) %如果没有定标板数据
    disp([filePath,' invalid VC file, Wb_VD will be manually assigned']);
    Wb_VD=[1.01327568727212,1,0.924242434049424];
else
    Wb_VD=mean(Wb_VC_all,1); %使用定标板数据得到的平均白平衡参数
end
figure();
for i=1:fileNum_VD
    basename=file_VD(i,:);
    sgtitle_name=['N',basename(end-7:end-5),'-D',num2str(day_file_n),'-VD-',COR_mod_name];
    filename_VD=[filePath,'\',basename];
    fID_VD = fopen(filename_VD);
    VD_sor0 = fread(fID_VD,'float')';
    fclose(fID_VD);
    VD_sor0=reshape(VD_sor0,[65536,100]);%
    [Wb_VD_out,VD_RGB,VD_Wb_VC]=CMOS2RGB(CMOS_WL,VD_sor0,Cor_ind,'mannual',Wb_VD);%mannual代表手动值（如来自每个月昼的定标板的白平衡数值）
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
    imagesc(Band_i);%第23波段是560nm
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
    
    %以自身进行自动白平衡
    subplot(2,2,4)
    [Wb_VD_auto,VD_RGB,VD_Wb_self]=CMOS2RGB(CMOS_WL,VD_sor0,Cor_ind,'self');%[Wb,C_RGB,C_RGB_Wb]=CMOS2RGB(CMOS_WL,CMOS_dat,Cor_ind,WB_mod,Wb_in),'self'代表以自身自动白平衡
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
    
    %亮度拉伸
    max_VD_RGB(i)=max(max(VD_RGB));
    if strcmpi(Day_normal,'Y')
        VD_none=VD_RGB./0.8;
        VD_Wb_VC_export=VD_Wb_VC./0.8;
        Band_i=Band_i./0.04;
    else
        VD_none=VD_RGB./max_VD_RGB(i);
        VD_Wb_VC_export=VD_Wb_VC./max_VD_Wb_VC(i);
        Band_i=Band_i./max(max(Band_i));
    end
    
    %存储CMOS的白平衡和非白平衡tiff图像
    VD_u16_1=uint16(VD_Wb_VC_export*65535);
    imwrite(VD_u16_1,[sgtitle_name,'-Wb_VC.tiff']) %存储白平衡后的图像
    VD_none=reshape(VD_none,[256,256,3]);
    VD_none=permute(VD_none,[2 1 3]);
    VD_u16_2=uint16(VD_none*65535);
    imwrite(VD_u16_2,[sgtitle_name,'-Wb_none.tiff']) %存储非白平衡的图像
    VD_u16_3=uint16(Band_i*65535);
    imwrite(VD_u16_3,[sgtitle_name,'-',num2str(CMOS_WL(23)),'nm','.tiff']) %save lossless photo
%figure();
% colormap gray ;
% for k=1:100
%     imagesc(VD_sor(:,:,k));
%     pause(0.2);
% end %播放每一个波段的画面 可以看到全局的坏点和800nm之后的左右不均匀性
end
end
toc


%输入：
%WL cmos的波长，CE4是450-2395@5nm
%CMOS_dat CMOS的数据，65536×100
%Cor_mod=[1 2 3]; 1代表3bands，2代表3色视觉，3代表CIE1931
%WB_mod='none/self/mannual'; none代表无白平衡，self代表以自身自动白平衡，mannual代表手动值（如来自每个月昼的定标板的白平衡数值）
%Wb_in=[Alpha_R,Alpha_G,Alpha_B]是手动设定的三个波段的白平衡值
%输出：
%Wb [R,G,B]的三个白平衡参数，以G为1
%C_RGB 未进行白平衡的结果
%C_RGB_Wb 进行白平衡之后的结果
function [Wb,C_RGB,C_RGB_Wb]=CMOS2RGB(CMOS_WL,CMOS_dat,Cor_ind,WB_mod,Wb_in) %
Cor_WL=Cor_ind(:,1);
CMOS_SG=smoothdata(CMOS_dat','sgolay',15);%由于波动噪声特别大，所以先进行一次滤波
CMOS_inte=interp1(CMOS_WL(7:100),CMOS_SG(7:100,:),Cor_WL,'linear','extrap');%将CMOS的波段插值至颜色模式的波段（短波方向需要外插，长波方向CMOS是冗余的）'extrap'是采用与内插相同的策略。CE-4的450-475波段有问题
CMOS_inte=CMOS_inte';
% figure();
% plot(CMOS_WL,CMOS_dat(1:5000:end,:),'Color',[0.7 0.7 0.7]);
% hold on
% CMOS_SG=CMOS_SG';
% plot(CMOS_WL,CMOS_SG(1:5000:end,:),'Color',[0.1 0.1 0.9]);
% plot(Cor_WL,CMOS_inte(1:5000:end,:),'Color',[0.9 0.4 0.4]);
% plot(Cor_WL,Cor_ind(:,2:4));
C_RGB(:,1)=CMOS_inte*Cor_ind(:,2);%红
C_RGB(:,2)=CMOS_inte*Cor_ind(:,3);%绿
C_RGB(:,3)=CMOS_inte*Cor_ind(:,4);%蓝
% CMOS_Eye_i=CMOS_Eye_i/max(max(CMOS_Eye_i));%归一化
% CMOS_Eye=reshape(CMOS_Eye_i,[256,256,3]);
if strcmpi(WB_mod,'none') %无白平衡
    C_RGB_Wb=C_RGB;
elseif strcmpi(WB_mod,'self') %自动白平衡
    R_avg=mean(C_RGB(:,1));
    G_avg=mean(C_RGB(:,2));
    B_avg=mean(C_RGB(:,3));
    K=(R_avg+G_avg+B_avg)/3;
    Wb=[K/R_avg,K/G_avg,K/B_avg];
    Wb=Wb./Wb(2);%以G波段为1
    C_RGB_Wb=C_RGB.*Wb;
elseif strcmpi(WB_mod,'mannual')%手动值
    Wb=Wb_in;
    Wb=Wb./Wb(2);%以G波段为1
    C_RGB_Wb=C_RGB.*Wb;
end

end

% CMOS_ant_sele=CMOS_sor(:,CMOS_show_band);
% CMOS_ant_sele(SWIR_ind_logi)=0;
% CMOS_ant_sele=reshape(CMOS_ant_sele,[256,256])';
% imagesc(CMOS_ant_sele,'AlphaData',0.6)
% images.roi.Circle(gca,'Center',[98 256-128],'Radius',53.8);

% system('rename *SD*.2B *.txt');%将SD/SWIR.2B文件更改后缀为txt，以便以readtable方法快速读取数据
% system('rename *SD*.txt *.2B');%将更改的文件后缀还原
%     [time, bands, exposure, ang, radi, qua] = textread(filename_SD2B,'%24s %2d %f %s %f %f');%会多读取空格为0值进来