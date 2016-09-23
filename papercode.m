function papercode()

% == Step1 Introduction ==
addpath('GUI');
addpath('MATLAB');
addpath('Confocal');
addpath('[ Original code ]');
warning('off', 'Images:initSize:adjustingMag');
p=p_setting();


% == Step2 Image conversion and pre-analysis processing ==
% Once done, you will have a subfolder called ?stacks? inside your data folder (e.g. ?Exp1?) containing all TIFF files, and another subfolder ?Results? for future analysis results.
if p.lsmconvert_run==1;lsmconvert(p);end
if p.flyimagetile3_run==1;flyimagetile3(p);end
if p.Dualprocess_run==1;Dualprocess(p);end


% == Step3. Nuclear segmentation ==
% To reconstruct the 3D shape of individual nuclei in the embryo, run the following scripts:,All segmentation results are stored in the subfolder ?masks? inside your data folder
if p.DAPI_seg3D2_run==1;DAPI_seg3D2(p);end        % For automatic segmentation
if p.nucle_manu3D_run==1;nuclei_manual3D;end        % For manual refinement of the segmentation result, run the MATLAB GUI 


% == Step4. smFISH spot analysis ==
% To identify smFISH spots and quantify their intensities, run the following scripts, Results of nascent and cytoplasmic mRNAs are stored in the subfolders ?Histogram? and?Histogram_A? respectively, inside your data folder.
%stack_RNA;      % For automatic spot quantification, Nascent and cytoplasmic mRNAs need to be identified separately by choosing different ?Analysis types? (?foci? or ?single?) in the GUI.
if p.stack_RNA_GUIfun_run==1;
    p.type_foci=1;stack_RNA_GUIfun(p);p.type_foci=0;
    p.type_single=1;stack_RNA_GUIfun(p);p.type_single=0;
end     % Re-write from stack_RNA
load matlab.mat;p.hist_on_GUIfun_run = 1;
if p.hist_on_GUIfun_run==1;hist_on_GUIfun(p,imname,imfolder,imstack,imstack0,mask_stack,bk_use,mask_out);end
if p.hist_save_GUIfun_run==1;hist_save_GUIfun(p);end
if p.hist_fit_run==1;hist_fit('Duallist.mat');end    % To extract the typical intensity of a single mRNA
if p.stack_RNA_check_run==1;stack_RNA_check;end        % To double-check and improve the identification of nascent mRNAs


% == Step5. IF spot analysis ==
% Results for anterior and posterior IF spots are stored in the subfolders ?Histogram_protein_A? and ?Histogram_protein_P? respectively, inside your data folder.
if p.stack_protei_run==1;stack_protein;end      % To identify IF spots, run the MATLAB GUI, Anterior and posterior spots need to be identified separately by choosing different ?Analysis types? (?Anterior? or ?Posterior?) in the GUI. 


% == Step6. Analysis of transcriptional regulation ==
% After completing steps 1-5, run ?Dualprocess6;? to analyze transcriptional regulation. The key results are saved as a MAT file, with the same name as your image file, and stored in the subfolder ?Results? inside your data folder.
if p.Dualprocess6_run==1;Dualprocess6;end

% == Step7. Fitting the nascent mRNA distribution ==
% To fit the nascent mRNA distribution to a two-state model, go to the folder ?FISHIF/simu_FISH?, and use MATLAB functions ?NSX_fitbinN?, ?NSX_fit3binN?, ?NSX_scan3binN?. Comments in these m-files provide detailed description of the inputs and outputs of the functions.
%NSX_fitbinN;
%NSX_fit3binN;
%NSX_scan3binN;

end

function p=p_setting()

%step2
p.outtiff=1;
p.lsmconvert_run = 0;
p.flyimagetile3_run = 0;
p.Dualprocess_run = 0;

% step 3
p.DAPI_seg3D2_run = 0;
p.nucle_manu3D_run = 0;

p.type_foci = 0;        % RNA channel
p.type_single = 0;
p.type_foci2 = 0;       % Signal2 channel
p.type_single2 = 0;


% step 4
p.stack_RNA_GUIfun_run = 0;     % run stack_RNA_GUIfun
p.hist_on_GUIfun_run = 1;       % run hist_on_GUIfun(p,imname,imfolder,imstack,imstack0,mask_stack,bk_use,mask_out);
p.hist_save_GUIfun_run = 0;     % run hist_save_GUIfun(p)
p.hist_fit_run = 0;             % run hist_fit
p.stack_RNA_check_run = 0;      % run stack_RNA_check GUI
p.presult_on = 0;
p.show_bk = 0;
p.background_on = 0;

% step5
p.stack_protei_run = 0;

% step 6
p.Dualprocess6_run = 0;


end

%% ======= Step1 Introduction =======
%% ======= Step2 Image conversion and pre-analysis processing =====
function lsmconvert(p)
%clear all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lf_name = 'Confocal/Duallist.mat';
out_folder = 'stacks/';
lsm_type = '*.lsm';
tif_name = 'stack';
figure_tail = '.tif';
match_file = 'matchlist.mat';
match_key = 'match';
compare_ratio = 1;

if ispc==1
    out_folder(findstr(out_folder, '/'))='\';
    lf_name(findstr(lf_name, '/'))='\';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.mat')
    list_name = lf_name;
    %[num_list, folder_list] = xlsread(list_name); 
    load(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
    % ispc=1 then use'\'
    if ispc==1
        folder_list{1}(findstr(folder_list{1}, '/'))='\';
    end
    
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);

for list_I = 1:N1
    input_folder = folder_list{list_I,1};
    if isempty(folder_list{list_I,2})
        match_list = [];
    else
        match_list = eval(folder_list{list_I,2});
    end
    all_channel = eval(folder_list{list_I,3});
    DAPI_channel = all_channel(1);
    WGA_channel = all_channel(2);
    signal1_channel = all_channel(3);
    signal2_channel = all_channel(4);
    if length(all_channel) > 4
        signal3_channel = all_channel(5:end);
    else
        signal3_channel = zeros(1);
    end
    match_channel = WGA_channel;
    all_other = eval(folder_list{list_I,4});
    Nbin = all_other(1);
    Mdim = all_other(2);
    
%% LSM file loading/resave: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lsm_name = dir([input_folder,lsm_type]);
    if exist([input_folder,out_folder]) ~= 7
        mkdir([input_folder,out_folder]);
    end

    file_name = cell(length(lsm_name),2);
    file_num = zeros(length(lsm_name),12);
    output_I = 0;
    output_switch = false;

    for I_file = 1:length(lsm_name)
        input_name = lsm_name(I_file).name;
        lsm_stack = tiffread([input_folder,input_name]); %%% Lsm file loading
        pack_num = 1;
        stack_num = 0;
        bin_num = 0;
        if  ~isempty(strfind(input_name,match_key))
            bin_size = Nbin*2-1+(Nbin == 1);
        else
            bin_size = 1;
        end
        if ~(ismember(input_name(find(input_name == '.',1,'last')-1),'Bb') && input_name(find(input_name == '.',1,'last')-2) == '_')
            output_I = output_I+1;
        end
        
        while stack_num(end) < length(lsm_stack)
            bin_num = bin_num+1;
            pack_num = pack_num+(bin_num > bin_size);
            output_I = output_I+(bin_num > bin_size);
            bin_num = mod(bin_num,bin_size);
            if bin_num == 0
                bin_num = bin_size;
            end
            stack_num = stack_num(end)+[1:lsm_stack(stack_num(end)+1).lsm.DimensionZ];
            
    %%% Output folder setup: %%% ------------------------------------------
            if  bin_size*lsm_stack(1).lsm.DimensionZ >= length(lsm_stack)
                pack_name = '';
            else
                pack_name = ['_',num2str(pack_num,'%03u')];
            end
            if ~isempty(strfind(input_name,match_key))
                match_name = ['_',char('B'-mod(bin_num,2))];
            else
                match_name = '';
            end
            output_name = [input_name(1:(find(input_name == '.',1,'last')-1)),pack_name,match_name,'/'];
            if exist([input_folder,out_folder,output_name]) ~= 7
                mkdir([input_folder,out_folder,output_name]);
            end
    %%% -------------------------------------------------------------------
    
            for I_layer = stack_num
                stack_raw = lsm_stack(I_layer).data;
                tiff_image = zeros(0);
                if iscell(stack_raw)
                    for I_color = 1:length(stack_raw)
                        tiff_image = cat(3,tiff_image,stack_raw{I_color});
                    end
                else
                    tiff_image = stack_raw;
                end
                
                if size(tiff_image,3) == 1
                    tiff_image(:,:,2) = tiff_image(:,:,1);
                end
                if size(tiff_image,3) == 2
                    tiff_image(:,:,3) = tiff_image(:,:,2);
                end
    %%% Image output: %%% -------------------------------------------------
                if length(lsm_stack) > 1
                    out_num = num2str(I_layer-stack_num(1)+1,'%02u');
                else
                    out_num = '';
                end
                out_stack = [tif_name,out_num,figure_tail];
                if p.outtiff==1
                    if bin_num <= 2
                        imwrite(tiff_image,[input_folder,out_folder,output_name,out_stack]);
                    else
                        temp_old = imread([input_folder,out_folder,output_name,out_stack]);
                        tiff_image = cat(Mdim,temp_old,tiff_image);
                        imwrite(tiff_image,[input_folder,out_folder,output_name,out_stack]);
                    end
                end
    %%% -------------------------------------------------------------------
            end
    
%%% matchlist.xls output: %%%==============================================
            if ismember(output_name(length(output_name)-1),'Bb') && output_name(length(output_name)-2) == '_'
                output_switch = true;
            else
                output_switch = false;
            end
            file_name{output_I,output_switch+1} = output_name;
            %match_layer = ceil(length(lsm_stack)/2);
            if isempty(match_list)
                match_layer = 1;
            else
                match_layer = match_list(output_I);
            end
            resolution = lsm_stack(1).lsm.VoxelSizeX/(1e-6);
            resolutionz = lsm_stack(1).lsm.VoxelSizeZ/(1e-6);
            file_num(output_I,:) = [match_layer,Nbin,Mdim,match_channel,compare_ratio,WGA_channel,DAPI_channel,signal1_channel,resolution,signal2_channel,resolutionz,signal3_channel];
%%% =======================================================================
        end
    end

%file_num(strmatch(else_name,{lsm_name.name}),:) = [match_layer,eNbin,eMdim,eWGA_channel,compare_ratio,eWGA_channel,eDAPI_channel,esignal_channel,resolution];    
    try
        %xlswrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
        match_folder=file_name;match_num=file_num;
        save([input_folder,out_folder,match_file],'match_folder','match_num');
    catch
        %xlwrite([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
        save([input_folder,out_folder,match_file],cat(2,file_name,num2cell(file_num)));
    end
end
end
function flyimagetile3(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fly embyro image tiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lf_name = 'Confocal/Duallist.mat';
list_sub = 'stacks/';
single_add = '_new/';
list_name = 'matchlist.mat';
standard_record = 'Calibration/Results/standard.mat';
standard_record_60X = 'Calibration3/Results/standard_60X.mat';
alternative_zoom = '_60X';
output_tail = '.xls';
image_tail = '*.tif';
figure_tail = '.fig';
th_parameter1 = 1.2;
th_parameter2 = 1.2;
if ispc==1
    lf_name(findstr(lf_name, '/'))='\';
    list_sub(findstr(list_sub, '/'))='\';
    single_add(findstr(single_add, '/'))='\';
    standard_record(findstr(standard_record, '/'))='\';
    standard_record_60X(findstr(standard_record_60X, '/'))='\';    
end

global standard_data standard_data_60X
standard_data = load(standard_record);
standard_data_60X = load(standard_record_60X);

if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.mat')
    list_name0 = lf_name;
    %[num_list, folder_list] = xlsread(list_name0);
    load(list_name0);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
    % ispc=1 then use'\'
    if ispc==1
        folder_list{1}(findstr(folder_list{1}, '/'))='\';
    end
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[NN1,NN2] = size(folder_list);

for list_I0 = 1:NN1
    folder_name = [folder_list{list_I0,1},list_sub];
    %channel_name = eval(folder_list{list_I0,5});
    channel_name = folder_list{list_I0,5};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[num_list, file_list] = xlsread([folder_name,list_name]);
    load([folder_name,list_name]);
    num_list = match_num;
    file_list = match_folder;
    if ispc==1
        file_list{1}(findstr(file_list{1}, '/'))='\';
    end
    
    [N1,N2] = size(file_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for list_I = 1:N1
        resolution = num_list(list_I,9);
        if ~isempty(strfind(file_list{list_I,1},alternative_zoom))
            change_zoom = true;
        else
            change_zoom = false;
        end
        if size(file_list,2) > 1 && (~isempty(file_list{list_I,2}))
    %%% Image list loading: %%%================================================
            imlist1 = dir([folder_name,file_list{list_I,1},image_tail]); %%% get the image list from image folder1
            imlist2 = dir([folder_name,file_list{list_I,2},image_tail]); %%% get the image list from image folder2
            out_name = [file_list{list_I,2}(1:(end-3)),'/']; %%% output folder name
            if length(imlist1) ~= length(imlist2)
                error(['Incompatible image stacks: ',folder_name,file_list{list_I,1},' and ',folder_name,file_list{list_I,2}])
            elseif length(imlist1) < num_list(list_I,1)
                error(['Not enough image layers (',folder_name,file_list{list_I,1},'): ',num2str(num_list(list_I,1)),'/',num2str(length(imlist1))])
            end
        %%% =======================================================================

    %%% Initialization of image test matching process: %%%=====================
            Mdim = num_list(list_I,3);
            match_channel = num_list(list_I,4);
            compare_ratio = num_list(list_I,5);
            prematch = 0;
    %%% =======================================================================

    %%% Image loading: %%%=====================================================
            image1 = imread([folder_name,file_list{list_I,1},imlist1(num_list(list_I,1)).name]);
            image2 = imread([folder_name,file_list{list_I,2},imlist2(num_list(list_I,1)).name]);
            Nbin1 = num_list(list_I,2);
            pixel_bin1 = size(image1,Mdim)/Nbin1;
            Nbin2 = Nbin1-1+(Nbin1 == 1);
            pixel_bin2 = size(image2,Mdim)/Nbin2;
        %%%%% Index setup/initialization: %%%%%--------------------------------
            IImatch = zeros(1,(Nbin1+Nbin2-1));
            input1_start = [1,1];
            input1_end = size(image1);
            input1_end(Mdim) = pixel_bin1;
            input2_start = [1,1];
            input2_end = size(image2);
            input2_end(Mdim) = pixel_bin2;
            outimage = image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:);
        %%%%% Intensity correction mask generation: ---------------------------
            corr1_size = size(image1);
            corr1_size(Mdim) = pixel_bin1;
            corr2_size = size(image2);
            corr2_size(Mdim) = pixel_bin2;
            imcorr1 = corr_mask(corr1_size,channel_name,resolution,change_zoom);
            imcorr2 = corr_mask(corr2_size,channel_name,resolution,change_zoom);
        %%%%% -----------------------------------------------------------------
            
        %%%%% Matching parameter search: %%%%%---------------------------------
            Ibin1 = 1;
            Ibin2 = 1;
            I_tile = 1;
            while (Ibin1 < Nbin1)||(Ibin2 <= Nbin2)
                if Ibin2 <= Nbin2
                    input2 = image2(input2_start(1):input2_end(1),input2_start(2):input2_end(2),:);
                    [outimage,IImatch(I_tile),mismatch] = tilematch(outimage,input2,Mdim,match_channel,prematch,compare_ratio);
                    I_tile = I_tile+1;
                    Ibin2 = Ibin2+1;
                    input2_start(Mdim) = input2_start(Mdim)+pixel_bin2;
                    input2_end(Mdim) = input2_end(Mdim)+pixel_bin2;
                end

                if Ibin1 < Nbin1
                    input1_start(Mdim) = input1_start(Mdim)+pixel_bin1;
                    input1_end(Mdim) = input1_end(Mdim)+pixel_bin1;
                    input2 = image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:);
                    [outimage,IImatch(I_tile),mismatch] = tilematch(outimage,input2,Mdim,match_channel,prematch,compare_ratio);
                    I_tile = I_tile+1;
                    Ibin1 = Ibin1+1;
                end
            end
            %IImatch
        %%%%% -----------------------------------------------------------------
    %%% =======================================================================

    %%% Tiling loops: %%%%%====================================================
            for I_layer = 1:length(imlist1)
                image1 = imread([folder_name,file_list{list_I,1},imlist1(I_layer).name]);
                image2 = imread([folder_name,file_list{list_I,2},imlist2(I_layer).name]);
                %%%%% Initialization:
                input1_start = [1,1];
                input1_end = size(image1);
                input1_end(Mdim) = pixel_bin1;
                input2_start = [1,1];
                input2_end = size(image2);
                input2_end(Mdim) = pixel_bin2;
                outimage = uint16(double(image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:))./imcorr1);

                %%%%% Tiling:
                Ibin1 = 1;
                Ibin2 = 1;
                I_tile = 1;
                while (Ibin1 < Nbin1)||(Ibin2 <= Nbin2)
                    if Ibin2 <= Nbin2
                        input2 = uint16(double(image2(input2_start(1):input2_end(1),input2_start(2):input2_end(2),:))./imcorr2);
                        [outimage,IImatch(I_tile),mismatch] = tilematch2(outimage,input2,Mdim,match_channel,IImatch(I_tile),compare_ratio);
                        I_tile = I_tile+1;
                        Ibin2 = Ibin2+1;
                        input2_start(Mdim) = input2_start(Mdim)+pixel_bin2;
                        input2_end(Mdim) = input2_end(Mdim)+pixel_bin2;
                    end

                    if Ibin1 < Nbin1
                        input1_start(Mdim) = input1_start(Mdim)+pixel_bin1;
                        input1_end(Mdim) = input1_end(Mdim)+pixel_bin1;
                        input2 = uint16(double(image1(input1_start(1):input1_end(1),input1_start(2):input1_end(2),:))./imcorr1);
                        [outimage,IImatch(I_tile),mismatch] = tilematch2(outimage,input2,Mdim,match_channel,IImatch(I_tile),compare_ratio);
                        I_tile = I_tile+1;
                        Ibin1 = Ibin1+1;
                    end
                end

                %%%%% Image output:
                %imshow(outimage);
                %title([folder_name,out_name,': ',num2str(I_layer),'/',num2str(length(imlist1))])
                %waitforbuttonpress
                if exist([folder_name,out_name]) ~= 7
                    mkdir([folder_name,out_name]);
                end
                imwrite(outimage,[folder_name,out_name,imlist1(I_layer).name]);
            end
    %%% =======================================================================

        else
            out_name = [file_list{list_I,1}(1:(end-1)),single_add]; %%% output folder name
            imlist1 = dir([folder_name,file_list{list_I,1},image_tail]); %%% get the image list from image folder1
            
            %%% Tiling loops: %%%%%====================================================
            for I_layer = 1:length(imlist1)
                image1 = imread([folder_name,file_list{list_I,1},imlist1(I_layer).name]);
                Nbin1 = num_list(list_I,2);
                Mdim = num_list(list_I,3);
                pixel_bin1 = size(image1,Mdim)/Nbin1;
                
            %%%%% Intensity correction mask generation: ---------------------------
                corr1_size = size(image1);
                corr1_size(Mdim) = pixel_bin1;
                imcorr1 = corr_mask(corr1_size,channel_name,resolution,change_zoom);
            %%%%% -----------------------------------------------------------------
                
                %%%%% Initialization:
                rep_size = [1,1];
                rep_size(Mdim) = Nbin1;
                outimage = uint16(double(image1)./repmat(imcorr1,rep_size));
                if exist([folder_name,out_name]) ~= 7
                    mkdir([folder_name,out_name]);
                end
                imwrite(outimage,[folder_name,out_name,imlist1(I_layer).name]);
            end
            %%% =======================================================================
            
        end

        file_list(list_I,3) = {out_name}; %%% Output folder name record
    end
    
    try
        %xlswrite([folder_name,list_name],cat(2,file_list,num2cell(num_list)));
        match_folder=file_list;match_num=num_list;
        save([folder_name,list_name],'match_folder','match_num');
    catch
        xlwrite([folder_name,list_name],cat(2,file_list,num2cell(num_list)));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

clear global
end
function Dualprocess(p)
clear all

%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Confocal/Duallist.mat';
in_folder = 'stacks/';
input_name = 'matchlist.mat';
standard_record = 'Calibration/Results/standard.mat';
image_type = '*.tif';
out_folder = 'Results/';
output_tail = '.mat';
figure_tail = '.fig';
mat_tail = '.mat';
seg_add = '_seg';
th_add = '_th';
foci_add = '_foci';
fish_add = '_fish';
num_add = '_num';
int_add = '_int';
cmp_add = '_cmp';
sel_add = '_sel';
protein_add = '_protein';
RNA_add = '_RNA';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
reg_add = '_regulation';
immu_add = '_immu';
fluc_add = '_fluc';
bino_add = '_bino';
emmask_add = '_emmask';
if ispc==1
    list_name(findstr(list_name, '/'))='\';
    in_folder(findstr(in_folder, '/'))='\';
    out_folder(findstr(out_folder, '/'))='\';
    standard_record(findstr(standard_record, '/'))='\';  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[num_list, folder_list] = xlsread(list_name);
    load(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
    % ispc=1 then use'\'
    if ispc==1
        folder_list{1}(findstr(folder_list{1}, '/'))='\';
    end


[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    %[sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    load([folder_list{list_I,1},in_folder,input_name]);
    sub_num = match_num;
    sub_list = match_folder;  
    
    [M1,M2] = size(sub_list);
    channel_name = folder_list{list_I,5};
    for list_J = 1:M1
        %image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,1}];  %George
        Mdim = sub_num(list_J,3);
        Nbin1 = sub_num(list_J,2);
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        RNA_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        signal2_channel = sub_num(list_J,12);
        resolution = sub_num(list_J,9);
        resolutionz = sub_num(list_J,11);
        info.WGA_channel=WGA_channel;
        info.DAPI_channel=DAPI_channel;
        info.RNA_channel=RNA_channel;
        info.protein_channel=protein_channel;
        info.signal2_channel=signal2_channel;

        %[seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_seg4(image_folder,image_type,WGA_channel,DAPI_channel,resolution);
        [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution);
        [mask_stack,signal_stack,~] = mask3D(seg_bw,protein_channel,DAPI_channel,RNA_channel,image_folder);
        em_mask = get_emmask(image_folder);
        figure(1)
            out_image = double(repmat(bwperim(em_mask),[1,1,3]));
            out_image(:,:,3) = out_image(:,:,3)+double(max_image(:,:,1))/max(max(double(max_image(:,:,1))));
            imshow(out_image)
            title(image_folder,'Interpreter','none')
            
%         max_size0 = size(max_image);
%         max_size0(Mdim) = size(max_image,Mdim)/Nbin1;
%         imcorr1 = corr_mask(max_size0,channel_name,resolution);
%         Ntile = ones(size(max_size0));
%         Ntile(Mdim) = Nbin1;
%         imcorr = repmat(imcorr1,Ntile);
%         max_image = uint16(double(max_image)./imcorr);
        
        [foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,RNA_channel,WGA_channel,resolution,image_folder,N_cycle);
        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,RNA_channel,resolution,image_folder,N_cycle);
%         [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile(seg_bw,cyto_bw,max_image,signal_channel,image_folder,N_cycle)
        [nucleus_protein_profile,cytoplasmic_protein_profile] = protein_profile3(seg_bw,cyto_bw,max_image,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution);
        clear mask_stack signal_stack
        dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle);

        figure(16)
            temp_image = zeros(size(max_image,1),size(max_image,2),3);
            temp_image(:,:,2) = double(max_image(:,:,protein_channel))./max(max(double(max_image(:,:,protein_channel))));
            imshow(temp_image);
            title(['Immunofluorescence signal: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none');
        figure(17)
            imshow(double(max_image(:,:,RNA_channel))./max(max(double(max_image(:,:,RNA_channel)))));
            title(['FISH signal: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none');
        
        
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        if ~exist([result_folder,sub_list{list_J,3}],'dir')
            mkdir([result_folder,sub_list{list_J,3}])
        end
        saveas(1,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),emmask_add,figure_tail]);
        saveas(51,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,figure_tail]);
        saveas(3,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
        saveas(4,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,cmp_add,figure_tail]);
        saveas(10,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,figure_tail]);
        saveas(12,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,figure_tail]);
        saveas(13,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,figure_tail]);
        saveas(15,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),immu_add,figure_tail]);
        saveas(17,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,figure_tail]);
        saveas(62,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,figure_tail]);
        save([result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','signal2_channel','resolution','seg_bw','cyto_bw','em_mask','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','psize0','g_threshold0','WGA_th0','DAPI_th0','resolutionz','info');
        
        if ~isempty(nucleus_RNA_profile)
            %xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
            save([result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],'nucleus_RNA_profile');
        end
        if ~isempty(cytoplasmic_RNA_profile)
            %xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
            save([result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],'cytoplasmic_RNA_profile');
        end
        if ~isempty(foci_RNA_profile)
            %xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
            save([result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],'foci_RNA_profile');
        end
        if ~isempty(nucleus_protein_profile)
            %xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],nucleus_protein_profile);
            save([result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],'nucleus_protein_profile');
        end
        if ~isempty(cytoplasmic_protein_profile)
            %xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],cytoplasmic_protein_profile);
            save([result_folder,sub_list{list_J,3},sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],'cytoplasmic_protein_profile');
        end
        sub_num(list_J,13) = N_cycle;
        
        %clear ('em_mask','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','WGA_th0','DAPI_th0');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
    %xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
    save([folder_list{list_I,1},in_folder,input_name],'match_folder','match_num');
end

close all

% clear global standard_data

end
   % sub
    function imcorr = corr_mask(imsize,channel_name,resolution,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to correct the distortion of the objective: %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global standard_data standard_data_60X

if ~isempty(varargin) && varargin{1}
    standard_data0 = standard_data_60X;
else
    standard_data0 = standard_data;
end
    
st_list = fieldnames(standard_data0);
imcorr = zeros(imsize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Channel matching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I_channel = 1:imsize(3)
    channel_match = find(strncmp(channel_name{I_channel},st_list,length(channel_name{I_channel})),1);
    if ~isempty(channel_match)
        eval(['current_standard = standard_data0.',st_list{channel_match},';']);
        
        %%% Fitting parameter collection: %%% =============================
        fit2D = current_standard.fit;
        resolution0 = current_standard.resolution;
        stsize = current_standard.lattice;
        %%% ===============================================================
        
        %%% Image correction mask generation: %%% =========================
        X_start = (stsize(1)*resolution0-imsize(1)*resolution)/2/resolution0+1;
        Y_start = (stsize(2)*resolution0-imsize(2)*resolution)/2/resolution0+1;
        [X,Y] = meshgrid((Y_start+[0:(imsize(2)-1)]*resolution/resolution0),(X_start+[0:(imsize(1)-1)]*resolution/resolution0));
        XY = [reshape(X,imsize(1)*imsize(2),1),reshape(Y,imsize(1)*imsize(2),1)];
        Z1 = polyvaln(fit2D,double(XY));
        imcorr(:,:,I_channel) = reshape(Z1,imsize(1),imsize(2));
        %%% ===============================================================
    else
        imcorr(:,:,I_channel) = ones(imsize(1:2));
    end
end
    end
    function [seg_bw,cyto_bw,max_image,N_cycle,WGA_th0,DAPI_th0] = nuclei_segDAPI(image_folder,image_type,DAPI_channel,resolution)
%clear
%% Segmentation of nuclei image stacks: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function recognize/segment fly embryo nuclei of a given stack of 
%%% images using WGA and DAPI signals, and output the nucleus segmentation 
%%% and cytoplasmic region segmentation in a binary image, as well as a 
%%% maximum projection image of the original stack for further processing.
%%% The results are also plotted on an overlaied image.

%%% seg_bw (image array): a binary image representing the nuclei recognition result.
%%% cyto_bw (image array): a binary image representing the cytoplasm recognition result.
%%% max_image (multi-channel image array): the maximum projection image of the original image stack.
%%% image_folder (string): name of the folder that contains all images in the stack.
%%% image_type (string): type of the image files in the folder (for example, '*.tif').
%%% WGA_channel (integer: 1,2,3): channel # of the WGA signal.
%%% DAPI_channel (integer: 1,2,3): channel # of the DAPI signal.
%%% resolution (real number): pixel size in um.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolution0 = 0.0915;
L_ratio = (resolution/resolution0);
%se = strel('disk',floor(20/L_ratio));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layer_name = dir([image_folder,image_type]);
if length(layer_name) < 1
    error(['Folder error: ',image_folder,'has no image!'])
end

for I_layer = 1:length(layer_name)
    temp_image = imread([image_folder,layer_name(I_layer).name]);   %%% get images from every layer
    %imstack(:,:,I_layer) = imfilter(double(temp_image(:,:,DAPI_channel)),fspecial('gaussian',10,3/L_ratio),'symmetric','conv');
    if ~exist('max_image')
        max_image = temp_image;
    else
        max_image = max(max_image,temp_image);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% imstack = double(imstack)./max(max(max(double(imstack))));   %%% Normalization
% immask = zeros(size(imstack));
% DAPI_th0 = graythresh(imstack)*1.2;
max_DAPI = max_image(:,:,DAPI_channel);
%max_DAPI = imtophat(max_image(:,:,DAPI_channel),strel('disk',floor(10/L_ratio)));
max_DAPI = imfilter(max_DAPI,fspecial('gaussian',10,3/L_ratio),'symmetric','conv');
DAPI_th0 = graythresh(max_DAPI);
WGA_th0 = DAPI_th0;
%% Raw mask calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for I_layer = 1:size(imstack,3)
%     immask2 = im2bw(imstack(:,:,I_layer),DAPI_th0);
%     immask2 = imclose(immask2,strel('disk',round(3/L_ratio)));
%     immask2 = imopen(immask2,strel('disk',round(20/L_ratio)));
%     immask2 = bwareaopen(immask2,round(1000/L_ratio/L_ratio));
%     immask(:,:,I_layer) = immask2;
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2D Sup mask fine segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maskall = logical(max(immask,[],3));
maskall = im2bw(max_DAPI,DAPI_th0*1);
maskall = imopen(maskall,strel('disk',round(3/L_ratio)));
maskall = imclose(maskall,strel('disk',round(3/L_ratio)));
maskall = imopen(maskall,strel('disk',round(20/L_ratio)));
maskall = bwareaopen(maskall,round(1000/L_ratio/L_ratio));

se_mask = regionprops(maskall,'Area');
se_area = [se_mask.Area];
area1 = geomean(se_area(se_area > 1000/L_ratio/L_ratio));   %%% calculate the mean area for a single nucleus
seg_bw = reseg(maskall,area1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cytoplasm region recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
embryo_region = bwconvhull(seg_bw);
cyto_bw = embryo_region & (~seg_bw);
cyto_bw = imerode(cyto_bw, strel('disk',floor(5/L_ratio)));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nuclei cycle estimation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_prop = regionprops(seg_bw,'Area');
N_cycle = round(log2(size(all_prop(:),1)))+2;   %%% Calculate the nuclear cycle number
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(51)
% bw_perim_WGA = bwperim(seg_bw);
% bw_perim_DAPI = bwperim(bw_DAPI);
% overlay = imoverlay(adapthisteq(max_image(:,:,DAPI_channel)), bw_perim_WGA, [0,1,0]);
% overlay = imoverlay(overlay, bw_perim_DAPI, [1,0,0]);
% imshow(overlay)   %%% show the embryo boundary recognition
% title([image_folder,' (green: nucleus recognition, red: DAPI signal), cycle = ',num2str(N_cycle)]);

bw_perim_DAPI = bwperim(seg_bw);
overlay = zeros([size(seg_bw),3]);
overlay(:,:,1) = 0;
overlay(:,:,2) = bw_perim_DAPI;
overlay(:,:,3) = double(max_image(:,:,DAPI_channel))./max(max(double(max_image(:,:,DAPI_channel))));
try
    imshow(overlay)
catch
    imshow_lxt(overlay)
end
title([image_folder,' (green: nucleus recognition, blue: DAPI signal), cycle = ',num2str(N_cycle)],'Interpreter','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear imstack
    end
    function [mask_stack,signal_stack,RNA_stack,varargout] = mask3D(seg_bw,signal_channel,DAPI_channel,RNA_channel,image_folder)

image_type = '*.tif';
imlist = dir([image_folder,image_type]); %%% get the image list from image folder
thresh3 = 0.5;
N_layer = length(imlist);
DAPI_stack = zeros([size(seg_bw),N_layer]);
signal_stack = zeros([size(seg_bw),N_layer]);
mask_stack = zeros([size(seg_bw),N_layer]);
RNA_stack = zeros([size(seg_bw),N_layer]);
z_DAPI = zeros(0);


%%% Image stack loading: %=================================================
for I_layer = 1:N_layer
    temp = imread([image_folder,imlist(I_layer).name]);
    DAPI_stack(:,:,I_layer) = temp(:,:,DAPI_channel);
    signal_stack(:,:,I_layer) = temp(:,:,signal_channel);
    RNA_stack(:,:,I_layer) = temp(:,:,RNA_channel);
    STATS = regionprops(seg_bw,DAPI_stack(:,:,I_layer),'MeanIntensity');
    z_DAPI(I_layer,:) = double([STATS.MeanIntensity]);
end
%%% =======================================================================

%%% 3D mask generation: ===================================================
z_DAPI = z_DAPI./repmat(max(z_DAPI,[],1),size(z_DAPI,1),1);
region3 = z_DAPI >= thresh3;
region3 = cat(2,zeros(N_layer,1),region3);
bwn = bwlabel(seg_bw);

for I_layer = 1:N_layer
    temp_DAPI = region3(I_layer,:);
    mask_stack(:,:,I_layer) = temp_DAPI(bwn+1);
end
%%% =======================================================================

varargout = {DAPI_stack};
    end
    function em_mask = get_emmask(image_folder,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to extract embryo area mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image_folder (input): embryo image folder;
%% em_mask (output): embryo area mask;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_tail = 'stack*.tif';
fname = dir([image_folder,im_tail]);
im0 = imread([image_folder,fname(end).name]);
if isempty(varargin)
    th0 = 0.6;
else
    th0 = varargin{1};
end
% % % im0 = max_image;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % mask0 = im2bw(im0(:,:,1),0.1*graythresh(im0(:,:,1)));
im00 = double(im0(:,:,1))/max(max(double(im0(:,:,1))));
im00(im2bw(im00,graythresh(im00))) = graythresh(im00);
mask0 = im2bw(im00,th0*graythresh(im00));
% mask0 = imopen(mask0,strel('disk',10));
mask1 = bwareaopen(mask0,1000);
mask1 = imfill(mask1,'holes');
mask1 = imopen(mask1,strel('disk',5));
prop0 = regionprops(mask1,'Area');
[~,Imax] = max([prop0.Area]);
em_mask = bwlabel(mask1) == Imax;
% em_mask = bwconvhull(em_mask);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    function EL_info = get_EL(em_mask)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to extract EL information from the embryo mask %%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EL_info (output): EL information (x0,y0,x1,y1,EL_length);
%% em_mask (input): embryo mask;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
em_mask = bwconvhull(em_mask);
[yall,xall] = find(bwperim(em_mask));
xy_all = [xall,yall];

if size(em_mask,2)/size(em_mask,1) > 1.4
    distXY = pdist2(xy_all,xy_all);
    axis_length = max(max(distXY));
    [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
    x0 = xy_all(extreme0,1);
    y0 = xy_all(extreme0,2);
    x1 = xy_all(extreme1,1);
    y1 = xy_all(extreme1,2);
    L2_extreme = (x1-x0)^2+(y1-y0)^2;
else
    I1 = xy_all(:,1) < size(mask_stack,2)/2;
    I2 = xy_all(:,1) >= size(mask_stack,2)/2;
    if std(xy_all(I1,2)) <= std(xy_all(I2,2))
        [x0,I0] = min(xy_all(:,1));
        y0 = xy_all(I0,2);
        x1 = 2*size(mask_stack,2)-x0;
        y1 = y0;
        axis_length = x1-x0;
        L2_extreme = axis_length^2;
    else
        [x0,I0] = max(xy_all(:,1));
        y0 = xy_all(I0,2);
        x1 = 2*1-x0;
        y1 = y0;
        axis_length = x1-x0;
        L2_extreme = axis_length^2;
    end
end

EL_info = [x0,y0,x1,y1,L2_extreme];
    end
    function [foci_bw,psize0,g_threshold0] = RNA_seg(seg_bw,max_image,signal_channel,WGA_channel,resolution,image_folder,N_cycle)

%% Transcription foci recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolution0 = 0.13;
L_ratio = (resolution/resolution0);
psize = 0.8:0.1:2;
Vpmax = zeros(size(psize));
Ngmax = zeros(size(psize));
Ng_th = 1e5;
max_image0 = double(max_image(:,:,signal_channel))/double(max(max(max_image(:,:,signal_channel))));
fit_radius = 3; %%% Radius of fit window for local gradient calculation (for point recognition threshold selection)
test_threshold = [0:0.002:0.5]; %%% Point recognition threshold scaning range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Point width optimization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for Ip = 1:length(psize)
    H = -fspecial('log',15,psize(Ip)/L_ratio);
    Ng = zeros(1,length(test_threshold)); %%% Initialization of Ng (# of recognized regions)
    Tg = zeros(1,length(test_threshold)); %%% Initialization of Tg (Threshold value)
    fit_th = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of threshold range (for point recognition threshold profile fit)
    fit_g = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of # of recognized regions under certain threshold (for point recognition threshold profile fit)

%%% Application of point recognition filter: %%%===========================
    g = imfilter(max_image0,H,'replicate');
    g(g<0) = 0;
%%% =======================================================================

%%% Profile (# of recognized regions under certain threshold) collection: %
    for Ig = 1:length(test_threshold)
        bw_g = im2bw(g,test_threshold(Ig));
        g_prop = regionprops(bw_g,'Area');
        Ng(Ig) = length(g_prop);
        Tg(Ig) = test_threshold(Ig);
    end
%%% =======================================================================

%%% Gradient of the profile using linear fit: %%%==========================
    for Ig = (1+fit_radius):(length(test_threshold)-fit_radius)
        fit_th(Ig-fit_radius) = test_threshold(Ig);
        temp_fit = polyfit(test_threshold((Ig-fit_radius):(Ig+fit_radius)),Ng((Ig-fit_radius):(Ig+fit_radius)),1);
        fit_g(Ig-fit_radius) = temp_fit(1);
    end
%%% =======================================================================

    %%% Selection of the threshold as the first local maxima of the gradient: %
    g_max = imregionalmin(abs(fit_g));
    g_max(1) = 0;
    II_max = intersect(find(g_max),find(Ng < Ng_th)-fit_radius);
    if ~isempty(II_max)
        Vpmax(Ip) = fit_g(II_max(1));
        Ngmax(Ip) = Ng(II_max(1)+fit_radius);
    else
        Vpmax(Ip) = -inf;
        Ngmax(Ip) = 0;
    end
end
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Selection of point spread width: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find the maxima of point recognition region and reinitialization: %%%==
figure(10)
[AX,H1,H2] = plotyy(psize,Vpmax,psize,Ngmax);
title(['Vmax, Ng vs pixel size (',image_folder,') cycle = ',num2str(N_cycle)],'Interpreter','none')
xlabel('Point size (pixel)')
set(get(AX(1),'Ylabel'),'String','Ng gradient') 
set(get(AX(2),'Ylabel'),'String','# of recognized regions') 

%[~,Ip0] = min(abs(Vpmax));
[~,Ip0] = max(Vpmax);
psize0 = psize(Ip0);
H = -fspecial('log',15,psize(Ip0)/L_ratio);
Ng = zeros(1,length(test_threshold)); %%% Initialization of Ng (# of recognized regions)
Tg = zeros(1,length(test_threshold)); %%% Initialization of Tg (Threshold value)
fit_th = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of threshold range (for point recognition threshold profile fit)
fit_g = zeros(1,length(test_threshold)-2*fit_radius); %%% Initialization of # of recognized regions under certain threshold (for point recognition threshold profile fit)
g=imfilter(max_image0,H,'replicate');
g(g<0) = 0;

for Ig = 1:length(test_threshold)
    bw_g = im2bw(g,test_threshold(Ig));
    g_prop = regionprops(bw_g,'Area');
    Ng(Ig) = length(g_prop);
    Tg(Ig) = test_threshold(Ig);
end
%%% =======================================================================

%%% Gradient of the profile using linear fit: %%%==========================
for Ig = (1+fit_radius):(length(test_threshold)-fit_radius)
    fit_th(Ig-fit_radius) = test_threshold(Ig);
    temp_fit = polyfit(test_threshold((Ig-fit_radius):(Ig+fit_radius)),Ng((Ig-fit_radius):(Ig+fit_radius)),1);
    fit_g(Ig-fit_radius) = temp_fit(1);
end
%%% =======================================================================

%%% Selection of the threshold as the first local maxima of the gradient: %
g_max = imregionalmin(abs(fit_g));
g_max(1) = 0;
II_max = intersect(find(g_max),find(Ng < Ng_th)-fit_radius);
g_threshold0 = fit_th(II_max(1));
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot of the profile with the selected threshold: %%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(Tg,Ng,g_threshold0*[1,1],[min(Ng),max(Ng)],'--')
xlabel('Threshold value (normalized)')
ylabel('# of recognized regions')
title(['Threshold - # profile and the selected threshold value (',image_folder,') ','cycle = ',num2str(N_cycle),', threshold = ',num2str(g_threshold0),', pixel size = ',num2str(psize(Ip0))],'Interpreter','none')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Transcription foci recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foci_bw = im2bw(g,g_threshold0);
foci_bw = imfill(foci_bw,'holes');

%%% Image output: %%%======================================================
foci_bw0 = imdilate(foci_bw,strel('disk',4));
bw_perim_g = bwperim(foci_bw0);
bw_perim_WGA = bwperim(seg_bw);
new_image(:,:,1:2) = max_image(:,:,[signal_channel,WGA_channel]);
new_image(:,:,2) = 0;
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);

figure(4)
try
    imshow(overlay)
catch
    imshow_lxt(overlay)
end
title([image_folder,' (cycle = ',num2str(N_cycle),', white: transcription foci recognition, blue: nucleus recognition)'],'Interpreter','none')
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    function [nucleus_profile,foci_profile,cytoplasmic_profile] = RNA_profile(seg_bw,cyto_bw,foci_bw,max_image,signal_channel,resolution,image_folder,N_cycle,varargin)

%% Nucleus/cytoplasmic RNA profile calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nucleus_profile = zeros(0);
%cytoplasmic_profile = zeros(0);
resolution0 = 0.09;
L_ratio = (resolution/resolution0);
f_radius = max(ceil(6/L_ratio),3);
Lmin = 0;
Lbin = 0.02;
Lmax = 1;
Lwindow = 0.02;
cyto_L = Lmin:Lbin:Lmax;
nucleus_bin = 0:0.01:1;%0.025:0.05:0.975;
average_radius = 0.03;
intensity_Nbin = 50;
cyto_mean = zeros(size(cyto_L));
max_image = double(max_image);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucleus_prop = regionprops(seg_bw,max_image(:,:,signal_channel),'MeanIntensity','Centroid'); %%% Nucleus mean intensity calculation
nucleusI_prop = regionprops(seg_bw,max_image(:,:,signal_channel).*(~imdilate(foci_bw,strel('disk',f_radius))),'MeanIntensity'); %%% Nucleus background total intensity calculation
nucleusA_prop = regionprops(seg_bw,(~imdilate(foci_bw,strel('disk',f_radius))),'MeanIntensity'); %%% Nucleus background area calculation
nucleus_background = [nucleusI_prop.MeanIntensity]./([nucleusA_prop.MeanIntensity]+([nucleusA_prop.MeanIntensity] == 0)); %%% Nucleus mean background calculation

%%% Calculate raw foci intensity: %%%======================================
Rfoci_prop = regionprops(foci_bw,max_image(:,:,signal_channel),'MaxIntensity','MeanIntensity','Area','Centroid'); %%%%% Raw foci intensity calculation
if ~isempty(varargin)
    Sfoci_prop = regionprops(foci_bw,varargin{1},'MaxIntensity'); %%%%% equivalent foci area calculation
    SSfoci = [Sfoci_prop.MaxIntensity];
    if length(varargin) >= 2
        Ifoci0 = varargin{2};
    else
        Ifoci0 = [1,1];
    end
else
    SSfoci = [Rfoci_prop.Area];
    Ifoci0 = [1,1];
end
%%% =======================================================================

%%% Calculate the foci properties of nuclei: %%%===========================
Nfoci_nucleus = zeros(length(nucleus_prop),1); %%%%% Initialization of nucleus foci # profile
Ifoci_nucleus = zeros(length(nucleus_prop),1); %%%%% Initialization of nucleus foci intensity profile
N_prop = regionprops(foci_bw,bwlabel(seg_bw),'MaxIntensity'); %%%%% Link foci to nuclei
background0 = [0,nucleus_background];
foci_intensity2 = round(([Rfoci_prop.MaxIntensity]-background0(([N_prop.MaxIntensity]+1)))./Ifoci0(2)); %%%%% Refined foci maximal intensity
foci_intensity = round(([Rfoci_prop.MeanIntensity]-background0(([N_prop.MaxIntensity]+1))).*SSfoci./Ifoci0(1)); %%%%% Refined foci total intensity

for I_nucleus = 1:length(nucleus_prop)
    Nfoci_nucleus(I_nucleus) = sum([N_prop.MaxIntensity] == I_nucleus); %%%%% Nucleus foci # profile
    Ifoci_nucleus(I_nucleus) = sum(([N_prop.MaxIntensity] == I_nucleus).*foci_intensity); %%%%% Nucleus foci intensity profile
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  Distance calculation and recording for nucleus: %%%===================
nucleus_distance = zeros(0);
nucleus_profile = zeros(0);
if size(nucleus_prop,1)>0
    xy = [nucleus_prop.Centroid];
    nucleus_xy = [xy(1:2:length(xy)-1);xy(2:2:length(xy))]';
    %%%%% Normalized coordinate calculation:
    if size(seg_bw,2)/size(seg_bw,1) > 1.4
        distXY = pdist2(nucleus_xy,nucleus_xy);
        axis_length = max(max(distXY));
        [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
        x0 = nucleus_xy(extreme0,1);
        y0 = nucleus_xy(extreme0,2);
        x1 = nucleus_xy(extreme1,1);
        y1 = nucleus_xy(extreme1,2);
        L2_extreme = (x1-x0)^2+(y1-y0)^2;
    else
        I1 = nucleus_xy(:,1) < size(seg_bw,2)/2;
        I2 = nucleus_xy(:,1) >= size(seg_bw,2)/2;
        if std(nucleus_xy(I1,2)) <= std(nucleus_xy(I2,2))
            [x0,I0] = min(nucleus_xy(:,1));
            y0 = nucleus_xy(I0,2);
            x1 = 2*size(seg_bw,2)-x0;
            y1 = y0;
            axis_length = x1-x0;
            L2_extreme = axis_length^2;
        else
            [x0,I0] = max(nucleus_xy(:,1));
            y0 = nucleus_xy(I0,2);
            x1 = 2*1-x0;
            y1 = y0;
            axis_length = x1-x0;
            L2_extreme = axis_length^2;
        end
    end
    nucleus_distance = dot((nucleus_xy-repmat([x0,y0],size(nucleus_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(nucleus_xy,1),1),2)/L2_extreme;

    %%%%% nucleus data recording:    
    nucleus_profile = [nucleus_distance,[nucleus_prop.MeanIntensity]',Nfoci_nucleus,Ifoci_nucleus,nucleus_background'];
end
%%% =======================================================================

%%%  Distance calculation and recording for foci: %%%======================
foci_distance = zeros(0);
foci_profile = zeros(0);
if size(Rfoci_prop,1)>0
    xy = [Rfoci_prop.Centroid];
    foci_xy = [xy(1:2:length(xy)-1);xy(2:2:length(xy))]';
    %%%%% Normalized coordinate calculation:
    foci_distance = dot((foci_xy-repmat([x0,y0],size(foci_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(foci_xy,1),1),2)/L2_extreme;
    %%%%% nucleus data recording:    
    foci_profile = [foci_distance,[N_prop.MaxIntensity]',foci_intensity',foci_intensity2'];
end
%%% =======================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ylim,xlim] = size(cyto_bw);
xmap = repmat([1:xlim],ylim,1)-x0;
ymap = repmat([1:ylim]',1,xlim)-y0;
Lmap = (xmap*(x1-x0)+ymap*(y1-y0))/L2_extreme;

for Lcenter = 1:length(cyto_L)
    cyto_map = (Lmap >= cyto_L(Lcenter)-Lwindow) & (Lmap <= cyto_L(Lcenter)+Lwindow) & cyto_bw & (~imdilate(foci_bw,strel('disk',f_radius)));
    cyto_mean(Lcenter) = sum(sum(cyto_map.*max_image(:,:,signal_channel)))/sum(sum(cyto_map));
end
cytoplasmic_profile = [cyto_L',cyto_mean'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
    bin_max = min(nucleus_bin+average_radius,1);
    bin_min = max(nucleus_bin-average_radius,0);
    fi0 = zeros(size(nucleus_bin));
    fi1 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fi0(I_bin) = mean(nucleus_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
        fi1(I_bin) = std(nucleus_profile((nucleus_distance >= bin_min(I_bin))&(nucleus_distance <= bin_max(I_bin)),4));
    end
    
    plot(nucleus_profile(:,1),nucleus_profile(:,2),'o',nucleus_profile(:,1),nucleus_profile(:,4),'+',nucleus_profile(:,1),nucleus_profile(:,5),'^',cytoplasmic_profile(:,1),cytoplasmic_profile(:,2),'*',nucleus_bin',fi0,nucleus_bin',fi1)
    title(['FISH profile: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Intensity (A.U.)')
    legend('Nucleus mean intensity','Nucleus foci intensity','Nucleus background mean intensity','Cytoplasmic mean intensity','Mean nucleus foci intensity','Std of nucleus foci intensity')

figure(6)
    bin_max = min(nucleus_bin+average_radius,1);
    bin_min = max(nucleus_bin-average_radius,0);
    fn0 = zeros(size(nucleus_bin));
    fn1 = zeros(size(nucleus_bin));
    fn2 = zeros(size(nucleus_bin));
    fn3 = zeros(size(nucleus_bin));
    for I_bin = 1:length(nucleus_bin)
        fn0(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 0) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 0) <= bin_max(I_bin)));
        fn1(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 1) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 1) <= bin_max(I_bin)));
        fn2(I_bin) = sum((nucleus_distance(Nfoci_nucleus == 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus == 2) <= bin_max(I_bin)));
        fn3(I_bin) = sum((nucleus_distance(Nfoci_nucleus  > 2) >= bin_min(I_bin)).*(nucleus_distance(Nfoci_nucleus  > 2) <= bin_max(I_bin)));
    end
    %fn0 = hist(nucleus_distance(Nfoci_nucleus == 0), nucleus_bin);
    %fn1 = hist(nucleus_distance(Nfoci_nucleus == 1), nucleus_bin);
    %fn2 = hist(nucleus_distance(Nfoci_nucleus == 2), nucleus_bin);
    %fn3 = hist(nucleus_distance(Nfoci_nucleus > 2), nucleus_bin);
    fn_all = fn0+fn1+fn2+fn3+(fn0+fn1+fn2+fn3 == 0);
    plot(nucleus_bin',(1-fn0./fn_all)*100,'m',nucleus_bin',fn1./fn_all*100,'b',nucleus_bin',fn2./fn_all*100,'g',nucleus_bin',fn3./fn_all*100,'r')
    title(['Nucleus foci distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('%')
    legend('nuclei with active foci','nuclei with 1 foci','nuclei with 2 foci','nuclei with >3 foci')

figure(7)
    foci_n = hist(foci_distance,nucleus_bin);
    plot(nucleus_bin',foci_n/(sum(foci_n)+(sum(foci_n) == 0))*100)
    title(['Foci position distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('%')

figure(8)
    [Ifoci_n,xout] = hist(foci_intensity,intensity_Nbin);
    bar(xout,Ifoci_n/(sum(Ifoci_n)+(sum(Ifoci_n) == 0))*100)
    title(['Foci intensity distribution: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('Intensity (A.U.)')
    ylabel('%')

figure(9)
    plot(foci_intensity,foci_intensity2,'.')
    title(['Foci maximal intensity vs total intensity: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel('Maximal intensity (A.U.)')
    ylabel('Total intensity (A.U.)')

clear nucleus_prop nucleusI_prop nucleusA_prop nucleus_background Rfoci_prop

figure(11)
    im0 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus == 0)));
    im1 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus == 1)));
    imfate(:,:,3) = double(im0 | im1);
    clear im1
    im2 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus == 2)));
    imfate(:,:,2) = double(im0 | im2);
    clear im2
    im3 = logical(ismember(bwlabel(seg_bw),find(Nfoci_nucleus > 2)));
    imfate(:,:,1) = double(im0 | im3);
    clear im3
    try
        imshow(imfate);
    catch
        imshow_lxt(imfate);
    end
    title(['Nuclei transcription status: ',image_folder,' (cycle = ',num2str(N_cycle),', white: 0 foci, blue: 1 foci, green: 2 foci, red: > 2 foci)'],'Interpreter','none');
    end
    function [nucleus_profile,cytoplasmic_profile,quanti_p] = protein_profile3(seg_bw,cyto_bw,max_image,signal_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution,varargin)


%% Nucleus/cytoplasmic protein profile calculation: %%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameter setting and initiation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
standard_dorsal = 'Standard protein/Bcd_dorsal.xls';
standard_ventral = 'Standard protein/Bcd_ventral.xls';
subtp = 'background subtraction';
bin2 = 1;
Nb2 = 8;
r_edge = 10;
bcd_dxy = xlsread(standard_dorsal);
bcd_vxy = xlsread(standard_ventral);
%nucleus_profile = zeros(0);
%cytoplasmic_profile = zeros(0);
Lmin = 0;
Lcbin = 0.02;
Lmax = 1;
Lcwindow = 0.02;
cyto_L = Lmin:Lcbin:Lmax;
cyto_mean = zeros(size(cyto_L));

Lmin = 0;
Lnbin = 0.05;
Lmax = 1;
Lnwindow = 0.05;
nu_L = Lmin:Lnbin:Lmax;
nu_mean = zeros(size(nu_L));
nu_std = zeros(size(nu_L));

max_image = double(max_image);
start_limit = 0.3;
end_limit = 1-start_limit;

if isempty(varargin)
    typetxt = 'Protein';
    n1 = 2;
    n2 = 14;
else
    typetxt = varargin{1};
    n1 = varargin{2};
    n2 = varargin{3};
end

sigma = 1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nucleus region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = zeros(0);
Inten = zeros(0);
Inten2 = zeros(0);
mask2 = bwlabel(seg_bw);

for I_layer = 1:size(mask_stack,3)
    sg_prop = regionprops((imerode(mask_stack(:,:,I_layer), strel('disk',r_edge))).*mask2,signal_stack(:,:,I_layer),'MeanIntensity','PixelValues');
    sg2_prop = regionprops((imerode(mask_stack(:,:,I_layer), strel('disk',r_edge))).*mask2,signal_stack(:,:,I_layer).^2,'MeanIntensity');
    if I_layer == 1
        temp0 = cell(size(mask_stack,3),max(max(mask2)));
        temp = zeros(size(mask_stack,3),max(max(mask2)));
        temp2 = zeros(size(mask_stack,3),max(max(mask2)));
    end
    if ~isempty(sg_prop)
        temp0(I_layer,1:length(sg2_prop)) = {sg_prop.PixelValues};  %%%%% Binomial partition preparation
        temp(I_layer,1:length(sg2_prop)) = [sg_prop.MeanIntensity];
        temp2(I_layer,1:length(sg2_prop)) = [sg2_prop.MeanIntensity]-[sg_prop.MeanIntensity].^2;
    end
end
sg_prop = regionprops((imerode(seg_bw, strel('disk',r_edge))).*mask2,'Centroid');
centr = cell2mat({sg_prop.Centroid}');
temp(isnan(temp)) = 0;
[Inten,max_index] = max(temp,[],1);
Inten2 = temp2(max_index+size(temp,1)*[0:(size(temp,2)-1)]);
temp00 = temp0(max_index+size(temp,1)*[0:(size(temp,2)-1)]);

%%%%% Binomial partition analysis: %%%%% ----------------------------------
dF = zeros(length(temp00),1);
FF = zeros(length(temp00),1);
SF = zeros(length(temp00),1);
for I_region = 1:length(temp00)
    Inten_value = temp00{I_region};
    size_region = size(Inten_value,1);
    if floor(size_region/2) < ceil(size_region/2)
        size_region = size_region-1;
    end
    dF(I_region) = ((sum(Inten_value(1:(size_region/2)))-sum(Inten_value((size_region/2+1):size_region)))/2).^2;
    FF(I_region) = sum(Inten_value(1:size_region));
    SF(I_region) = size_region;
end
%%%%% ---------------------------------------------------------------------


%%% Normalized coordinate calculation: %%%=================================
if ~isempty(centr)
    xy = centr;
    nucleus_xy = xy(:,1:2);
    if size(seg_bw,2)/size(seg_bw,1) > 1.4
        distXY = pdist2(nucleus_xy,nucleus_xy);
        axis_length = max(max(distXY));
        [extreme0,extreme1] = find(distXY >= axis_length,1); %%% find the anterior and posterior pole nucleus
        x0 = nucleus_xy(extreme0,1);
        y0 = nucleus_xy(extreme0,2);
        x1 = nucleus_xy(extreme1,1);
        y1 = nucleus_xy(extreme1,2);
        L2_extreme = (x1-x0)^2+(y1-y0)^2;
    else
        I1 = nucleus_xy(:,1) < size(seg_bw,2)/2;
        I2 = nucleus_xy(:,1) >= size(seg_bw,2)/2;
        if std(nucleus_xy(I1,2)) <= std(nucleus_xy(I2,2))
            [x0,I0] = min(nucleus_xy(:,1));
            y0 = nucleus_xy(I0,2);
            x1 = 2*size(seg_bw,2)-x0;
            y1 = y0;
            axis_length = x1-x0;
            L2_extreme = axis_length^2;
        else
            [x0,I0] = max(nucleus_xy(:,1));
            y0 = nucleus_xy(I0,2);
            x1 = 2*1-x0;
            y1 = y0;
            axis_length = x1-x0;
            L2_extreme = axis_length^2;
        end
    end
    nucleus_distance = dot((nucleus_xy-repmat([x0,y0],size(nucleus_xy,1),1)),repmat([(x1-x0),(y1-y0)],size(nucleus_xy,1),1),2)/L2_extreme;
    
    nucleus_profile = [nucleus_distance,Inten'];
    fluc_profile = [nucleus_distance,Inten2'];   %%% intensity variance profile
end
%%% =======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cytoplasmic region recognition/processing: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ylim0,xlim0] = size(cyto_bw);
xmap = repmat([1:xlim0],ylim0,1)-x0;
ymap = repmat([1:ylim0]',1,xlim0)-y0;
clear xlim0 ylim0
Lmap = (xmap*(x1-x0)+ymap*(y1-y0))/L2_extreme;

for Lcenter = 1:length(cyto_L)
    cyto_map = (Lmap >= cyto_L(Lcenter)-Lcwindow) & (Lmap <= cyto_L(Lcenter)+Lcwindow) & cyto_bw;
    cyto_mean(Lcenter) = sum(sum(cyto_map.*max_image(:,:,signal_channel)))/sum(sum(cyto_map));
end
cytoplasmic_profile = [cyto_L',cyto_mean'];

for Lcenter = 1:length(nu_L)
    nu_map = (nucleus_profile(:,1) >= nu_L(Lcenter)-Lnwindow) & (nucleus_profile(:,1) <= nu_L(Lcenter)+Lnwindow);
    nu_mean(Lcenter) = mean(nucleus_profile(nu_map,2));
    nu_std(Lcenter) = std0(nucleus_profile(nu_map,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Standard profile reshaping: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = -3;
nu_x = nucleus_profile(:,1);
standard_x = [bcd_dxy(:,1);bcd_vxy(:,1)];
if mean(nucleus_profile(nucleus_profile(:,1) <= start_limit,2)) < mean(nucleus_profile(nucleus_profile(:,1) >= end_limit,2))
    nu_x = 1-nucleus_profile(:,1);
    bcd_dxy(:,1) = 1-bcd_dxy(:,1);
    bcd_vxy(:,1) = 1-bcd_vxy(:,1);
end
%beta1 = nlinfit(nu_x,nucleus_profile(:,2),@exp_C,[(max(nucleus_profile(:,2))-min(nucleus_profile(:,2))),b2,min(nucleus_profile(:,2))]);
%protein_min = beta1(3);
%protein_max = beta1(1)+beta1(3);
protein_min = min(nu_mean);
protein_max = max(nu_mean);

%beta2 = nlinfit(standard_x,[bcd_dxy(:,2);bcd_vxy(:,2)],@exp_C,[(max([bcd_dxy(:,2);bcd_vxy(:,2)])-min([bcd_dxy(:,2);bcd_vxy(:,2)])),b2,min([bcd_dxy(:,2);bcd_vxy(:,2)])]);
%standard_min = beta2(3);
%standard_max = beta2(1)+beta2(3);
standard_min = min(min(bcd_dxy(:,2)),min(bcd_vxy(:,2)));
standard_max = max(max(bcd_dxy(:,2)),min(bcd_vxy(:,2)));

bcd_dxy(:,2) = (bcd_dxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
bcd_vxy(:,2) = (bcd_vxy(:,2)-standard_min)*(protein_max-protein_min)/(standard_max-standard_min)+protein_min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(n1)
    clf
    if isempty(subtp);
        sub_nu = 0;
        sub_cyto = 0;
    else
        sub_nu = min(nucleus_profile(:,2));
        sub_cyto = min(cytoplasmic_profile(:,2));
    end
    plot(nucleus_profile(:,1),nucleus_profile(:,2)-sub_nu,'o',cytoplasmic_profile(:,1),cytoplasmic_profile(:,2)-sub_cyto,'*',bcd_dxy(:,1),bcd_dxy(:,2)-sub_nu,'X',bcd_vxy(:,1),bcd_vxy(:,2)-sub_nu,'X')
    hold on
    errorbar(nu_L,nu_mean-sub_nu,nu_std,'k');
    hold off
    title([typetxt,' concentration profile: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'],'Interpreter','none')
    xlabel('AP axis (normalized)')
    ylabel('Concentration (A.U.)')
    legend('Nucleus concentration','Cytoplasmic concentration','Standard dorsal profile','Standard ventral profile','Averaged nuclear profile')

figure(n2)
    clf
    nulabel = bwlabel(seg_bw)+1;
    if max(nulabel(:)) > 20
        nuintensity = [-inf;log2(nucleus_profile(:,2)-sub_nu)];
        imagesc(nuintensity(nulabel))%,[min(nucleus_profile(:,2)),max(nucleus_profile(:,2))])
        %axis off
        prolim = max(nuintensity);
        %colormap(jet(1024));
        hcb = colorbar;
        set(hcb,'YTickMode','manual') 
        set(hcb,'YTick',prolim+[-(Nb2-1):1:0]*bin2) 
        nylabel = cell(1,Nb2);
        for ib2 = 1:Nb2
            nylabel{ib2} = ['1/',num2str(2^((Nb2-ib2)*bin2))];
        end
        set(hcb,'YTickLabel',nylabel);
        title(['Nucleus ',typetxt,' concentration plot: ',image_folder,', cycle = ',num2str(N_cycle),' (',subtp,')'],'Interpreter','none')
    end
    

%%% Noise output: %%%======================================================
if size(nucleus_profile,1) >= 2
    p = polyfit(nucleus_profile(:,2),fluc_profile(:,2),1);
else
    p = [0,0];
end
figure(62)
    clf
    plot(nucleus_profile(:,2),fluc_profile(:,2),'b.')
    hold on
    x0 = xlim;
    plot(x0,p(1)*x0+p(2),'r')
    hold off
    xlabel('Mean Intensity');
    ylabel('Intensity varience');
    Cmax = (max(nucleus_profile(:,2))-min(nucleus_profile(:,2)))./p(1)./resolution./resolution./6.02e8./4./pi./sigma.^2;
    title(['Protein Intensity noise plot (',image_folder,'): Var = ',num2str(p(1)),' * I + ',num2str(p(2)),', Cmax = ',num2str(Cmax),' M'],'Interpreter','none')
    legend('Data points','Linear fit')
    quanti_p = [p(1).*resolution.*resolution.*6.02e8.*4.*pi.*sigma.^2,p(2)/p(1)];
%%% =======================================================================

%%% Binomial plot: %%% ====================================================
[FF,IX] = sort(FF);
dF = dF(IX);
bin_size = 20;
if length(FF) >= bin_size
    xbin = zeros(1,floor(length(FF)/bin_size)-2);
    ybin = zeros(1,floor(length(FF)/bin_size)-2);
    sbin = zeros(1,floor(length(FF)/bin_size)-2);
    for Ibin = 0:(floor(length(FF)/bin_size)-3)
        ybin(Ibin+1) = mean(dF((Ibin*bin_size+1):((Ibin+1)*bin_size)));
        xbin(Ibin+1) = mean(FF((Ibin*bin_size+1):((Ibin+1)*bin_size)));
        sbin(Ibin+1) = std0(dF((Ibin*bin_size+1):((Ibin+1)*bin_size)));
    end

    pF = polyfit(xbin,ybin,1);
    %pF(1) = (mean(xbin.*ybin)-mean(xbin).*mean(ybin))/(std(xbin,1)^2);
    %pF(2) = 0;
else
    xbin = FF;
    ybin = dF;
    sbin = zeros(size(FF));
    pF = [0,0];
end

figure(63)
    clf
    plot(FF,dF,'b.')
    hold on
    x0 = xlim;
    plot(x0,pF(1)*x0+pF(2),'r')
    
    errorbar(xbin,ybin,sbin,'g.')
    
    hold off
    xlabel('Total cell intensity (FF)');
    ylabel('Binomial partition variation (dF^2)');
    Cmax = max((FF-pF(2)/pF(1))./SF)./7.26./pF(1)./resolution./resolution./6.02e8./sigma.^2;
    title(['Protein binomial partition noise plot (',image_folder,'): dF^2 = ',num2str(pF(1)),' * FF + ',num2str(pF(2)),', Cmax = ',num2str(Cmax),' M'],'Interpreter','none')
    legend('Binomial partition noise','Linear fit of binomial partition noise')
    

%%% =======================================================================


    nucleus_profile = [nucleus_profile,fluc_profile(:,2)];

    
    function y = std0(x)
        y = std(x)/sqrt(length(x));
    end   
    end
    function dual_profile(nucleus_protein_profile,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Protein - RNA regulation curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_pro = 100;
w_pro = 0.05;
reg_min = 0.25;
reg_max = 0.75;
subtp = 'background subtraction';
single_add = '_new';

if ~isempty(varargin) && size(nucleus_protein_profile,2) >= 4
    z_size = varargin{1};
    reg_I = (nucleus_protein_profile(:,1) >= reg_min) & (nucleus_protein_profile(:,1) <= reg_max) & (nucleus_protein_profile(:,4) > 1) & (nucleus_protein_profile(:,4) < z_size);
else
    reg_I = (nucleus_protein_profile(:,1) >= reg_min) & (nucleus_protein_profile(:,1) <= reg_max);
end
if length(varargin) > 1
    unit1 = varargin{2}{1};
    unit2 = varargin{2}{2};
else
    unit1 = 'A.U.';
    unit2 = 'A.U.';
end
if length(varargin) > 2
    figure_n = varargin{3};
else
    figure_n = [12,15,13];
end


if isempty(subtp)
    pmin = 0;
else
    pmin = min(nucleus_protein_profile(reg_I,2));
end
pmax = max(nucleus_protein_profile(reg_I,2));
pmax = pmax-pmin;
pmin = 0;
pbin = (pmax-pmin)/N_pro;
pwindow = (pmax-pmin)*w_pro;
pro_center = (pmin+pbin/2):pbin:(pmax-pbin/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Background subtraction: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro_center = pro_center-min(nucleus_protein_profile(:,2));
p_profile = nucleus_protein_profile(:,2)-min(nucleus_protein_profile(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regulation curve plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(figure_n(1))
clf
RNA_mean = zeros(size(pro_center));
RNA_std = zeros(size(pro_center));
h1 = zeros(0);
h2 = zeros(0);
    for Ip = 1:length(pro_center)
        RNA_mean(Ip) = mean(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),4));
        RNA_std(Ip) = std(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),4));
    end
    plot(p_profile(reg_I),nucleus_RNA_profile((reg_I),4),'bo')
    hold(gca,'on')
    errorbar(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),RNA_std((RNA_std)>0),'.k')
    if sum((RNA_std)>0) > 5% && sum(nucleus_RNA_profile((reg_I),3)) > 20 && pro_center(end) > 10000 && isempty(strfind(image_folder((end-4):end), single_add))
        %%% Nonlinear fitting 1: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I_max] = max(RNA_mean);
        pro_I = find(pro_center <= pro_center(I_max));
        [~,I_middle] = min(abs(RNA_mean(pro_I)-max(RNA_mean)/2));
        beta0 = [4,pro_center(pro_I(I_middle)),max(RNA_mean)];
        try
            [beta,r,~,~,~] = nlinfit(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),@Hill,beta0);
            plot(pro_center(~isnan(RNA_std)),Hill(beta,pro_center(~isnan(RNA_std))),'r')
            h1 = beta(1);
        catch err
            plot([],[],'r')
            h1 = NaN;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Nonlinear fitting 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            [beta,r,~,~,~] = nlinfit(p_profile(reg_I),nucleus_RNA_profile((reg_I),4),@Hill,beta0);
            plot(pro_center(~isnan(RNA_std)),Hill(beta,pro_center(~isnan(RNA_std))),'--g')
            h2 = beta(1);
        catch err
            plot([],[],'--g')
            h2 = NaN;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        plot([],[],'r')
        plot([],[],'--g')
        h1 = NaN;
        h2 = NaN;
    end
    ax1 = gca;
    ax2=axes('position',get(ax1,'position'),'yaxislocation','right','color','none','YColor','m','XLim',get(ax1,'XLim'),'XTick',get(ax1,'XTick'),'XTickLabel','');
    hold(ax2,'on')
    plot(ax2,pro_center(~isnan(RNA_std)),RNA_std(~isnan(RNA_std))./RNA_mean(~isnan(RNA_std)),'-*m')
    
    hold(ax1,'off')
    hold(ax2,'off')
    set(ax1,'Box','off');
    title(ax1,['Nucleus regulation (RNA level) curve: ',image_folder,', cycle = ',num2str(N_cycle),', h1 = ',num2str(h1),', h2 = ',num2str(h2)],'Interpreter','none')
    xlabel(ax1,['Protein level (',unit1,')'])
    ylabel(ax1,['RNA level (',unit2,')'])
    ylabel(ax2,['Std/mean of RNA level (',unit2,')'])
    h1 = get(ax1,'Children');
    h2 = get(ax2,'Children');
    legend([h1(end:-1:1);h2(end:-1:1)],'Individual nuclei','Averaged profile','Fitting curve for averaged data','Fitting curve for raw data')

figure(figure_n(2))
clf
RNA_mean = zeros(size(pro_center));
RNA_std = zeros(size(pro_center));
    for Ip = 1:length(pro_center)
        RNA_mean(Ip) = mean(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),3));
        RNA_std(Ip) = std(nucleus_RNA_profile(reg_I&(p_profile >= (pro_center(Ip)-pwindow/2))&(p_profile <= (pro_center(Ip)+pwindow/2)),3));
    end
    plot(p_profile(reg_I),nucleus_RNA_profile(reg_I,3),'bo')
    hold(gca,'on')
    errorbar(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),RNA_std((RNA_std)>0),'.k')
    if sum((RNA_std)>0) > 5% && sum(nucleus_RNA_profile((reg_I),3)) > 20 && pro_center(end) > 10000 && isempty(strfind(image_folder((end-4):end), single_add))
        %%% Nonlinear fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I_max] = max(RNA_mean);
        pro_I = find(pro_center <= pro_center(I_max));
        [~,I_middle] = min(abs(RNA_mean(pro_I)-max(RNA_mean)/2));
        beta0 = [4,pro_center(pro_I(I_middle)),max(RNA_mean)];
        try
            [beta,r,~,~,~] = nlinfit(pro_center,RNA_mean,@Hill,beta0);
            plot(pro_center((RNA_std)>0),Hill(beta,pro_center((RNA_std)>0)),'r')
            h1 = beta(1);
        catch err
            plot([],[],'r')
            h1 = NaN;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        plot([],[],'r')
        h1 = NaN;
    end
    ax1 = gca;
    ax2=axes('position',get(ax1,'position'),'yaxislocation','right','color','none','YColor','m','XLim',get(ax1,'XLim'),'XTick',get(ax1,'XTick'),'XTickLabel','');
    hold(ax2,'on')
    plot(ax2,pro_center((RNA_std)>0),RNA_std((RNA_std)>0)./RNA_mean((RNA_std)>0),'-*m')
    
    hold(ax1,'off')
    hold(ax2,'off')
    set(ax1,'Box','off');
    title(['Nucleus regulation (foci #) curve: ',image_folder,', cycle = ',num2str(N_cycle),', h = ',num2str(h1)],'Interpreter','none')
    xlabel(ax1,['Protein level (',unit1,')'])
    ylabel(ax1,'# of foci per nucleus')
    ylabel(ax2,'Std/mean of RNA #')
    h1 = get(ax1,'Children');
    h2 = get(ax2,'Children');
    legend([h1(end:-1:1);h2(end:-1:1)],'Individual nuclei','Averaged profile','Fitting curve','Std/mean of averaged profile')
    
figure(figure_n(3))
clf
RNA_mean = zeros(size(pro_center));
RNA_std = zeros(size(pro_center));
    if ~isempty(foci_RNA_profile)
        I_foci = find(foci_RNA_profile(:,2) > 0);
        I_foci = I_foci((nucleus_protein_profile(foci_RNA_profile(I_foci,2),1) >= reg_min)&(nucleus_protein_profile(foci_RNA_profile(I_foci,2),1) <= reg_max));
        for Ip = 1:length(pro_center)
            RNA_mean(Ip) = mean(foci_RNA_profile(I_foci((p_profile(foci_RNA_profile(I_foci,2)) >= (pro_center(Ip)-pwindow/2))&(p_profile(foci_RNA_profile(I_foci,2)) <= (pro_center(Ip)+pwindow/2))),3));
            RNA_std(Ip) = std(foci_RNA_profile(I_foci((p_profile(foci_RNA_profile(I_foci,2)) >= (pro_center(Ip)-pwindow/2))&(p_profile(foci_RNA_profile(I_foci,2)) <= (pro_center(Ip)+pwindow/2))),3));
        end
        plot(p_profile(foci_RNA_profile(I_foci,2)),foci_RNA_profile(I_foci,3),'o')
        hold on
        errorbar(pro_center((RNA_std)>0),RNA_mean((RNA_std)>0),RNA_std((RNA_std)>0),'k')
        hold off
        legend('Individual nuclei','Averaged profile')
    else
        I_foci = [];
        plot([],[],'o')
    end
    title(['Foci regulation curve: ',image_folder,', cycle = ',num2str(N_cycle)],'Interpreter','none')
    xlabel(['Protein level (',unit1,')'])
    ylabel(['RNA level (',unit2,')'])
    
    function y = std0(x)
    NN = size(x);
    if NN(1) == 1;
        N1 = NN(2);
    else
        N1 = NN(1);
    end
    y = std(x)/sqrt(N1);
    end           
    function y = Hill(beta,x)
     y = beta(3)*x.^beta(1)./(x.^beta(1)+beta(2)^beta(1));
    end
        
    end 
    function bwout = reseg(immask2,area1)    



    se_mask = regionprops(immask2,'Area');
    se_area = [se_mask.Area];
    bwout = false(size(immask2));   %%% Output mask initialization

    r1 = sqrt(area1/pi);   %%% estimated radius of a single nucleus

%% Partition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bwraw = immask2;
%bwraw = bwareaopen(bwraw,round(area1*0.5));   %%% remove small areas

%%% Estimation of the number of nuclei included in each recognized regions
se_mask = regionprops(bwraw,'Area');   %%% area calculation
nmask = round([se_mask.Area]./area1);   %%% calculate the estimated # of nuclei for each recognized area
lmask = bwlabel(bwraw);   %%% mask labeling

%%% Filter out the merged nuclei
bwout(ismember(lmask,find(nmask <= 1))) = true;   %%% Add the single nuclei areas to the output mask
bwraw(ismember(lmask,find(nmask <= 1))) = false;   %%% Remove the single nuclei areas from the raw mask

if any(any(bwraw))
    %%% Watersheding to chop merged nuclei
    D= bwdist(~bwraw);
    g = imfilter(double(bwraw),double(getnhood(strel('disk',round(r1*0.75)))),'symmetric','corr');
    g2 = imimposemin(-D,imdilate(imregionalmax(g),strel('disk',5)));
    L2 = watershed(g2);
    bwnew = bwraw & (~(L2 == 0));
    bwout = bwout | bwnew;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convex the recognized regions: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% se_mask = regionprops(bwout,'Area');
% se_area = [se_mask.Area];
% bwout = ismember(bwlabel(bwout),find(se_area > 1000));

%bwout0 = bwout;
se_mask = regionprops(bwout,'ConvexImage','BoundingBox');
bwout = zeros(size(bwout));
bwp = zeros(size(bwout));
if size(se_mask(:),1)>0
    for I_temp = 1:length(se_mask)
        x1 = uint16(se_mask(I_temp).BoundingBox(2));
        x2 = uint16(x1+se_mask(I_temp).BoundingBox(4)-1);
        y1 = uint16(se_mask(I_temp).BoundingBox(1));
        y2 = uint16(y1+se_mask(I_temp).BoundingBox(3)-1);
        bwout(x1:x2,y1:y2) = bwout(x1:x2,y1:y2) + double(se_mask(I_temp).ConvexImage);
        bwp(x1:x2,y1:y2) = bwp(x1:x2,y1:y2) | bwperim(se_mask(I_temp).ConvexImage);
    end
end

bwout = logical(bwout == 1) & (~bwp);
bwout = imerode(bwout,strel('disk',5));
bwout = bwconvhull(bwout,'objects');
bwout = bwmorph(bwout,'thicken',5);
end



%% ===== Step3. Nuclear segmentation =====
function DAPI_seg3D2(p)
%clear all
close all
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lf_name = 'Confocal/Duallist.mat';
in_folder = 'stacks/';
out_folder = 'masks/';
input_name = 'matchlist.mat';
mask_name = 'mask.mat';
standard_record = 'Calibration/Results/standard.mat';
image_type = '*.tif';
%out_folder = 'Results/';

% imblur0 = 10;   %%% Initial bluring mask radius
imblur0 = 3;   %%% Initial bluring mask radius
blur_fit = [0.0033,-7.5];   %%% Bluring mask radius recalculation parameters
bad_ratio = 1.5;   %%% Threshold for being a merged (bad) nuclei center mask region
I_max = 256;   %%% threshold scanning step for circular mask algorithm
% % low_th = 1500;    %%% lower area limit for circular mask algorithm
low_th = 1000;    %%% lower area limit for circular mask algorithm
% low_th = 100;    %%% lower area limit for circular mask algorithm
high_th = 50000;   %%% higher area limit for circular mask algorithm
cir_th = 0.75;   %%% circularity threshold for circular mask algorithm
% conv_th = 0.9;   %%% convexity threshold for circular mask algorithm
sigma_th = 100;  %%% neighboring mask threshold correlation decay length for supplemental mask recognition
merge_th = 100000;   %%% higher area limit (nuclei merge limit) for supplemental mask recognition
th_ratio = 1.0;   %%% threshold resetting ratio for mask refinement on z direction
back_th = 1.5;
% back_th = 0.6;
r_refine = 10;
% r_refine = 5;
if ispc==1
    lf_name(findstr(lf_name, '/'))='\';
    in_folder(findstr(in_folder, '/'))='\';
    out_folder(findstr(out_folder, '/'))='\';
    standard_record(findstr(standard_record, '/'))='\';
end

% Advanced Guassian filtering    
Ex = fspecial('gaussian',100,20);
Ix = fspecial('gaussian',100,30);
% Ex = fspecial('gaussian',100,10);
% Ix = fspecial('gaussian',100,20);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global standard_data
% standard_data = load(standard_record);

if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.mat')
    list_name = lf_name;
    %[num_list, folder_list] = xlsread(list_name);
    load(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
    % ispc=1 then use'\'
    if ispc==1
        folder_list{1}(findstr(folder_list{1}, '/'))='\';
    end

elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);

for list_I = 1:N1
    %[sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    load([folder_list{list_I,1},in_folder,input_name]);
    sub_num = match_num;
    sub_list = match_folder;
    
    [M1,M2] = size(sub_list);
     channel_name = folder_list{list_I,5};
    for list_J = 1:M1
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        Mdim = sub_num(list_J,3);
        Nbin1 = sub_num(list_J,2);
        WGA_channel = sub_num(list_J,6);
        DAPI_channel = sub_num(list_J,7);
        RNA_channel = sub_num(list_J,10);
        protein_channel = sub_num(list_J,8);
        resolution = sub_num(list_J,9);
        imblur = imblur0;
        
        imlist = dir([folder_list{list_I,1},in_folder,sub_list{list_J,3},image_type]); %%% get the image list from image folder1
%         max_temp = zeros(1,length(imlist));
% matlabpool
%         parfor image_I = 1:length(imlist)
%             raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
%             max_temp(image_I) = max(max(raw_im(:,:,DAPI_channel)));
%         end
% matlabpool close
%         clear raw_im
%         im_max = max(max_temp);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% segmentation test and mean nuclear size measurement: %%%%%%%%%%%%%%%%%%%
        raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(round(length(imlist)/2)).name]);
        input_im = raw_im(:,:,DAPI_channel);
        size_im = size(input_im);
        H = fspecial('disk',imblur); % Filter Kernel  
        fil_im = imfilter(input_im,H,'same','conv');
        im_max = max(fil_im(:));
        im_temp = fil_im;
        bw_im = false(size(fil_im));

        for I = 1:(I_max-1)
            th_I = I/I_max*double(im_max)/65535;
            bw_temp = imfill(im2bw(im_temp,th_I),'holes');
            bw_prop = regionprops(bw_temp,'Area','Perimeter');%,'ConvexArea');
            bw_area = [bw_prop.Area];
%             bw_conv_area = [bw_prop.ConvexArea];
            bw_perim = [bw_prop.Perimeter];
            ind_true = find(bw_area >= low_th & bw_area <= high_th & 4*pi*bw_area./bw_perim.^2 >= cir_th);% & bw_area./bw_conv_area >= conv_th);
            bw_true = ismember(bwlabel(bw_temp),ind_true);
            bw_im = bw_im | bw_true;
            im_temp(bw_true) = 0;
        end
        bw_im = imopen(bw_im,strel('disk',10));
        temp_prop = regionprops(bw_im,'Area');
        mean_area = mean([temp_prop.Area]);
        imblur = max(3,round(mean_area*blur_fit(1)+blur_fit(2)));
%         r_mean = max(ceil(sqrt(mean_area/pi)),28);
% 
%         Ex = fspecial('gaussian',100,round(r_mean*0.4));
%         Ix = fspecial('gaussian',100,round(r_mean*0.6));
        
        clear raw_im input_im fil_im im_temp bw_im bw_temp bw_prop bw_area bw_perim bw_true temp_prop
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


        raw3D = false([size_im([1:2]),length(imlist)]); %%% Raw 3D mask
        th_3D = zeros([size_im([1:2]),length(imlist)]); %%% Raw 3D threshold
        fil3D = zeros([size_im([1:2]),length(imlist)]); %%% Smoothed image
        
parpool
        parfor image_I = 1:length(imlist)
            raw_im = imread([folder_list{list_I,1},in_folder,sub_list{list_J,3},imlist(image_I).name]);
            input_im = raw_im(:,:,DAPI_channel);
        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Targeting nuclear centers: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            H = fspecial('disk',imblur); % Filter Kernel       
            I_blur = imfilter(input_im,H,0); %Apply Filter
            I_blur =  adapthisteq(I_blur); %Step 1: enhances the contrast of the grayscale image   

            %%% Gaussian difference filtering and segmentation  
            outE = imfilter(single(I_blur),Ex,'replicate'); 
            outI = imfilter(single(I_blur),Ix,'replicate'); 
            outims = outE - outI;  
            W = watershed(max(outims(:))-outims);

            outims(outims < 0) = 0;
            outims = outims/max(outims(:));

            bw_diff = im2bw(outims,1.2*graythresh(outims)); 
            temp2_props = regionprops(bw_diff,'Area');
            mean_area_diff = mean([temp2_props.Area]);
            bad_nu_I = find([temp2_props.Area] > bad_ratio*mean_area_diff);
            bad_bw0 = ismember(bwlabel(bw_diff),bad_nu_I);
            bad_bw = imerode(bad_bw0 & (W ~= 0),strel('disk',5));
            bw_diff(bad_bw0) = bad_bw(bad_bw0);

            %%% Recognize the center points of nuclei
            diff_props = regionprops(bw_diff,'Centroid');
            center_xy = round(cell2mat({diff_props.Centroid}'));
            bw_center = sub2ind(size(bw_diff),center_xy(:,2),center_xy(:,1));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Circular segmentation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            H = fspecial('disk',imblur); % Filter Kernel  
            fil_im = imfilter(input_im,H,'same','conv');
            th_2D = zeros(size_im([1:2])); %%% Raw 2D threshold
            im_max = max(fil_im(:));
            im_temp = fil_im;
            bw_im = false(size(fil_im));
%%% Cicle search
            for I = 1:(I_max-1)
                th_I = I/I_max*double(im_max)/65535;
                bw_temp = imfill(im2bw(im_temp,th_I),'holes');
                bw_prop = regionprops(bw_temp,'Area','Perimeter');
                bw_area = [bw_prop.Area];
                bw_perim = [bw_prop.Perimeter];
                ind_true = find(bw_area >= low_th & bw_area <= high_th & 4*pi*bw_area./bw_perim.^2 >= cir_th);
                bw_true = ismember(bwlabel(bw_temp),ind_true);
                bw_im = bw_im | bw_true;
                im_temp(bw_true) = 0;
                th_2D(bw_true) = th_I;
            end
%             bw_im = imopen(bw_im,strel('disk',10));
            bw_out = bw_im;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Supplement segmentation based on center recognition: %%%%%%%%%%%%%%%%%%%
            bw_props = regionprops(bw_im,'Centroid');
            if ~isempty(bw_props)
                center_xy2 = round(cell2mat({bw_props.Centroid}'));
                bw_center2 = sub2ind(size(bw_im),center_xy2(:,2),center_xy2(:,1));
                th_center2 = th_2D(bw_center2);

                bw_center0 = bw_center(~bw_im(bw_center));
                bw_center0_area = ismember(bwlabel(bw_diff),find((~bw_im(bw_center))));
                [center_y0,center_x0] = ind2sub(size(bw_diff),bw_center0);
                center_xy0 = [center_x0,center_y0];

                for I_center = 1:size(bw_center0,1)
                    center_dist2 = (center_xy2(:,1)-center_xy0(I_center,1)).^2+(center_xy2(:,2)-center_xy0(I_center,2)).^2;
                    th_estimate = sum(th_center2.*exp(-center_dist2/2/sigma_th^2))/sum(exp(-center_dist2/2/sigma_th^2));
                    bw_th = imfill(bwareaclose(im2bw(im_temp,th_estimate),merge_th),'holes');
%                     bw_th = im2bw(im_temp,th_estimate);
                    bw_th_label = bwlabel(bw_th);
                    bw_I = (bw_th_label == bw_th_label(bw_center0(I_center))) & bw_th;
                    bw_out = bw_out | bw_I;
                    th_2D(bw_I) = max(th_2D(bw_I),th_estimate);
                end

                bw_cut = false(size(fil_im));
    %             bw_cut([bw_center2;bw_center0]) = true;
                bw_cut(bw_center2) = true;
                bw_cut = bw_cut | (bw_center0_area & bw_out);

                temp_label = bwlabel(bw_out);
                temp_sel = ismember(temp_label,temp_label(bw_cut));
                D= bwdist(~bw_out);
                g2 = imimposemin(-D,bw_cut);
                bw_cut = imdilate(watershed(g2) == 0,strel('disk',1)) & (temp_sel);
                bw_out = bw_out & (~ bw_cut);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Mask area refinement and convexation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                im_refine = imfilter(input_im,fspecial('disk',5),'same','conv');
                Inten_back = sort(im_refine(~bw_out));
                back_thresh = back_th*Inten_back(round(length(Inten_back)*0.9));
                bw_refine0 = bwareaopen(bw_out & (im_refine >= back_thresh),low_th);
                bw_refine0 = imfill(bw_refine0,'holes');
                bw_refine0 = imopen(bw_refine0,strel('disk',r_refine));
                bw_refine = bwconvhull(bw_refine0,'objects');

                D= bwdist(~bw_refine);
                g2 = imimposemin(-D,bw_refine0);
                bw_separate = imdilate(watershed(g2) == 0,strel('disk',1));
                bw_out = bw_refine & (~ bw_separate);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Raw segmentation output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            raw3D(:,:,image_I) = bw_out;
            th_3D(:,:,image_I) = th_2D;
            fil3D(:,:,image_I) = fil_im;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
% % % %         parfor image_I = 1:size(raw3D,3)
% % % %             raw3D(:,:,image_I) = imopen(raw3D(:,:,image_I),strel('disk',10));
% % % % %             raw3D(:,:,image_I) = imerode(raw3D(:,:,image_I),strel('disk',2));
% % % %         end
delete(gcp) 
        
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D mask tiling and refinement: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         mask_label = label3D(raw3D);   %%% mask tiling
        mask_stack = label3D(raw3D);   %%% mask tiling
        clear raw3D
    %%% mask rethresholding to refine z distribution: %%% =================
%         mask_stack = zeros(size(raw3D));
%         parfor I_nuclei = 1:max(mask_label(:))
%             temp_label = mask_label == I_nuclei;
%             mask_stack(temp_label) = I_nuclei*im2bw(fil3D(temp_label),max(th_3D(temp_label)));
%         end
%         label_prop = regionprops(mask_label,th_3D,'MaxIntensity');
%         th_max = [1,[label_prop.MaxIntensity]];+
%         mask_stack = mask_label.*(fil3D > th_ratio*65535*th_max(mask_label+1));
%         mask_stack = mask_label;
%         clear mask_label th_max
    %%% ===================================================================
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resign the label #: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        new_label = zeros(1,max(mask_stack(:))+1);
        new_label(sort(unique(mask_stack))+1) = [0:(length(unique(mask_stack))-1)];
        mask_stack = new_label(mask_stack+1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation result output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist([folder_list{list_I,1},out_folder,sub_list{list_J,3}],'dir')
            mkdir([folder_list{list_I,1},out_folder,sub_list{list_J,3}])
        end
        save([folder_list{list_I,1},out_folder,sub_list{list_J,3},mask_name],'mask_stack','-v7.3')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end

toc
end


%% ===== Step4. smFISH spot analysis =====
function stack_RNA_GUIfun(p)

%  ==================== Function 1. load_file_Callback  ==============================


%% Setup folder and load images: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout imname imfolder bk_label bk_use mask_out stack_th mask_stack
imstack = zeros(0);
imstack0 = zeros(0);
immask = false(0);
max_spot = zeros(0);
n0 = zeros(0);
xout0 = zeros(0);
n = zeros(0);
xout = zeros(0);
bk_use = cell(0);
%Idim = 2;

% -------- parameters from stack_RNA_OpeningFcn
dlayer = 1;
dcontrast = 1;
dcontrmin = 0;
dfsize = 1.5;
dimthresh = 0.01;
dmask = true;
dpeak = true;
stack_th0 = 0;
stack_th = stack_th0;
bk_label = 'OreR';
mstatus = {'Mask on','Mask off'};
resolution0 = 0.13;
limsize = [0,10];
limthresh = [0,1];
limcontrast = [0,1];
channel2_add = '_RNA2';
standard_record = 'Calibration/Results/standard.mat';
image_tail = '*.tif';
% --------------

% load duallist.mat and matchlist.mat files
dual_name = 'Confocal/Duallist.mat';
match_name = 'matchlist.mat';
mask_sub = 'masks/';
Results_sub = 'Results/';
stack_RNAsave_s = 'stack_RNA_single.mat';
stack_RNAsave_f = 'stack_RNA_foci.mat';
recent_folder=pwd;
F=findstr(recent_folder,'FISHIF');
if size(recent_folder,2)>F+5
    root_folder = [recent_folder(1:F+5) '/'];  
elseif size(recent_folder,2)==F+5
    root_folder = [recent_folder '/'];
else
    display('error folder')
end
if ispc==1
    standard_record(findstr(standard_record, '/'))='\';
    root_folder(findstr(root_folder, '/'))='\';
    dual_name(findstr(dual_name, '/'))='\';
    mask_sub(findstr(mask_sub, '/'))='\';
    Results_sub(findstr(Results_sub, '/'))='\';
end
load([root_folder dual_name]);
standard_record=[root_folder standard_record];
[N1,~] = size(folder_list);

for list_I = 1:N1
    % ========================================
    stacks_folder = [root_folder folder_list{list_I,1} 'stacks/'];
    if ispc==1
        folder_list{list_I,1}(findstr(folder_list{list_I,1}, '/'))='\';
        stacks_folder(findstr(stacks_folder, '/'))='\';
    end
    load([stacks_folder,match_name]);
    [M1,~] = size(match_folder);
    num_list=match_num;
    
    for list_J = 1:M1
        imfolder = stacks_folder;
        if size(match_folder,2)>2
            imname =match_folder{list_J,3};
        else
            imname =match_folder{list_J,1};
        end
        %if imname(end)=='/'||imname(end)=='\'
        %    imname=imname(1:end-1);
        %end
        if ispc==1
            imname(findstr(imname, '/'))='\';
        end
        standard_data = load(standard_record);


        %% Image stack loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if isempty(varargin)
        %    imfolder = uigetdir('','Select a stack to load');   %%% Pick up a image stack folder
            %imfolder = stacks_folder;
            %imname =match_folder{1,1}; 
            %imname = imfolder(find((imfolder == '/') | (imfolder == '\'),1,'last')+1:end);
            %imfolder = imfolder(1:find((imfolder == '/') | (imfolder == '\'),1,'last'));
        %end
        %[num_list, file_list] = xlsread([imfolder,list_name]);
        %[num_list, file_list] = xlsread('E:\not to modify\George vision\original code for mRNA quantification\FISHIF\Confocal\Exp1\stacks\matchlist.xls'); 
        %I_file = strcmp([imname,'/'],file_list(:,1));
        %I_file = I_file | strcmp([imname,'\'],file_list(:,1));
        I_file=list_J;

        lsm_stack = dir([imfolder,imname,image_tail]); %%% image file list loading
        resolution = num_list(list_J,9);
        if p.type_foci ==1 || p.type_single ==1
            RNA_channel = num_list(I_file,10); %%% RNA channel
            % RNA_channel = num_list(I_file,8);
        else
            RNA_channel = num_list(I_file,12); %%% Signal2 channel
        end
        L_ratio = 1; %(resolution/resolution0);

        immax = length(lsm_stack);
        tempmax = zeros(1,immax);

        for I_layer = 1:immax
            temp0 = imread([imfolder,imname,lsm_stack(I_layer).name]);
            temp = temp0(:,:,RNA_channel);
            imstack0(:,:,I_layer) = temp;
            imstack(:,:,I_layer) = imfilter(temp,fspecial('gaussian',3,1),'symmetric','conv');
            tempmax(I_layer) = max(max(imstack(:,:,I_layer)));
        end

        %imstack0 = imstack0(ceil(size(imstack0,1)*1/4):ceil(size(imstack0,1)*3/4),ceil(size(imstack0,2)*1/4):ceil(size(imstack0,2)*3/4),:);
        %imstack  = imstack(ceil(size(imstack,1)*1/4):ceil(size(imstack,1)*3/4),ceil(size(imstack,2)*1/4):ceil(size(imstack,2)*3/4),:);

        if p.type_single == 1 || p.type_single2 == 1
            try
                %load([imfolder(1:find((imfolder(1:end-1) == '/') | (imfolder(1:end-1) == '\'),1,'last')),'masks/',imname,'/mask.mat']);
                load([root_folder folder_list{list_I,1} mask_sub imname 'mask.mat']);
            catch
                mask_stack = true(size(imstack));
            end
        %     mask_stack = true(size(imstack));
            mask1D = max(max(logical(mask_stack),[],3),[],1);
            ELmin = find(mask1D,1);
            ELmax = find(mask1D,1,'last');

            rxmin1 = 2/16;
            rxmax1 = 6/16;
            rymin1 = 2/8;
            rymax1 = 6/8;

            xmin1 = round(ELmin+(ELmax-ELmin)*rxmin1);
            xmax1 = round(ELmin+(ELmax-ELmin)*rxmax1);
            ymin1 = round(rymin1*size(imstack,1));
            ymax1 = round(rymax1*size(imstack,1));


            rxmin2 = 1-rxmax1;%6/16;
            rxmax2 = 1-rxmin1;%2/16;
            rymin2 = rymin1;%2/8;
            rymax2 = rymax1;%6/8;

            xmin2 = round(ELmin+(ELmax-ELmin)*rxmin2);
            xmax2 = round(ELmin+(ELmax-ELmin)*rxmax2);
            ymin2 = round(rymin2*size(imstack,1));
            ymax2 = round(rymax2*size(imstack,1));

            if mean(mean(mean(imstack0(ymin1:ymax1,xmin1:xmax1,:)))) >= mean(mean(mean(imstack0(ymin2:ymax2,xmin2:xmax2,:))))
                imstack0 = imstack0(ymin1:ymax1,xmin1:xmax1,:);
                imstack  = imstack(ymin1:ymax1,xmin1:xmax1,:);
                mask_stack  = mask_stack(ymin1:ymax1,xmin1:xmax1,:);
            else
                imstack0 = imstack0(ymin2:ymax2,xmin2:xmax2,:);
                imstack  = imstack(ymin2:ymax2,xmin2:xmax2,:);
                mask_stack  = mask_stack(ymin2:ymax2,xmin2:xmax2,:);
            end
        else
            mask_stack = true(size(imstack));
        end
        % mask_stack = true(size(imstack));


        temp0 = fspecial('gaussian',5,1);
        kernel0 = zeros(1,1,5);
        kernel0(1,1,:) = temp0(3,:)./sum(temp0(3,:));
        imstack = imfilter(imstack,kernel0,'symmetric','conv');

        imstack = double(imstack)/65535;   %%% Normalization
        imstack0 = double(imstack0);
        immask = false(size(immask));
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%

        %% initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %vlayer = dlayer;
        [~,vlayer] = max(tempmax);
        vcontrast = dcontrast;
        vcontrmin = dcontrmin;
        vfsize = dfsize;
        vimthresh = dimthresh;
        vmask = dmask;
        vpeak = dpeak;
        %{
        set(handles.file_name,'String',[imfolder,imname]);
        set(handles.layer_slide,'Max',immax);
        set(handles.layer_slide,'Min',1);
        set(handles.layer_slide,'SliderStep',[1/(immax-1),1/(immax-1)]);
        set(handles.layer_slide,'Value',vlayer);
        set(handles.layer_value,'String',num2str(vlayer));
        set(handles.max_layer,'String',num2str(immax));
        set(handles.contrast_min,'String',num2str(vcontrmin));
        set(handles.contrast_value,'String',num2str(vcontrast));
        set(handles.fsize,'String',num2str(vfsize));
        set(handles.imthresh,'String',num2str(vimthresh));
        set(handles.Mask_status,'Value',vmask);
        set(handles.Mask_status,'String',mstatus{vmask+1});
        set(handles.Peak_on,'Value',vpeak);
        %}
        % if ~isempty(strfind(imname,bk_label))
        %     set(handles.background_on,'Value',true);
        % else
        %     set(handles.background_on,'Value',false);
        % end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Image output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mask_out = spmask(imstack,1,stack_th);   %%% 3D local maximal mask generation (Primary filter)
        %mask_out = mask_out & (imstack*65535 >= 7000);   %%% Prefilter using peak intensity threshold
        [overlay,g] = cf_show(imstack,vlayer,vfsize,L_ratio,vimthresh,vpeak,vcontrast,vcontrmin,vmask,mask_out);

        figure;imshow(overlay);
        % =======================================================
        
        hist_on_GUIfun(p,imname,imfolder,imstack,imstack0,mask_stack,bk_use,mask_out);
        
        

        % =======================================================
        % save results
        if p.type_single == 1 || p.type_single2 == 1
            save([root_folder folder_list{list_I,1} Results_sub imname stack_RNAsave_s],'mask_out','overlay','imname','imfolder','imstack','imstack0','mask_stack','bk_use','g');
        else
            save([root_folder folder_list{list_I,1} Results_sub imname stack_RNAsave_f],'mask_out','overlay','imname','imfolder','imstack','imstack0','mask_stack','bk_use','g');
        end
    end
end


end
function hist_save_GUIfun(p)
% hObject    handle to hist_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global imstack0 imstack immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast imname imfolder n0 xout0 n xout bk_name sg_name bk_use sub_folder max_spot m_data stack_th lim t0 dt0
hist_tail = '.fig';
tif_tail = '.tif';
xls_tail = '.mat';
hist_name = [imname,hist_tail];
tif_name = [imname,tif_tail];
xls_name0 = [imname,'_all',xls_tail];
xls_name1 = [imname,'_low',xls_tail];
xls_name2 = [imname,'_raw',xls_tail];
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));

%% Save images and data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([imfolder0,sub_folder],'dir')
    mkdir([imfolder0,sub_folder]);
end
hgsave([imfolder0,sub_folder,hist_name]);
saveas(gcf,[imfolder0,sub_folder,tif_name]);
if ~isempty(n0)
    xlswrite([imfolder0,sub_folder,xls_name0],[n0',xout0']);
    xlswrite([imfolder0,sub_folder,xls_name1],[n',xout']);
    xlswrite([imfolder0,sub_folder,xls_name2],max_spot);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Save the record: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if get(handles.background_on,'Value')
        %%% Save background record: %%% =======================================
        if exist([imfolder0,bk_name],'file')
            [bk_data,bk_record,~] = xlsread([imfolder0,bk_name]);
        else
            bk_record = {};
            bk_data = [];
        end
        if ~any(strcmp(xls_name1,bk_record))
            I_bk = size(bk_record,1)+1;
        else
            I_bk = find(strcmp(xls_name1,bk_record),1);
        end
        bk_record{I_bk,1} = xls_name1;
        bk_data(I_bk,:) = m_data;
        
        if size(bk_data,1) < size(bk_record,1)
            bk_data([(size(bk_data,1)+1):size(bk_record,1)],:) = zeros(size(bk_record,1)-size(bk_data,1),size(bk_data,2));
        end
        xlswrite([imfolder0,bk_name],cat(2,bk_record,num2cell(bk_data)));
        %%% ===================================================================
    else
        %%% Save data record: %%% =============================================
        if exist([imfolder0,sg_name],'file')
            [sg_data,sg_record,~] = xlsread([imfolder0,sg_name]);
        else
            sg_record = {};
            sg_data = [];
        end
        if ~isempty(sg_record)
            I_record = find(strcmp(imname,sg_record(:,1)),1);
        else
            I_record = [];
        end
        if isempty(I_record)
            I_record = size(sg_record,1)+1;
        end
        sg_record{I_record,1} = imname;
        sg_data(I_record,:) = m_data;
        if ~isempty(bk_use)
            sg_record(I_record,2:(length(bk_use)+1)) = reshape(bk_use,1,length(bk_use));
        end
        if (length(bk_use)+1) < size(sg_record,2)
            sg_record(I_record,(length(bk_use)+2):size(sg_record,2)) = cell(1,size(sg_record,2)-(length(bk_use)+1));
        end
        
        if size(sg_data,1) < size(sg_record,1)
            sg_data([(size(sg_data,1)+1):size(sg_record,1)],:) = zeros(size(sg_record,1)-size(sg_data,1),size(sg_data,2));
        end
        xlswrite([imfolder0,sg_name],cat(2,sg_record,num2cell(sg_data)));
        %%% ===================================================================
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%% Foci recognition output: %%% ==========================================
foci_add = '_foci';
seg_add = '_seg';
figure_tail = '.fig';
out_folder = 'Results/';
result_folder = [imfolder(1:find((imfolder(1:(end-1)) == '/')|(imfolder(1:(end-1)) == '\'),1,'last')),out_folder];
foci_bw = false(size(imstack0,1),size(imstack0,2));

for ii = 1:size(max_spot,1)
    x_spot = min(max(round(max_spot(ii,6)),1),size(imstack,1));
    y_spot = min(max(round(max_spot(ii,7)),1),size(imstack,2));
    foci_bw(x_spot,y_spot) = true;
end

foci_bw0 = imdilate(foci_bw,strel('disk',4));
bw_perim_g = bwperim(foci_bw0);
%bw_perim_WGA = bwperim(seg_bw);
new_image(:,:,1) = max(imstack0,[],3)./max(max(max(imstack0,[],3)));
new_image(:,:,2) = 0;
new_image(:,:,3) = 0;
overlay = imoverlay(new_image, bw_perim_g, [1,1,1]);
%overlay = imoverlay(overlay,bw_perim_WGA,[0,0,1]);

figure(4)
imshow(overlay)
title([result_folder,imname,', white: transcription foci recognition'],'Interpreter','none')
saveas(4,[imfolder0,sub_folder,imname,foci_add,seg_add,figure_tail]);
close(4)
%%% =======================================================================


%%% Foci recognition mask selection: %%% ==================================
th_add = '_th';
figure(3)
[AX,H1,H2] = plotyy(lim,t0,lim,dt0);
hold on
plot(stack_th*[1,1],ylim,'--')
title(['Threshold - # profile and the selected threshold value (',imfolder0,') threshold = ',num2str(stack_th)],'Interpreter','none')
xlabel('Threshold value (A.U.)')
set(get(AX(1),'Ylabel'),'String','#') 
set(get(AX(2),'Ylabel'),'String','#/#') 
%legend('# of recognized regions','#(n)/#(n-1)','Threshold value')
saveas(3,[imfolder0,sub_folder,imname,foci_add,th_add,figure_tail]);
close(3)
%%% =======================================================================


% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function hist_fit(lf_name)
%clear all
close all

% sub_folder = 'Histogram_A/';
sub_folder = 'Histogram_A_gal4/';
r_th = 0;
sub_folder_nullo = 'Histogram_default/';
in_folder = 'stacks/';
input_name = 'matchlist.mat';
data_tail = '_raw.mat';
output_add = '_spot_fit';
mat_tail = '.mat';
fig_tail = '.fig';
if ispc==1
    sub_folder(findstr(sub_folder, '/'))='\';
    sub_folder_nullo(findstr(sub_folder_nullo, '/'))='\';
    in_folder(findstr(in_folder, '/'))='\'; 
end


hist_min = 0; 
hist_max = 4e5;
hist_bin = 2e3;
fit_initial = [0.035,0.001,0.0001,0.0001,1e4,1e4,1];
fit_lower = [0,0,0,0,0,0,0];

%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(lf_name,'char') && strcmp(lf_name(end-3:end),'.mat')
    list_name = lf_name;
    %[num_list, folder_list] = xlsread(list_name);
    load(list_name);
    folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
    % ispc=1 then use'\'
    if ispc==1
        folder_list{1}(findstr(folder_list{1}, '/'))='\';
    end
    
    
elseif isa(lf_name,'cell')
    folder_list = lf_name;
else
    error('Incorrect input!')
end
[N1,N2] = size(folder_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    %[sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    load([folder_list{list_I,1},in_folder,input_name]);
    sub_num = match_num;
    sub_list = match_folder;   
    
    [M1,M2] = size(sub_list);
    im_folder = folder_list{list_I,1};
        
    for list_J = 1:M1
        sub_list{list_J,3}=sub_list{list_J,1}; % George
        data_name = [sub_list{list_J,3}(1:end-1),data_tail];
        raw_data = xlsread([im_folder,sub_folder,data_name]);
        raw_data = raw_data(sqrt(raw_data(:,2).*raw_data(:,3)) >= r_th,:);
        Inten_spot = raw_data(:,1).*raw_data(:,2).*raw_data(:,3)*2*pi;
        [n_hist,x_hist] = hist(Inten_spot,[hist_min:hist_bin:(hist_max+hist_bin)]);
        x_hist = x_hist(1:end-1);
        y_hist = n_hist(1:end-1)/sum(n_hist(1:end-1));
        mgau = @(a1,a2,a3,a4,b,c,d,x) a1*exp(-(x-b).^2./2./c.^2)+a2*exp(-(x-2*b).^2./2./2./c.^2)+a3*exp(-(x-3*b).^2./3./2./c.^2)+a4*exp(-(x-4*b).^2./4./2./c.^2)+d;
        spot_fit = fit( x_hist',y_hist',mgau, 'StartPoint', fit_initial);
        b = spot_fit.b;
        save([im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
        figure;maximize(gcf);
        bar(x_hist,y_hist,'b')
        hold on
        x_fit = [hist_min:(hist_max-hist_min)/1000:hist_max];
        y_fit = mgau(spot_fit.a1,spot_fit.a2,spot_fit.a3,spot_fit.a4,spot_fit.b,spot_fit.c,spot_fit.d,x_fit);
        plot(x_fit,y_fit,'r')
        xlabel('Spot intensity (A.U.)')
        ylabel('Frequency')
        title(['Spot intensity histogram fit: ',data_name,', b = ',num2str(b)],'Interpreter','none')
        saveas(gcf,[im_folder,sub_folder,data_name(1:find(data_name == '_',1,'last')-1),output_add,fig_tail])
        
        b = 0;
        if exist([im_folder,sub_folder_nullo]) ~= 7
            mkdir([im_folder,sub_folder_nullo]);
        end
        save([im_folder,sub_folder_nullo,data_name(1:find(data_name == '_',1,'last')-1),output_add,mat_tail],'spot_fit','x_hist','y_hist','b');
    end
end
end
   % sub
    function hist_on_GUIfun(p,imname,imfolder,imstack,imstack0,mask_stack,bk_use,mask_out)
%global imstack imstack0 immask dlayer dcontrast dcontrmin dfsize dimthresh dmask dpeak mstatus resolution0 L_ratio g vlayer vcontrast vcontrmin vfsize vimthresh vmask vpeak limsize limthresh limcontrast max_spot n0 xout0 n xout sg_name bk_use imname imfolder sub_folder mask_out m_data stack_th stack_th0 RNA_mask lim t0 dt0 mask_stack
%p,imname,imfolder,imstack,imstack0,mask_stack,bk_use,mask_out


n = zeros(0);
xout = zeros(0);
low_th = 0;
% lim = [0:10:500];
rxy = 10;
%set(handles.hist_on,'String','Wait')
%guidata(hObject, handles);
xls_tail = '.mat';
xls_name2 = [imname,'_raw',xls_tail];
k0 = strfind(imfolder,'stacks');
imfolder0 = imfolder(1:(k0(end)-1));

if p.type_single==1 || p.type_single2==1
    lim = [0:10:2000];
    decrease_th = 0.9^((lim(2)-lim(1))/20);
    [stack_th,t0,dt0] = th_find_low(imstack,spmask1(imstack,1),lim,decrease_th,low_th);
else
	lim = [500:100:8000];
    decrease_th = 0.85^((lim(2)-lim(1))/200);
    [stack_th,t0,dt0] = th_find(imstack,spmask1(imstack,1),lim,decrease_th,low_th);
end

% [~,temp_I] = findpeaks(dt0);
% stack_th = lim(temp_I(1));
% % [n_hist,~] = hist(imstack(mask_out)*65535,lim);
% % [~,temp_I] = findpeaks(-n_hist);
% % if isempty(temp_I)
% %     temp_I = 1;
% % end
% % stack_th = lim(temp_I(1));

% % % % if (stack_th > lim(end)) ||(stack_th < lim(1)) || isempty(stack_th) || isnan(stack_th)
% % % %     stack_th = stack_th0;
% % % % end
% % % % if get(handles.type_single,'Value')
%     stack_th = 200;
% % % % end
% stack_th = 2000;
% stack_th = 900;
% stack_th = 640;
% stack_th = 120;

sub_folder=[];
if p.presult_on==1 && exist([imfolder0,sub_folder,xls_name2],'file')
    [max_spot,~,~] = xlsread([imfolder0,sub_folder,xls_name2]);
%     max_spot = spfilter(max_spot,rxy,rxy,1);   %%% Fine filter to get rid of noise
else
    % for I_layer = 1:get(handles.layer_slide,'Max')
    %     immask(:,:,I_layer) = im2bw(g(:,:,I_layer),vimthresh);
    % end
    % STATS = regionprops(immask,imstack*65535,'MaxIntensity');
    % max_spot = imstack(imregionalmax(g)&immask)*65535;
    [mask_out,mask_out2D] = spmask1(imstack,1,stack_th);   %%% 3D local maximal mask generation (Primary filter)
    %mask_out = mask_out & (imstack*65535 >= 4000);   %%% Prefilter using peak intensity threshold
    % max_spot = imstack(mask_out)*65535;

    embryo_mask = repmat(imdilate(bwconvhull(max(logical(mask_stack),[],3)),strel('disk',20)),[1,1,size(mask_out,3)]);
    mask_out = mask_out & embryo_mask;
    mask_out2D = mask_out2D & embryo_mask;
%     mask_out(:,:,1:8) = false;
%     mask_out2D(:,:,1:8) = false;

    peak_parameter = ptrack(imstack0,imstack,mask_out,mask_out2D);   %%% Track/fit mRNA spots
    max_spot = spfilter(peak_parameter,rxy,rxy,1);   %%% Fine filter to get rid of noise
end
M_spot = max_spot(:,1).*max_spot(:,2).*max_spot(:,3)*2*pi;

RNA_mask = false(size(imstack,1),size(imstack,2),size(imstack,3));
for ii = 1:size(max_spot,1)
    x_spot = min(max(round(max_spot(ii,6)),1),size(imstack,1));
    y_spot = min(max(round(max_spot(ii,7)),1),size(imstack,2));
    RNA_mask(x_spot,y_spot,round(max_spot(ii,8))) = true;
end

%% Background calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p.show_bk==1 && p.background_on==1
    if isempty(bk_use)
        if exist([imfolder,sg_name],'file')
            [~,sg_record,~] = xlsread([imfolder0,sg_name]);
            I_record = find(strcmp(imname,sg_record(:,1)),1);
            if (~isempty(I_record)) && (size(sg_record,2) > 1)
                bk_use = sg_record(I_record,2:end);
            end
        end
    end
    if isempty(bk_use)
        sel_bk_Callback(hObject, eventdata, handles)
    end
end

n_bk = zeros(0);   %%% Initialization of background distribution
xout_bk = zeros(0);

if ~isempty(bk_use)
    for I_bk = 1:length(bk_use)
        if ~isempty(bk_use{I_bk})
            [bk_temp,~,~] = xlsread([imfolder0,sub_folder,bk_use{I_bk}]);
            n_bk = cat(2,n_bk,bk_temp(:,1));
        end
    end
    xout_bk = bk_temp(:,2);
    n_bk = mean(n_bk,2);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Histogram plot: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ---------------------------------------------------------------------
%axes(handles.hist_show);
figure;
%[n0,xout0] = hist([STATS.MaxIntensity],60);
[n0,xout0] = hist(M_spot,60);
%n0 = n0/sum(n0)*100;
bar(xout0,n0)
if p.show_bk==1 && ~isempty(n_bk)
    [~,Imax] = max(n_bk);
    %[~,Imin] = min(abs(xout0-xout_bk(Imax)));
    [~,Imin] = max(n0);
    nlevel = n0(Imin);
    hold on
    plot(xout_bk,n_bk/max(n_bk)*nlevel,'r')
    hold off
end
set(gca,'YScale','log')
ylabel('#')
ylim0 = ylim;
ylim0(1) = 1;
ylim(ylim0);
xlim0 = xlim;
%%%%% ---------------------------------------------------------------------
%axes(handles.sigma_show);
% plot(max_spot(:,1),max_spot(:,2),'ro',max_spot(:,1),max_spot(:,3),'gx',max_spot(:,1),max_spot(:,4),'k+');
% ylabel('Sigma (pixel)')
%set(gca,'YAxisLocation','right','Color','none')
% xlim(xlim0);
% legend('Sigma x','Sigma y','Sigma z')
%%% =======================================================================
%axes(handles.hist_zoom);
figure;
% low_th = 60000;
% low_bin = 500;
% low_th2 = 50000;
% low_bin2 = 5000;
low_th = 1.01e5;
low_bin = 1e3;
low_th2 = 1e5;
low_bin2 = 5e3;

[n,xout] = hist(M_spot,[0:low_bin:low_th]);
%[n,xout] = hist(M_spot(:,1),[0:1000:60000]);
%n = n/sum(n)*100;
bar(xout,n)
if p.show_bk==1 && ~isempty(n_bk)
    [~,Imax] = max(n_bk);
    %[~,Imin] = min(abs(xout-xout_bk(Imax)));
    [~,Imin] = max(n);
    nlevel = n(Imin);
    hold on
    plot(xout_bk,n_bk/max(n_bk)*nlevel,'r')
    hold off
end
xlim([0,low_th2]);
%xlim([0,10000]);
set(gca,'XTick',[0:low_bin2:low_th2])%,'YScale','log')
%set(gca,'XTick',[0:5000:50000],'YScale','log')

ylabel('#')
ylim0 = ylim;
if ylim0(1) > 1 
    ylim0(1) = 1;
end
ylim(ylim0);
xlim0 = xlim;

%%%%% Statistics:
m_temp = M_spot(M_spot(:,1)<low_th,1);
m_mean = mean(m_temp);
m_std = std(m_temp);
m_median = median(m_temp);
m_peak = xout(find(n == max(n(1:(end-1))),1));
m_p =  1./(48/(m_mean/m_std)^2+1);
m_data = [m_mean,m_std,m_median,m_peak,m_p];
%set(handles.text_input,'String',['Low intensity zoom in (mean = ',num2str(m_mean),', std = ',num2str(m_std),', median = ',num2str(m_median),', peak = ',num2str(m_peak),', p0 = ',num2str(m_p)])

%%%%% ---------------------------------------------------------------------
%axes(handles.sigma_zoom);
% plot(max_spot(:,1),max_spot(:,2),'ro',max_spot(:,1),max_spot(:,3),'gx',max_spot(:,1),max_spot(:,4),'k+');
% ylabel('Sigma (pixel)')
%set(gca,'YAxisLocation','right','Color','none')
% xlim(xlim0);
% legend('Sigma x','Sigma y','Sigma z')
%%%%% ---------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(handles.hist_on,'String','Histogram plot')
% Update handles structure
%guidata(hObject, handles);

end
    function [mask_out,varargout] = spmask(imstack,Nr,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% A function for spot recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mask_out: Output spot recognition result;
%% imstack: input image stack;
%% Nr: spot spread half width (on z direction);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max3 = imregionalmax(imstack);   %%% 3D local maxima
max2 = false(size(max3));   %%% 2D local maxima
max20 = max2;
N_layer = size(imstack,3);
for I_layer = 1:N_layer
    max20(:,:,I_layer) = imregionalmax(imstack(:,:,I_layer));
    max2(:,:,I_layer) = imdilate(max20(:,:,I_layer),strel('square',3));
end

Nr = ceil(Nr);
Ar = -Nr:Nr;
zmin = max(1,(1+Ar));
zmax = min(N_layer,(N_layer+Ar));
for Ir = 1:length(Ar)
    max3(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
end

mask_out = max3;

% -Edit by George 0919--------
%for I_layer = 1:N_layer
%    mask_out(:,:,I_layer) = max2(:,:,max(I_layer-1,1)) & max20(:,:,I_layer) & max2(:,:,min(I_layer+1,N_layer));
%end
% ----------------------------

if ~isempty(varargin)
    mask_out = mask_out & (imstack*65535 >= varargin{1});
    max20 = max20 & (imstack*65535 >= varargin{1});
end

varargout = {max20 & (~mask_out)};
    end
    function [mask_out,varargout] = spmask1(imstack,Nr,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% A function for spot recognition: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mask_out: Output spot recognition result;
%% imstack: input image stack;
%% Nr: spot spread half width (on z direction);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nr30 = 2;   %%% default 3D maxima extension range
max3 = imregionalmax(imstack);   %%% 3D local maxima
max3(:,:,[1,end]) = false;
max30 = false(size(max3));   %%% 2D local maxima of 3D local maxima residues
max2 = false(size(max3));   %%% 2D local maxima
max20 = max2;
N_layer = size(imstack,3);
for I_layer = 1:N_layer
    max20(:,:,I_layer) = imregionalmax(imstack(:,:,I_layer));
    max2(:,:,I_layer) = imdilate(max20(:,:,I_layer),strel('square',3));
end

Nr = ceil(Nr);
Ar = -Nr:Nr;
zmin = max(1,(1+Ar));
zmax = min(N_layer,(N_layer+Ar));
for Ir = 1:length(Ar)
    max3(:,:,zmin(Ir):zmax(Ir)) = max3(:,:,zmin(Ir):zmax(Ir)) & max2(:,:,zmin(end-Ir+1):zmax(end-Ir+1));
end

if length(varargin) > 1
    Nr3 = ceil(varargin{2});
else
    Nr3 = Nr30;
end
Ar3 = -Nr3:Nr3;
zmin3 = max(1,(1+Ar3));
zmax3 = min(N_layer,(N_layer+Ar3));
for Ir = 1:length(Ar3)
    max30(:,:,zmin3(Ir):zmax3(Ir)) = max30(:,:,zmin3(Ir):zmax3(Ir)) | imdilate(max3(:,:,zmin3(end-Ir+1):zmax3(end-Ir+1)),strel('square',3));
end
max30 = max20 & max30;

mask_out = max3;

%{
if ~isempty(varargin)
    mask_out = mask_out & (imstack*65535 >= varargin{1});
    max30 = max30 & (imstack*65535 >= varargin{1});
end
%}

varargout = {max30 & (~mask_out)};
    end
    function [overlay,g] = cf_show(max_image0,vlayer,psize0,L_ratio,g_threshold0,peak_on,vcontrast,vcontrmin,vmask,mask_out,varargin)

g = zeros(size(max_image0));
H = -fspecial('log',15,psize0/L_ratio);
for I_layer = 1:size(max_image0,3)
    g(:,:,I_layer) = imfilter(max_image0(:,:,I_layer),H,'replicate');
end
g(g<0) = 0;

if ~isempty(varargin)
    temp = varargin{1};
else
    temp = false;
end
overlay = cf_show1(mask_out,max_image0,g,vlayer,g_threshold0,peak_on,vcontrast,vcontrmin,vmask,temp);



% foci_bw = im2bw(g(:,:,vlayer),g_threshold0);
% foci_bw = imfill(foci_bw,'holes');
% %foci_bw0 = imdilate(foci_bw,strel('disk',4));
% foci_bw0 = bwmorph(foci_bw,'thicken',4);
% bw_perim_g = bwperim(foci_bw0);
% 
% if peak_on
%     bw_peak = imregionalmax(g(:,:,vlayer)) & foci_bw;
%     overlay(:,:,1) = (max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin)+double(bw_peak);
%     overlay(:,:,2) = vmask*double(bw_perim_g)+double(bw_peak);
%     overlay(:,:,3) = 0;
% else
%     overlay(:,:,1) = (max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin);
%     overlay(:,:,2) = vmask*double(bw_perim_g);
%     overlay(:,:,3) = 0;
% end
    end
    function overlay = cf_show1(mask_out,max_image0,g,vlayer,g_threshold0,peak_on,vcontrast,vcontrmin,vmask,varargin)

%global mask_out RNA_mask

foci_bw = im2bw(g(:,:,vlayer),g_threshold0);
foci_bw = imfill(foci_bw,'holes');
%%foci_bw0 = imdilate(foci_bw,strel('disk',4));
%foci_bw0 = bwmorph(foci_bw,'thicken',4);
%bw_perim_g = bwperim(foci_bw0);
bw_perim_g = foci_bw;
%peak_on = true;
if peak_on
    %bw_peak = imregionalmax(max_image0(:,:,vlayer));% & foci_bw;
    %bw_peak = imregionalmax(g);
    %bw_peak = bw_peak(:,:,vlayer);
    bw_peak = imdilate(mask_out(:,:,vlayer),strel('disk',1));
    RNA_mask=[];
    if ~isempty(RNA_mask) && ~isempty(varargin) && varargin{1}
        bw_peak2 = bwperim(imdilate(RNA_mask(:,:,vlayer),strel('disk',5)));
    else
        bw_peak2 = false(size(bw_peak));
    end
    
    overlay(:,:,1) = max((max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin),bw_peak);
    overlay(:,:,2) = max(vmask*double(bw_perim_g),bw_peak | bw_peak2);
    overlay(:,:,3) = double(bw_peak);
else
    overlay(:,:,1) = (max_image0(:,:,vlayer)-vcontrmin)/(vcontrast-vcontrmin);
    overlay(:,:,2) = vmask*double(bw_perim_g);
    overlay(:,:,3) = 0;
end
    end
    function [peak_parameter,varargout] = ptrack(imstack0,imstack,mask_out,varargin)
tic
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% A function to track/fit local maxima to gaussian shapes: %%%%%%%
%% peak_parameter: fitting parameter values for each local maxima with the 
%%                 form: (I0,sigma_x,sigma_y,sigma_z,A,x0,y0,z0,theta,residual)
%% imstack0: Raw image stack
%% imstack: Prefiltered image stack
%% mask_out: 3D local maxima from the primary filtering step
%% varargin: {mask_out2D} 2D-only maxima
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global all_offset xdata ydata xdata1 ydata1 raw_parameter isfit gau3 X0 Y0 Z0 xf fit_range Ipeak Nfit Jfit

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xrange = 3;  %%% Peak fitting range on x direction
yrange = 3;  %%% Peak fitting range on y direction
zrange = 0;  %%% Peak fitting range on z direction
sigmax0 = 2;   %%% Initial value of sigma_x
sigmay0 = 2;   %%% Initial value of sigma_y
sigmaz0 = 1;   %%% Initial value of sigma_z
crossxy = 0;   %%% Initial value of cross coefficient for "xy" terms
gau2_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,1)-x(',num2str(n),',3)).^2+x(',num2str(n),',4)*(xdata(:,2)-x(',num2str(n),',5)).^2)+x(',num2str(n),',6)*(xdata(:,1)-x(',num2str(n),',3)).*(xdata(:,2)-x(',num2str(n),',5)))+x(',num2str(n),',7)'];   %%% 2D single-gaussian model function text generator
gau1_gen = @(n) ['x(',num2str(n),',1)*exp(-(x(',num2str(n),',2)*(xdata(:,3)-x(',num2str(n),',3)).^2))+x(',num2str(n),',4)'];   %%% 1D single-gaussian model function text generator
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Coordinate matrix generation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(varargin)
    mask_out2D = varargin{1};
else
    mask_out2D = false(size(mask_out));
end
clear varargin
dim0 = size(imstack0);
pdim0 = prod(dim0);
%[XC,YC,ZC] = ndgrid(1:dim0(1),1:dim0(2),1:dim0(3));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak sorting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_out(:,:,[1,end]) = false;
Ipeak = find(mask_out);
[~,Itr] = sort(imstack(Ipeak),'descend');
Ipeak = Ipeak(Itr);
[X0,Y0,Z0] = ind2sub(dim0,Ipeak);   %%% Peak coordinates list
% X0 = XC(Ipeak);   %%% Peak coordinates x list
% Y0 = YC(Ipeak);   %%% Peak coordinates y list
% Z0 = ZC(Ipeak);   %%% Peak coordinates z list
vpeak = imstack0(Ipeak);   %%% Peak value list
isfit = false(size(X0));   %%% Peak fitting status

mask_out2D(:,:,[1,end]) = false;
Ipeak2D = find(mask_out2D);
[X02D,Y02D,Z02D] = ind2sub(dim0,Ipeak2D);   %%% Peak coordinates list
vpeak2D = imstack0(Ipeak2D);   %%% Peak value list
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Peak fitting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_parameter = zeros(length(Ipeak),10);
xf_all = zeros(0);
exstatus = false(length(Ipeak),1);
all_offset = zeros(pdim0,1);
% lb = [-1e9,-1e9,-1e9,-1e9,-1e9,-1e9,-1e9];
% ub = [1e9,1e9,1e9,1e9,1e9,1e9,1e9];

lb = [-inf,-inf,-inf,-inf,-inf,-inf,-inf];
ub = [inf,inf,inf,inf,inf,inf,inf];

options = optimset('Display','off');
IIfit = 0;
num_peak = length(Ipeak);
for Ifit = 1:num_peak
    if ~isfit(Ifit)   %%% check whether the peak has been fitted
        Jfit = Ifit;   %%% instant peak in process
        Nfit = [];   %%% total peaks in fitting
        Jfit2D = [];   %%% instant 2D peak in process
        Nfit2D = [];   %%% total 2D peaks in fitting
        %%% Fitting range/peak arrangement: %%%============================
%         fit_range = false(size(XC));
%         fit_range1 = false(size(XC));
        fit_range0 = zeros(0,3);
        while Jfit
            Nfit = cat(1,Nfit,Jfit);
            if (~isempty(Nfit2D)) && (~isempty(Jfit2D))
                Nfit2D = cat(1,Nfit2D,Jfit2D);
            end
            
            for I_temp = 1:length(Jfit)
                [Xtemp,Ytemp,Ztemp] = ndgrid(max((X0(Jfit(I_temp))-2*xrange),1):min((X0(Jfit(I_temp))+2*xrange),dim0(1)),max((Y0(Jfit(I_temp))-2*yrange),1):min((Y0(Jfit(I_temp))+2*yrange),dim0(2)),max((Z0(Jfit(I_temp))-2*zrange),1):min((Z0(Jfit(I_temp))+2*zrange),dim0(3)));
                fit_range0 = union(fit_range0,[Xtemp(:),Ytemp(:),Ztemp(:)],'rows');
            end
%             fit_range([max(1,X0(Jfit)-xrange):min(dim0(1),X0(Jfit)+xrange)],[max(1,Y0(Jfit)-yrange):min(dim0(2),Y0(Jfit)+yrange)],Z0(Jfit)) = true;
%             fit_range1(X0(Jfit),Y0(Jfit),[max(1,Z0(Jfit)-zrange):min(dim0(3),Z0(Jfit)+zrange)]) = true;
%             [~, ~, Jpeak] = intersect(find(fit_range),Ipeak);
            [~,~,Jpeak] = intersect(fit_range0,[X0,Y0,Z0],'rows');
            Jfit = setdiff(setdiff(Jpeak,Nfit),find(isfit));
            if ~isempty(X02D)
                [~,~,Jpeak2D] = intersect(fit_range0,[X02D,Y02D,Z02D],'rows');
                Jfit2D = setdiff(Jpeak2D,Nfit2D);
            else
                Jpeak2D = zeros(0);
                Jfit2D = zeros(0);
            end
        end
        Nfit2D = cat(1,Nfit2D,Jfit2D);
       
        %%% ===============================================================
        
        %%% Fitting preparation: %%%=======================================
        %%%%% 3D multi-gaussian fit function generation:
%         textfun = 'gau3 = @(x,xdata) ';
%         for n = 1:length(Nfit)
%             textfun = [textfun,gau3_gen(n),'+'];
%         end
        textfun = 'gau2 = @(x,xdata) ';
        textfun1 = 'gau1 = @(x,xdata) ';
        for n = 1:(length(Nfit)+length(Nfit2D))
            textfun = [textfun,gau2_gen(n),'+'];
            textfun1 = [textfun1,gau1_gen(n),'+'];
        end
        textfun = [textfun(1:(end-1)),';'];
        textfun1 = [textfun1(1:(end-1)),';'];
        eval(textfun);
        eval(textfun1);
        %%%%% Data points collection: 
%         range1D = find(fit_range);
%         range11D = find(fit_range1);
%         xdata = [XC(range1D),YC(range1D),ZC(range1D)];   %%%%% Collect the coordinates of data points
%         ydata = imstack0(range1D)-all_offset(range1D);   %%%%% Collect the intensity values of data points
        xdata = zeros(0,3);
        ydata = zeros(0,3);
        for I_temp = 1:length(Nfit)
            [Xtemp,Ytemp,Ztemp] = ndgrid(max((X0(Nfit(I_temp))-xrange),1):min((X0(Nfit(I_temp))+xrange),dim0(1)),max((Y0(Nfit(I_temp))-yrange),1):min((Y0(Nfit(I_temp))+yrange),dim0(2)),max((Z0(Nfit(I_temp))-zrange),1):min((Z0(Nfit(I_temp))+zrange),dim0(3)));
            xdata = union(xdata,[Xtemp(:),Ytemp(:),Ztemp(:)],'rows');
        end
        for I_temp = 1:length(Nfit2D)
            [Xtemp,Ytemp,Ztemp] = ndgrid(max((X02D(Nfit2D(I_temp))-xrange),1):min((X02D(Nfit2D(I_temp))+xrange),dim0(1)),max((Y02D(Nfit2D(I_temp))-yrange),1):min((Y02D(Nfit2D(I_temp))+yrange),dim0(2)),max((Z02D(Nfit2D(I_temp))-zrange),1):min((Z02D(Nfit2D(I_temp))+zrange),dim0(3)));
            xdata = union(xdata,[Xtemp(:),Ytemp(:),Ztemp(:)],'rows');
        end
        ydata =imstack0(sub2ind(dim0,xdata(:,1),xdata(:,2),xdata(:,3)));
%         xdata1 = [XC(range11D),YC(range11D),ZC(range11D)];   %%%%% Collect the coordinates of data points
%         ydata1 = imstack(range11D)-all_offset(range11D);   %%%%% Collect the intensity values of data points

        %%%%% Fitting parameter initialization:
        xf0 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit)), X0(Nfit), 1./2./sigmay0.^2.*ones(size(Nfit)), Y0(Nfit), crossxy.*ones(size(Nfit)), min(ydata)*ones(size(Nfit))];
        xf10 = [vpeak(Nfit)-all_offset(Ipeak(Nfit))-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit)), Z0(Nfit), min(ydata)*ones(size(Nfit))];

        xf02D = [vpeak2D(Nfit2D)-all_offset(Ipeak2D(Nfit2D))-min(ydata), 1./2./sigmax0.^2.*ones(size(Nfit2D)), X02D(Nfit2D), 1./2./sigmay0.^2.*ones(size(Nfit2D)), Y02D(Nfit2D), crossxy.*ones(size(Nfit2D)), min(ydata)*ones(size(Nfit2D))];
        xf102D = [vpeak2D(Nfit2D)-all_offset(Ipeak2D(Nfit2D))-min(ydata), 1./2./sigmaz0.^2.*ones(size(Nfit2D)), Z02D(Nfit2D), min(ydata)*ones(size(Nfit2D))];
        %%% ===============================================================
        
        %%% Multi-Gaussian fitting: %%%====================================
        [xf,resnorm,~,exitflag,~] = lsqcurvefit(gau2,[xf0;xf02D],xdata,ydata,lb,ub,options);
        %[xf,resnorm,~,exitflag,~] = lsqcurvefit(gau2,xf0,xdata,ydata);
        %[xf1,~,~,exitflag1,~] = lsqcurvefit(gau1,xf10,xdata1,ydata1);
        xf1 = xf10;
        xf_all = cat(1,xf_all,xf);
        xf = xf(1:length(Nfit),:);
        exitflag1 = 1;
        raw_parameter(Nfit,:) = [xf(:,1:6),xf1(:,2:3),xf(:,7),sqrt(resnorm./length(ydata)*ones(size(Nfit)))];
        %raw_parameter(Nfit,:) = [xf(:,1:6),repmat(Z0(Ifit),size(xf,1),1),xf1(:,3),xf(:,7),sqrt(resnorm./length(ydata)*ones(size(Nfit)))];
        exstatus(Nfit) = (exitflag > 0) & (exitflag1 > 0) & ((xf(:,2)+xf(:,4)) >= 0) & (4*xf(:,2).*xf(:,4) >= xf(:,6).*xf(:,6)) & (xf(:,3) > 0) & (xf(:,3) < dim0(1)) & (xf(:,5) > 0) & (xf(:,5) < dim0(2)) & (xf1(:,3) > 0) & (xf1(:,3) < dim0(3));   %%% exit status
        %xf_offset = [xf(:,1:8),zeros(size(xf,1),1)];
        %all_offset = all_offset+gau3(xf_offset,[XC([1:pdim0]'),YC([1:pdim0]'),ZC([1:pdim0]')]);
        %%% ===============================================================
        isfit(Nfit) = true;
    end
    IIfit = IIfit+1;
    if IIfit >= 100
        disp(['Finishing: ',num2str(Ifit),'/',num2str(num_peak),', time: ',num2str(toc)])
        IIfit = 0;
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fitting parameter output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_parameter = raw_parameter(exstatus,:);
peak_parameter(:,[1,5:8,10]) = raw_parameter(:,[1,9,3,5,8,10]);
peak_parameter(:,2) = 1./sqrt(raw_parameter(:,2)+raw_parameter(:,4)+sqrt((raw_parameter(:,2)-raw_parameter(:,4)).^2+raw_parameter(:,6).^2));
peak_parameter(:,3) = 1./sqrt(raw_parameter(:,2)+raw_parameter(:,4)-sqrt((raw_parameter(:,2)-raw_parameter(:,4)).^2+raw_parameter(:,6).^2));
peak_parameter(:,4) = 1./sqrt(2*raw_parameter(:,7));
peak_parameter(:,9) = atan2(raw_parameter(:,6),(raw_parameter(:,2)-raw_parameter(:,4)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout = {xf_all};
%clear XC YC ZC
toc
    end
    function max_spot = spfilter(peak_parameter,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% A function to filter out low intensity/large eccentricity noise: %%%
%% max_spot: spot characteristics (Intensity,sigma_x,sigma_y,sigma_z)
%% peak_parameter: fitting parameter values for each local maxima with the 
%%                 form: (I0,sigma_x,sigma_y,sigma_z,A,x0,y0,z0,theta,residual)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmin = 0.4;

true_spot = true(size(peak_parameter,1),1);
true_spot = true_spot & (peak_parameter(:,1) <= 1e5) & (peak_parameter(:,1) > 0);
true_spot = true_spot & (peak_parameter(:,2) > 0) & (peak_parameter(:,3) > 0);
if isempty(varargin) || length(varargin) < 3
    true_spot = true_spot & (sqrt(1-(peak_parameter(:,2)./peak_parameter(:,3)).^2) <= 0.9);
else
    true_spot = true_spot & (sqrt(1-(peak_parameter(:,2)./peak_parameter(:,3)).^2) <= varargin{3});
end

if isempty(varargin)
    true_spot = true_spot & (peak_parameter(:,2) >= rmin) & (peak_parameter(:,2) <= 3);
    true_spot = true_spot & (peak_parameter(:,3) >= rmin) & (peak_parameter(:,3) <= 3);
else
    true_spot = true_spot & (peak_parameter(:,2) >= rmin) & (peak_parameter(:,2) <= varargin{1});
    true_spot = true_spot & (peak_parameter(:,3) >= rmin) & (peak_parameter(:,3) <= varargin{2});
end
%true_spot = true_spot & (peak_parameter(:,4) > 0.25) & (peak_parameter(:,3) <2);


%max_spot(:,1) = peak_parameter(true_spot,1);
%max_spot(:,1) = peak_parameter(:,1).*peak_parameter(:,2).*peak_parameter(:,3)*sqrt(2*pi);
%max_spot(:,2:10) = peak_parameter(true_spot,2:10);
max_spot = peak_parameter(true_spot,:);
    end
    function mask_label = label3D(raw3D)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A function to label the stitched 2D segmentation mask %%%%%%%%%%%%%%%%%%
%% raw3D: stitched 2D mask.                              %%%%%%%%%%%%%%%%%%
%% mask_label: labeled 3D mask.                          %%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_label = zeros(size(raw3D));
max_label = 0;   %%% maximal label #
if size(raw3D,3) <= 1
   thick_thresh = 1;   %%% thickness threshold
else
    thick_thresh = 3;   %%% thickness threshold
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% For the first layer: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_label(:,:,1) = bwlabel(raw3D(:,:,1));
max_label = max(mask_label(:));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2D mask stitching: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I_layer = 1:(size(raw3D,3)-1)
    temp_label1 = mask_label(:,:,I_layer);
    temp_label2 = zeros(size(raw3D(:,:,1)));
    bwlabel1 = bwlabel(raw3D(:,:,I_layer));
    bwlabel2 = bwlabel(raw3D(:,:,I_layer+1));
    prop1 = regionprops(raw3D(:,:,I_layer),'Centroid');
    prop2 = regionprops(raw3D(:,:,I_layer+1),'Centroid');
    
    if ~isempty(prop1)
        xy_center0 = round(cell2mat({prop1.Centroid}'));   %%% center coordinates matrix for layer I
        ind_center0 = sub2ind(size(bwlabel1),xy_center0(:,2),xy_center0(:,1));   %%% center linear indices matrix for layer I
        bw2temp_prop = regionprops(bwlabel1,temp_label1,'MaxIntensity');   %%% convert bwlabel1 labels to temp_label1 labels
        bw2temp = [bw2temp_prop.MaxIntensity];
        
        for I_area = 1:length(prop2)

            x_center = round(prop2(I_area).Centroid(1));
            y_center = round(prop2(I_area).Centroid(2));
            if temp_label1(y_center,x_center)
                temp_label2(bwlabel2 == I_area) = temp_label1(y_center,x_center);
            end

            ind_list = find(bwlabel2(ind_center0) == I_area);
            current_ind = temp_label2(find(bwlabel2 == I_area,1));
            if (~current_ind) && (~isempty(ind_list))
                current_ind = min(bw2temp(ind_list));
                temp_label2(bwlabel2 == I_area) = current_ind;
            elseif (~current_ind) && isempty(ind_list);
                max_label = max_label+1;
                current_ind = max_label;
                temp_label2(bwlabel2 == I_area) = current_ind;
            end
            ind_label = setdiff(bw2temp(ind_list),current_ind);
            if ~isempty(ind_label)
                temp_label1(ismember(temp_label1,ind_label)) = current_ind;
                mask_label(ismember(mask_label,ind_label)) = current_ind;
                bw2temp(ind_list) = current_ind;
            end
        end
        mask_label(:,:,I_layer+1) = temp_label2;
    else
        mask_label(:,:,I_layer+1) = (bwlabel2+max_label).*raw3D(:,:,I_layer+1);
        max_label = max(mask_label(:));
    end
%     figure; imshow(mask_label(:,:,1));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Remove the fake regions that are too thin: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thickness_region = zeros(1,max(mask_label(:)));
for I_layer = 1:size(mask_label,3)
    region_plus = setdiff(unique(mask_label(:,:,I_layer)),[0]);
    thickness_region(region_plus) = thickness_region(region_plus)+1;
end
thin_region = find(thickness_region < thick_thresh);
mask_label(ismember(mask_label,thin_region)) = 0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Resign the label #: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_label = zeros(1,max_label+1);
new_label(sort(unique(mask_label))+1) = [0:(length(unique(mask_label))-1)];
mask_label = new_label(mask_label+1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    function [stack_th,t0,dt0] = th_find(imstack,mask_out,lim,decrease_th,varargin)

if ~isempty(varargin)
    low_th = varargin{1};
else
    low_th = 0;
end
t0 = zeros(size(lim));
all_spot = imstack(mask_out)*65535;
for I_th = 1:length(lim)
    t0(I_th) = nnz(all_spot >=lim(I_th));
end
dt0 = exp(diff(log(t0)));
dt0(end+1) = nan;
stack_th1 = lim(find((dt0 == max(dt0(lim >= low_th)) & (lim >= low_th)),1));

[~,I0] = min(dt0(t0 > 10));
I1 = find((dt0(I0:end) >= decrease_th) & (lim(I0:end) >= low_th),1)+I0-1;
stack_th2 = lim(I1);

if ~isempty(stack_th2)
    stack_th = min(stack_th1,stack_th2);
else
    stack_th = stack_th1;
end
    end
    function [stack_th,t0,dt0] = th_find_low(imstack,mask_out,lim,decrease_th,varargin)

if ~isempty(varargin)
    low_th = varargin{1};
else
    low_th = 0;
end
t0 = zeros(size(lim));
all_spot = imstack(mask_out)*65535;
for I_th = 1:length(lim)
    t0(I_th) = nnz(all_spot >=lim(I_th));
end
dt0 = diff(t0);
dt0(end+1) = nan;

H = fspecial('gaussian',10,4);
[~,I0] = findpeaks(conv(dt0,H(5,:)/sum(H(5,:)),'same'));
[~,I1] = findpeaks(conv(diff(conv(dt0,H(5,:)/sum(H(5,:)),'same')),H(5,:)/sum(H(5,:)),'same'));
[~,I2] = findpeaks(-conv(diff(conv(dt0,H(5,:)/sum(H(5,:)),'same')),H(5,:)/sum(H(5,:)),'same'));
lim1 = lim(I1(find(lim(I1) >= low_th,1))+2);
if isempty(lim1)
    lim1 = low_th;
end
lim2 = lim(I2(find(lim(I2) >= lim1,1))+2);
if ~isempty(lim2)
    lim0 = lim(I0(find(lim(I0) >= lim1 & lim(I0) <= lim2,1))+1);
else
    lim0 = zeros(0);
end

if ~isempty(lim0)
    stack_th = lim0;
else
    stack_th = lim1;
end
end



% ===== Step5. IF spot analysis =====


% ===== Step6. Analysis of transcriptional regulation =====
function Dualprocess6()
clear all
close all

tic
%% Parameter setting: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_name = 'Duallist.mat';
in_folder = 'stacks/';
input_name = 'matchlist.mat';
mismatch_name = 'mismatch.xls';
mismatch_name2 = 'mismatch_60X.xls';
xymismatch_name = 'xymismatch.mat';
xymismatch_name2 = 'xymismatch_60X.mat';
flip_list = {'Cad'};
image_type = '*.tif';
mask_folder = 'masks/';
mask_name = 'mask.mat';
out_folder = 'Results/';
% out_folder0 = 'Results_3Dz1/';
out_folder0 = '';
hist_folder = 'Histogram/';
fit_folder = 'Histogram_A/';
hist_folder2 = 'Histogram_RNA2/';
fit_folder2 = 'Histogram_default/';
fit_add = '_spot_fit';
hist_tail = '_raw.xls';
N_thresh = 3;
output_tail = '.xls';
figure_tail = '.fig';
mat_tail = '.mat';
seg_add = '_seg';
th_add = '_th';
foci_add = '_foci';
fish_add = '_fish';
num_add = '_num';
int_add = '_int';
cmp_add = '_cmp';
sel_add = '_sel';
protein_add = '_protein';
RNA_add = '_RNA';
signal2_add = '_RNA2';
nu_add = '_nucleus';
cyto_add = '_cytoplasm';
fate_add = '_fate';
reg_add = '_regulation';
fluc_add = '_fluc';
bino_add = '_bino';
local_add = '_local';
hist_add = '_hist';
fake_add = '_fake';
abs_add = '_abs';
rad_add = '_rad';
D3_add = '_3D';
compare_add = '_compare';
DAPI_add = '_DAPI';
noise_add = '_noise';
sub_pos = [3,3];
cycle_pos0 = zeros(1,sub_pos(1));
cycle_range = [11,12,13];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data loading: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_list, folder_list] = xlsread(list_name);
% folder_list = folder_list(strcmpi('b',folder_list(:,6)) | strcmpi('ba',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3',folder_list(:,6)),:);
folder_list = folder_list(strcmpi('T',folder_list(:,6)),:);
% folder_list = folder_list(strcmpi('b3h',folder_list(:,6)),:);
[N1,N2] = size(folder_list);
% run_list = {'[1]','[2]'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Image analysis: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for list_I = 1:N1
    [sub_num, sub_list] = xlsread([folder_list{list_I,1},in_folder,input_name]);
    [M1,M2] = size(sub_list);
    if ~isempty(out_folder0)
        copyfile([folder_list{list_I,1},out_folder],[folder_list{list_I,1},out_folder0]);
    end
    
    for list_J = 1:M1%eval(run_list{list_I})
        if isempty(strfind(sub_list{list_J,3},'_60X'))
            [~,~,mismatch_matrix] = xlsread(mismatch_name);
            load(xymismatch_name);
        else
            [~,~,mismatch_matrix] = xlsread(mismatch_name2);
            load(xymismatch_name2);
        end
        
        flip_axis = any(cellfun(@(x) ~isempty(strfind(folder_list{list_I,1},x)),flip_list));   %%% Check whether the protein profile has an opposite direction (to Bcd)
        
        image_folder = [folder_list{list_I,1},in_folder,sub_list{list_J,3}];
        result_folder = [folder_list{list_I,1},out_folder];
        resolutionz = sub_num(list_J,11);
        signal2_channel = sub_num(list_J,12);
        Nbin = sub_num(list_J,2);
        Mdim = sub_num(list_J,3);
        load([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail]);

        load([folder_list{list_I,1},mask_folder,sub_list{list_J,3},mask_name]);   %%% load 3D mask
        z_size = size(mask_stack,3);
        load([folder_list{list_I,1},fit_folder,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh = b*N_thresh;   %%% set foci intensity threshold
        single_Inten = b;
        load([folder_list{list_I,1},fit_folder2,sub_list{list_J,3}(1:end-1),fit_add,mat_tail]);   %%% load spot intensity fitting result
        Inten_thresh2 = b*N_thresh;   %%% set foci intensity threshold
        N_cycle = round(log2(max(mask_stack(:))))+2;   %%% Calculate the nuclear cycle number
%         N_cycle = sub_num(list_J,13);

        all_color = eval(folder_list{list_I,5});
        protein_color = all_color{protein_channel};
        RNA_color = all_color{RNA_channel};
        protein_RNA_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(RNA_color,mismatch_matrix(:,1))};
        protein_RNA_xymismatch = {eval([proteindual_profile_color,'_',RNA_color]),eval([protein_color,'_',RNA_color,'_con']),eval([protein_color,'_',RNA_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};
        signal2_color = all_color{signal2_channel};
        protein_signal2_mismatch = mismatch_matrix{strcmp(protein_color,mismatch_matrix(1,:)),strcmp(signal2_color,mismatch_matrix(:,1))};
        protein_signal2_xymismatch = {eval([protein_color,'_',signal2_color]),eval([protein_color,'_',signal2_color,'_con']),eval([protein_color,'_',signal2_color,'_x0']),Nbin,Mdim,resolution,resolution_mismatch};

        [signal_stack,RNA_stack,DAPI_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,RNA_channel,image_folder,protein_RNA_mismatch);   %%% load 3D image stacks
        [nucleus_DAPI_profile,DNA_mask] = DAPI_profile3D(mask_stack,DAPI_stack,image_folder,N_cycle);
        clear DAPI_stack
        
        
        [foci_bw3D,max_image00,SS] = modified_foci3D([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,Inten_thresh);
        foci_bw2D = max(foci_bw3D,[],3);
        max_image00 = max_image00/single_Inten;
        EL_info = get_EL(em_mask);   %%% get EL information (extreme points, EL length) from the embryo mask
        [nucleus_protein_profile,cytoplasmic_protein_profile,quanti_p,nucleus_protein_ab,flip_EL] = protein_profile3D(nucleus_DAPI_profile,max_image,EL_info,protein_channel,mask_stack,signal_stack,image_folder,N_cycle,resolution,resolutionz,[],[],[],[],[],[],flip_axis);

%         [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,cyto_bw,foci_bw3D,max_image00,EL_info,SS,RNA_stack,resolution,image_folder,N_cycle,[],flip_EL);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
        [nucleus_RNA_profile,foci_RNA_profile,cytoplasmic_RNA_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,em_mask,foci_bw3D,max_image00,EL_info,SS,RNA_stack,resolution,image_folder,N_cycle,[],flip_EL);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used

%         signal_stack = (signal_stack+quanti_p(2))/quanti_p(1);
%         nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),(nucleus_protein_profile(:,2)+quanti_p(2))/quanti_p(1),nucleus_DAPI_profile];
%        nuclei_variance = [0;nucleus_protein_profile(:,3)];
%        variance_stack = nuclei_variance(mask_stack+1);
        nucleus_protein_profile_ab = [nucleus_protein_profile(:,1),nucleus_protein_ab,nucleus_protein_profile(:,4),nucleus_protein_profile(:,4)];
%         pro_temp = [0;nucleus_protein_ab];
%         signal_stack = pro_temp(mask_stack+1);
        signal_stack = (signal_stack+quanti_p(2))/quanti_p(1);
%         signal_stack = signal_stack/quanti_p(1);
        nucleus_protein_profile = [nucleus_protein_profile,nucleus_DAPI_profile];
        RNA_DNA(nucleus_protein_profile(:,5),nucleus_RNA_profile(:,[3,4]),image_folder,N_cycle,100);
        TXnoise_plot(nucleus_RNA_profile,N_cycle,image_folder,[101,102]);
        
        [foci_bw3D,max_image00,SS] = modified_foci3D2([folder_list{list_I,1},hist_folder,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,RNA_channel,mask_stack,N_cycle,image_folder,protein_RNA_mismatch,protein_RNA_xymismatch,Inten_thresh);
        max_image00 = max_image00/single_Inten;
        foci_bw00 = find(foci_bw3D);
        max_image00 = max_image00.*SS;
        max_image_list = max_image00(foci_bw00);
        foci_mask = logical(max(foci_bw3D,[],3));
        clear foci_bw3D max_image00 SS 
        [foci_data,fake_data,h,t_absolute,r_size] = dual_local_local3D3(foci_bw00,max_image_list,mask_stack,DNA_mask,signal_stack,RNA_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz);
        clear RNA_stack
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if exist([folder_list{list_I,1},hist_folder2])
            [~,signal2_stack] = stack3D(imclearborder(seg_bw),protein_channel,DAPI_channel,signal2_channel,image_folder,protein_signal2_mismatch);   %%% load 3D image stacks

            [foci_bw3D2,max_image002,SS2] = modified_foci3D([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,signal2_channel,mask_stack,N_cycle,image_folder,protein_signal2_mismatch,Inten_thresh2,44);
            [nucleus_signal2_profile,foci_signal2_profile,cytoplasmic_signal2_profile] = RNA_profile3D(mask_stack,nucleus_DAPI_profile,em_mask,foci_bw3D2,max_image002,EL_info,SS2,signal2_stack,resolution,image_folder,N_cycle,[45,46,47,48,49,411],flip_EL);   %%% Input the foci area matrix from modified_foci.m if gaussian fitting method is used
            nucleus_signal2_profile(:,1) = nucleus_protein_profile(:,1);

            [foci_bw3D2,max_image002,SS2] = modified_foci3D2([folder_list{list_I,1},hist_folder2,sub_list{list_J,3}(1:(end-1)),hist_tail],max_image,signal2_channel,mask_stack,N_cycle,image_folder,protein_signal2_mismatch,protein_signal2_xymismatch,Inten_thresh2,44);
            foci_bw002 = find(foci_bw3D2);
            max_image002 = max_image002.*SS2;
            max_image_list2 = max_image002(foci_bw002);
            clear foci_bw3D2 max_image002 SS2 
            [foci_data2,fake_data2,h2,t_absolute2,r_size2] = dual_local_local3D3(foci_bw002,max_image_list2,mask_stack,DNA_mask,signal_stack,signal2_stack,nucleus_protein_profile_ab,image_folder,N_cycle,resolution,resolutionz,foci_mask,[473:477,480,481]);
            clear mask_stack signal_stack signal2_stack% nucleus_protein_profile_ab
        end

            dual_profile(nucleus_protein_profile_ab,nucleus_RNA_profile,foci_RNA_profile,image_folder,N_cycle,z_size,{'M','#'});
            %enrichment_compare(foci_data,fake_data,r_size,foci_data2,fake_data2,signal2_add(2:end),image_folder,N_cycle,true)
            
%%% Output: %%%============================================================
        result_folder = [folder_list{list_I,1},out_folder];
        if exist(result_folder) ~= 7
            mkdir(result_folder);
        end
        %saveas(51,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,seg_add,figure_tail]);
        saveas(2,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,D3_add,figure_tail]);
        %saveas(3,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,figure_tail]);
        saveas(4,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,seg_add,D3_add,figure_tail]);
        saveas(5,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,int_add,D3_add,figure_tail]);
        saveas(6,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,num_add,D3_add,figure_tail]);
        saveas(7,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,num_add,D3_add,figure_tail]);
        saveas(8,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,figure_tail]);
        saveas(9,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,int_add,D3_add,cmp_add,figure_tail]);
        %saveas(10,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,th_add,sel_add,figure_tail]);
        saveas(11,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,fate_add,D3_add,figure_tail]);
        
        saveas(12,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,reg_add,D3_add,figure_tail]);
        saveas(13,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,reg_add,D3_add,figure_tail]);
        saveas(14,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fate_add,D3_add,figure_tail]);
        saveas(15,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,foci_add,reg_add,D3_add,figure_tail]);
        saveas(16,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,DAPI_add,D3_add,figure_tail]);
        
        saveas(62,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,fluc_add,D3_add,figure_tail]);
        saveas(63,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,bino_add,D3_add,figure_tail]);
        
        saveas(100,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),DAPI_add,RNA_add,figure_tail]);
        saveas(101,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,noise_add,figure_tail]);

        saveas(73,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,D3_add,figure_tail]);
        saveas(74,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,D3_add,figure_tail]);
        saveas(75,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,D3_add,figure_tail]);
        saveas(76,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,figure_tail]);
        saveas(77,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,D3_add,figure_tail]);
        saveas(80,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,D3_add,figure_tail]);
        saveas(81,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,figure_tail]);

        if exist([folder_list{list_I,1},hist_folder2])
            saveas(44,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,seg_add,D3_add,figure_tail]);
            saveas(45,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,int_add,D3_add,figure_tail]);
            saveas(46,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),fish_add,signal2_add,num_add,D3_add,figure_tail]);
            saveas(47,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,num_add,D3_add,figure_tail]);
            saveas(48,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,figure_tail]);
            saveas(49,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),foci_add,signal2_add,int_add,D3_add,cmp_add,figure_tail]);
            saveas(411,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),nu_add,signal2_add,fate_add,D3_add,figure_tail]);

            saveas(473,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,D3_add,signal2_add,figure_tail]);
            saveas(474,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,hist_add,D3_add,signal2_add,figure_tail]);
            saveas(475,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,RNA_add,D3_add,signal2_add,figure_tail]);
            saveas(476,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,signal2_add,figure_tail]);
            saveas(477,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,hist_add,D3_add,signal2_add,figure_tail]);
            saveas(480,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,rad_add,D3_add,signal2_add,figure_tail]);
            saveas(481,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,signal2_add,figure_tail]);
        end

        %saveas(576,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,foci_add,D3_add,compare_add,figure_tail]);
        %saveas(581,[result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,abs_add,rad_add,D3_add,compare_add,figure_tail]);

        %save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','-append');
        save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'image_folder','image_type','WGA_channel','DAPI_channel','RNA_channel','protein_channel','resolution','seg_bw','cyto_bw','max_image','foci_bw2D','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','h','t_absolute','r_size','resolutionz','quanti_p','nucleus_protein_profile_ab','-append');

        if exist([folder_list{list_I,1},hist_folder2])
            save([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),mat_tail],'foci_data2','fake_data2','h2','t_absolute2','r_size2','-append')
        end
        
        if ~isempty(nucleus_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,nu_add,output_tail],nucleus_RNA_profile);
        end
        if ~isempty(cytoplasmic_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,cyto_add,output_tail],cytoplasmic_RNA_profile);
        end
        if ~isempty(foci_RNA_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),RNA_add,foci_add,output_tail],foci_RNA_profile);
        end
        if ~isempty(nucleus_protein_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,nu_add,output_tail],nucleus_protein_profile);
        end
        if ~isempty(cytoplasmic_protein_profile)
           xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,cyto_add,output_tail],cytoplasmic_protein_profile);
        end
        for I_data = 1:length(foci_data)
            if ~isempty(foci_data{I_data})
                xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,foci_add,output_tail],foci_data{I_data},I_data);
            end
        end
        for I_data = 1:length(fake_data)
            if ~isempty(fake_data{I_data})
                xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,fake_add,foci_add,output_tail],fake_data{I_data},I_data);
            end
        end
        
        if exist([folder_list{list_I,1},hist_folder2])
            if ~isempty(nucleus_signal2_profile)
               xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,nu_add,output_tail],nucleus_signal2_profile);
            end
            if ~isempty(cytoplasmic_signal2_profile)
               xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,cyto_add,output_tail],cytoplasmic_signal2_profile);
            end
            if ~isempty(foci_signal2_profile)
               xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),signal2_add,foci_add,output_tail],foci_signal2_profile);
            end
            for I_data = 1:length(foci_data2)
                if ~isempty(foci_data2{I_data})
                    xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,signal2_add,foci_add,output_tail],foci_data2{I_data},I_data);
                end
            end
            for I_data = 1:length(fake_data2)
                if ~isempty(fake_data2{I_data})
                    xlswrite([result_folder,sub_list{list_J,3}(1:(length(sub_list{list_J,3})-1)),protein_add,local_add,signal2_add,fake_add,foci_add,output_tail],fake_data2{I_data},I_data);
                end
            end
        end

        
        sub_num(list_J,13) = N_cycle;
        
        clear ('seg_bw','cyto_bw','max_image','foci_bw','nucleus_RNA_profile','foci_RNA_profile','cytoplasmic_RNA_profile','nucleus_protein_profile','cytoplasmic_protein_profile','N_cycle','foci_data','fake_data','foci_signal2_profile','cytoplasmic_signal2_profile','foci_data2','fake_data2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end
%     try
%         xlswrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
%     catch
% %         xlwrite([folder_list{list_I,1},in_folder,input_name],cat(2,sub_list,num2cell(sub_num)));
%     end
end
toc


end


% ===== Step7. Fitting the nascent mRNA distribution =====