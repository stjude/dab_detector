%% Detect DAB spots:
% [1] convert czi images to tif 
% [2] run ilastik pixel classification
% [3] post-process to generate final DAB spot segmentation
% [4] split images into top/bottom to separate hemispheres
% [4b] split images manually
% [5] analysis to determine(DAB spot area)/(total segmented area)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%% ----- Modify lines in this section --------%%%%%%%%%
dir_src = '<directory-with-images-to-be-processed>';
dir_tif = '<directory-to-save-converted-tif-files>';  % target dir for tif files
dir_ilastik_out = '<directory-to-store-ilastik-segmentation>'; % target ilastik out
dir_final_seg_out = '<directory-to-store-processed-segmentation>'; % target dir for final segmentation for analysis

% Configure ilastik
ilastik_command_headless = 'LAZYFLOW_THREADS=1 LAZYFLOW_TOTAL_RAM_MB=20000 <path-to-run_ilastik.sh> --headless --readonly READONLY'; % where is ilastik?
ilastik_project_1 = '<ILP-classifier-location>'; % pretrained ilastik project: at same scale and for isotropic voxels

% Tell MATLAB where your Bio-Formats Toolbox installation
dir_bfmatlab = '<directory-containing-bfmatlab-including-jar-file>';

% Determine which steps of the pipeline to run: set these variables to 1 to run, 0 to omit
do_conversion = 1;  %[1] conversion from czi to tif
do_ilastik = 0;     %[2] segmentation with ilastik
do_spot_seg= 0;    %[3] post processing for final DAB spot identification
do_split = 0;       %[4] split images into top/bottom to separate hemispheres
do_split_manual = 0; %[4b] manually adjust the splitting mask
do_area_analysis = 0; %[5] spot area analysis

% Choose specific files to process according to their index in the 
% directory. Leave empty as [] to process all files. Note that if you are
% running area analysis, ALL available files will be processed regardless 
% of what is selected below. 
specific_files_to_process = []; %e.g. [1,2,3,5:10]

% Force procedure even if the corresponding output file already exists,
% e.g. if something must be redone. Set to 1 to force procedure, 0 if
% files that have already been processed can be skipped. Note that manual 
% image splitting is always forced if specified.
force_process = 0; 



%%%%%%%%%%%%%%%%% SEGMENTATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Only modify if you are changing labels in the ilastik project
BKGRD_label = 1;   %label in the ilastik project for background
DAB_label = 5;     %label corresponding to (good) DAB staining

%%%% Parameters for post-processing: only modify if scale or change in acquisition settings demands it
dilation_r = 35;   %disk radius for background dilation to remove tissue edge labelling for spot identifn
area_thresh = 50000; %minimum area threshold to remove small artifacts around main tissue sample



%% %%%%%%%%%%%%%%%%%%% ----- DO NOT MODIFY BELOW THIS LINE -----%%%%%%
%Add external functions to the path
addpath(genpath(dir_bfmatlab));

%Initialize dataset object
dataset = DAB_spot_detector('dir_src',dir_src,'dir_tif',dir_tif,'dir_ilastik_out',dir_ilastik_out, ...
    'dir_final_seg_out',dir_final_seg_out,'ilastik_command_headless',ilastik_command_headless,...
    'ilastik_project_1',ilastik_project_1,'BKGRD_label',BKGRD_label,'DAB_label',DAB_label,...
    'dilation_r',dilation_r,'area_thresh',area_thresh);

%Create directories for processing if they do not exist
dataset.create_dirs;

%Parse to determine the files to be processed
if isempty(specific_files_to_process)
    files_to_process = @(x) 1:numel(x);
else
    files_to_process = @(x) reshape(specific_files_to_process,1,numel(specific_files_to_process));
end

%% [1] convert czi images to tif

if do_conversion
    
fprintf('Converting CZI in: %s\n',dataset.dir_src); %#ok<UNRCH>
disp('-----------');

fnames = dir([dir_src filesep '*.czi']);
f2p = files_to_process(fnames);
did_conversion = false(numel(f2p),1); %Flag to keep track of failed conversions

parfor ix = 1:numel(f2p)
    fprintf('Processing file %d of %d ...\n',ix,numel(f2p))
    try
        did_conversion(ix) = dataset.convert_file(fnames(f2p(ix)).name,force_process);
    catch
        fprintf('\tError encountered during conversion of file %d: %s\n',f2p(ix),fullfile(dataset.dir_src,fnames(ix).name));
    end
end
fprintf('Finished conversion from CZI to TIF.\n')

end


%% [2] perform ilastik segmentation on all converted images in folder

if do_ilastik

%Check for Zen CZI to TIF conversion 
D = dir(dataset.dir_tif); %#ok<UNRCH>
zen_used = all([D.isdir]);
if zen_used
    fnames = dataset.gen_zen_filenames(dataset.dir_tif);
else
    fnames = dir([dir_tif filesep '*.tif']);
end

f2p = files_to_process(fnames);
did_ilastik = false(numel(f2p),1); %Flag to track successful ilastik segmentation

fprintf('Performing ilastik-based segmentation\n'); 
fprintf('ilastik project file: %s',ilastik_project_1);
disp('------------');

parfor ix = 1:numel(f2p)    
    fprintf('Processing file %d of %d ...\n',ix,numel(f2p))
    try
        did_ilastik(ix) = dataset.run_ilastik(fnames(f2p(ix)).name,force_process); %#ok<PFBNS>
    catch
        fprintf('\tFile Error encountered during ilastik segmentation of file %d: %s\n',f2p(ix),fnames(f2p(ix)).name);
    end
end
fprintf('Finished ilastik-based segmentation.\n')

end

%% [3] post-segmentation processing on all ilastik-segmented images

if do_spot_seg
    
fnames = dir(fullfile(dataset.dir_ilastik_out,'*.tif')); %#ok<UNRCH>
f2p = files_to_process(fnames);
did_spot_seg = false(numel(f2p),1); 

fprintf('Isolating segmented plaques in folder: %s\n',dataset.dir_ilastik_out);
disp('------------');

parfor ix = 1:numel(f2p)
    fprintf('Processing file %d of %d ...\n',ix,numel(f2p))
    try
        did_spot_seg(ix) = dataset.run_spot_segmentation(fnames(f2p(ix)).name,force_process) %#ok<PFBNS>
    catch
        fprintf('\Error encountered during spot segmentation of file %d: %s\n',f2p(ix),fnames(f2p(ix)).name);
    end
end
fprintf('Finished post-processing plaque segmentation.\n');

end


%% [4] split images along symmetry line to separate L/R hemispheres

if do_split

dataset.clean_tif_filenames(dataset.dir_final_seg_out); %#ok<UNRCH>
fnames = dir(fullfile(dataset.dir_ilastik_out,'*.tif')); 
f2p = files_to_process(fnames);
did_split = false(numel(f2p),1);

fprintf('Calculating mask to separate hemispheres in folder: %s\n',dataset.dir_ilastik_out);
disp('------------');

for ix = 1:numel(f2p)
    fprintf('Processing file %d of %d ...\n',ix,numel(f2p))
    try
        did_split(ix) = dataset.run_split(fnames(f2p(ix)).name,force_process);
    catch
        fprintf('\tFile Error encountered during hemisphere splitting in file %d: %s\n',f2p(ix),fnames(f2p(ix)).name);
    end
end

end

%% [4b] split images manually to separate L/R hemispheres

if do_split_manual

dataset.clean_tif_filenames(dataset.dir_final_seg_out); %#ok<UNRCH>    
fnames = dir(fullfile(dataset.dir_ilastik_out,'*.tif')); 
f2p = files_to_process(fnames);
did_split_manual = false(numel(f2p),1);

fprintf('Drawing mask to separate hemispheres in folder: %s\n',dataset.dir_ilastik_out);
disp('------------');

for ix = 1:numel(f2p)
    fprintf('Processing file %d of %d ...\n',ix,numel(f2p))
    try
        did_split_manual(ix) = dataset.run_split_manual(fnames(f2p(ix)).name);
    catch
        fprintf('\tFile Error encountered during hemisphere splitting in file %d: %s\n',f2p(ix),fnames(f2p(ix)).name);
    end
end

end
%% [5] analysis: spot area/total segmented area using the ilastik labeled images

if do_area_analysis

dataset.clean_tif_filenames(dataset.dir_final_seg_out); %#ok<UNRCH>
fnames = dir([dir_ilastik_out filesep '*.tif']); 
C = cell(numel(fnames),5,2); %cell to store analysis results
C_header = cell(5,1);
did_area_analysis = false(numel(fnames),1);

for ix = 1:numel(fnames)
    fprintf('Analyzing image %d of %d\n',ix,numel(fnames))
    try
        [Cix,C_header,f] = dataset.run_area_analysis(fnames(ix).name);
        C(ix,:,1) = Cix(1,:);
        C(ix,:,2) = Cix(2,:);
        did_area_analysis(ix) = f;
    catch
        fprintf('\tError encountered during area analysis of file %d: %s\n',ix,fnames(ix).name);
    end
end

C = [C(:,:,1); C(:,:,2)];
is_any_C = any(cellfun(@(x) ~isempty(x),C),2);

C = [C_header; C(is_any_C,:)];
writecell(C,fullfile(dataset.dir_final_seg_out,'area_analysis.csv'))

end
