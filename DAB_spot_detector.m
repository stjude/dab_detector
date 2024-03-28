classdef DAB_spot_detector
    %Class containing main functions for d'Azzo group DAB spot detection.
    %--------------------------------------------------------------------
    %
    %Authors: Mia Panlilio and Khaled Khairy (2021)
    %
    %2021-10-21 MP. Changed TIF write method from imwrite to Tiff class
    %               write to accommodate larger files using BigTIFF format.
    %               Note that imwrite is limited to (2^32 - 1) bytes. 
    %2021-09-27 MP. Fixed ROI-finding bug for conversion from CZI to TIF.
    %
    %--------------------------------------------------------------------
    
    properties (Access = 'public')
        dir_src = '';
        dir_tif = '';
        dir_ilastik_out = '';
        dir_final_seg_out = '';
        ilastik_command_headless = 'LAZYFLOW_THREADS=1 LAZYFLOW_TOTAL_RAM_MB=20000 <path-to-run_ilastik.sh> --headless --readonly READONLY'; % where is ilastik?
        ilastik_project_1 = '<ILP-classifier-location>'; % pretrained ilastik project: at same scale and for isotropic voxels   
        BKGRD_label = 1;
        DAB_label = 5;
        dilation_r = 35;
        area_thresh = 50000;    
        use_legacy2tiff = false;    
    end
    
    properties (Access = 'private')
        dir_split_mask = '';
        dir_split_mask_vis = '';
        dir_total_seg_area = '';
    end
    
    methods (Static = false)
        
        function self = DAB_spot_detector(varargin)
        %DAB_SPOT_DETECTOR('Key','Value')
        %Class containing key functions for DAB spot detection and area
        %analysis. 
        %
        %Key-value input combinations to process a specific dataset
        %
        % 'dir_src' : directory containing source files (CZI)
        %
        % 'dir_tif' : directory for TIFs converted from source
        %
        % 'dir_ilastik_out' : directory for Ilastik segmentation output
        %
        % 'dir_final_seg_out' : directory to contain final segmentation
        %                       (spot and whole tissue sample)
        %
        % 'BKGRD_label' : label in segmented images corresponding to the
        %                 background. Default is 1 (from Ilastik).
        %
        % 'DAB_label' : label in segmented images corresponding to the DAB
        %               plaques. Default is 5.
        %
        % 'dilation_r' : radius of disk used to dilate background tissue
        %                mask to remove edge pixels incorrectly categorized
        %                as plaques. Default is 35. 
        %
        % 'area_thresh' : minimum object size (in pixels) below which
        %                 segmented objects are removed. Default is 50000.
        %
        %Key-value input combinations for Ilastik parameterization
        %
        % 'ilastik_command_headless' : system command to run Ilastik
        % (default is 'LAZYFLOW_THREADS=1 LAZYFLOW_TOTAL_RAM_MB=20000
        % /research/sharedresources/cbi/data_exchange/dazzogrp/Jayce_Plaque_Quantification/Khaled/tools/ilastik-1.3.2-Linux/run_ilastik.sh --headless --readonly READONLY')
        %
        %  'ilastik_project_1' : pretrained Ilastik project (default is
        %  '/research/sharedresources/cbi/data_exchange/dazzogrp/Jayce_Plaque_Quantification/Khaled/initial_classifier.ilp')
        %
    
        %Parse inputs
        pp = inputParser;
        props = properties(self);
        for k = 1:numel(props)
            pp.addParameter(props{k},self.(props{k}));
        end
        pp.parse(varargin{:})
        
        for k = 1:numel(props)
            self.(props{k}) = pp.Results.(props{k});
        end
        
        self.dir_split_mask = fullfile(self.dir_final_seg_out,'split_mask');
        self.dir_split_mask_vis = fullfile(self.dir_final_seg_out,'split_mask_vis');
        self.dir_total_seg_area = fullfile(self.dir_final_seg_out,'total_seg_area');
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function czi2tif(self,czi_src,dir_dest)
        %CZI2TIF(SRC_FILE,DEST_FILE)
        %Converts CZI named by src_file to TIF in the destination directory
        %dir_dest.
        if ~self.use_legacy2tiff
            bf2tiff(czi_src,dir_dest)
        else
            self.czi2tif_legacy(czi_src,dir_dest)
        end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = czi2tif_legacy(self,czi_src,dir_dest)
        %out = CZI2TIF_legacy(SRC_FILE,DEST_FILE)
        %Converts the CZI indicated by src_file to TIF in the destination
        %directory dir_dest. If requested, output out is the bfopen-ed CZI 
        %file.
        %2021-10-21 MP. Updated to use Tiff class write method to BigTiff
        %               due to increases in image size.
        %
        %2023-04-10 MP. DEPRECATED. Class now defaults to use singularity
        %               container with python-bioformats, using wrapper
        %               bf2tiff.m

        [~,fname,~] = fileparts(czi_src);
        im = bfopen(czi_src);
        omeMeta = im{1,4};
        
        roi_ix = self.get_unique_roi(omeMeta);
        for k = 1:numel(roi_ix)
            im1 = im{roi_ix(k)};
            imR = im1{1,1};
            imG = im1{2,1};
            imB = im1{3,1};
            im_out = zeros(size(imR,1), size(imR,2), 3);
            im_out(:,:,1) = imR;
            im_out(:,:,2) = imG;
            im_out(:,:,3) = imB;
            
            %Write array to BigTIFF:
            im_out = im2uint8(mat2gray(im_out));
            fn_out = fullfile(dir_dest,sprintf('%s_r%d.tif',fname,k));
            
            %Use Tiff class to write to BigTIFF format
            t = Tiff(fn_out,'w8'); 
            t.setTag('ImageLength',size(im_out,1)); 
            t.setTag('ImageWidth', size(im_out,2));
            t.setTag('Photometric', Tiff.Photometric.RGB);
            t.setTag('BitsPerSample', 8); 
            t.setTag('SamplesPerPixel', size(im_out,3)); 
            t.setTag('TileWidth', 128); 
            t.setTag('TileLength', 128);
            t.setTag('Compression', Tiff.Compression.LZW); 
            t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky); 
            t.setTag('Software', 'MATLAB'); 
            t.write(im_out); 
            t.close(); 
        end
        
        if nargout==1
            out = im;
        end
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function create_dirs(self)
        %CREATE_DIRS(SELF)
        %Creates, if they do not exist, all directories for processed
        %output
        
        dstr = {self.dir_tif, self.dir_ilastik_out, self.dir_final_seg_out, ...
            self.dir_split_mask, self.dir_split_mask_vis,self.dir_total_seg_area};

        for k = 1:numel(dstr)
            if ~isfolder(dstr{k})
                mkdir(dstr{k})
            end
        end
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = convert_file(self,fname_czi,override)
        %f = CONVERT_FILE(SELF,FNAME_CZI,OVERRIDE).
        %Converts the CZI referencecd in fname to TIF. If the converted TIF
        %already exists, conversion will not be attempted unless override
        %is set to 1 (true)
        %
        % fname_czi : must be relative to dir_src
        %
        % override : logical flag to convert even if the file already
        %            exists. Default is false.
        
        if nargin<3; override = false; end
        
        [~,fname,~] = fileparts(fullfile(self.dir_src,fname_czi));
        
        already_converted = exist(fullfile(self.dir_tif,[fname,'.tif']),'file') || exist(fullfile(self.dir_tif,[fname,'_r1.tif']),'file');
        if override || ~already_converted
            %Convert CZI to TIF if it has not already been done or if user
            %is not forcing the conversion
            fprintf('\tConverting CZI to TIF: %s.\n',fullfile(self.dir_src,fname_czi));
            self.czi2tif(fullfile(self.dir_src,fname_czi),self.dir_tif)
        else
            fprintf('\tAlready converted: %s.\n',fullfile(self.dir_src,fname_czi));
        end
        f = true; 
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = run_ilastik(self,fname_tif,override)
        %f = RUN_ILASTIK(SELF,FNAME_TIF,OVERRIDE)
        %Performs segmentation using the pretrained Ilastik network
        %indicated by the class property ilastik_project_1 and the
        %parameters defined by the command in ilastik_command_headless.
        %
        % fname_tif : must be relative to dir_tif
        %
        % override : logical flag to force segmentation even if a 
        %            segmented file already exists. Default is false.
        
        if nargin<3; override = false; end
        
        [~,base_fname,~] = fileparts(fname_tif); %in case TIFs are in subdirectory from Zen
        base_fname = [base_fname '.tif'];
        
        already_segmented = exist(fullfile(self.dir_ilastik_out,base_fname),'file');
        if override || ~already_segmented
            %Run Ilastik segmentation if it has not already been done or if
            %the user is not forcing segmentation
            fprintf('\tRunning ilastik segmentation for: %s\n',fname_tif)
            fn_in = fullfile(self.dir_tif,fname_tif);
            fn_out = fullfile(self.dir_ilastik_out,base_fname);
            self.segment_with_ilastik(fn_in,fn_out);
        else
            fprintf('\tFile has already been segmented: %s\n',fname_tif)
        end
        f = true;
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function segment_with_ilastik(self,fname_in,fname_out)
        %SEGMENT_WITH_ILASTIK(SELF,FNAME_IN,FNAME_OUT)
        %Using the parameters set by the class properties ilastik_project_1
        %and ilastik_command_headless, runs ilastik headlessly to segment a
        %the file given by its full path name, fname.
        
        command = [self.ilastik_command_headless ...
            ' --project=' self.ilastik_project_1 ...
            ' --export_source="Simple Segmentation" --output_format="tif" --output_filename_format=' ...
            fname_out ' ' fname_in];
        system(command);
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = run_spot_segmentation(self,fname_tif,override)
        %f = RUN_SPOT_SEGMENTATION(SELF,FNAME_TIF,OVERRIDE)
        %Performs post-processing on segmented image to isolate DAB spots.
        
        if nargin<3; override = false; end
        
        already_segmented = exist(fullfile(self.dir_final_seg_out,fname_tif),'file');
        if override || ~already_segmented
            fprintf('\tPost-processing to isolate spots: %s.\n',fname_tif);
            fn_in = fullfile(self.dir_ilastik_out,fname_tif);
            self.segment_spots(fn_in,self.dir_final_seg_out,self.DAB_label,...
                self.BKGRD_label,self.dilation_r)
        else
            fprintf('\tFile has already been spot segmented: %s.\n',fname_tif)
        end
        f = true;
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = run_split(self,fname_tif,override)
        %f = RUN_SPLIT(SELF,FNAME_TIF,OVERRIDE)
        %Calculates and saves the mask separating brain hemispheres using
        %split_image.m. 
        
        if nargin<3; override = false; end
        
        already_split = exist(fullfile(self.dir_split_mask,fname_tif),'file');
        if override || ~already_split
            fprintf('\tCalculating split mask: %s.\n',fname_tif);
            fn_in = fullfile(self.dir_ilastik_out,fname_tif);
            fn_out = fullfile(self.dir_split_mask,fname_tif);
            
            B_ = self.rm_tissue_background(fn_in,self.BKGRD_label,self.area_thresh);
            B = double(B_);
            B(B~=0) = 255;
            M = split_image(B);
            imwrite(uint8(mat2gray(M)*255),fn_out)
            
            %Save a copy of the segmented image with the split mask overlay
            fn_out_vis = fullfile(self.dir_split_mask_vis,fname_tif);
            imwrite(uint8(mat2gray(M+double(B_))*255),fn_out_vis)
        else
            fprintf('\tFile has already been split: %s.\n',fname_tif)
        end
        f = true;
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f = run_split_manual(self,fname_tif)
        %f = RUN_SPLIT_MANUAL(SELF,FNAME_TIF)
        %Presents the user with the image in fname_tif with the option to
        %draw a polygon indicating the splitting mask to separate different
        %brain hemispheres.
        
        f = false;
        I = imread(fullfile(self.dir_tif,fname_tif));
        
        hf = figure('name','Manual splitting mask selection');
        hp1 = uipanel(hf,'position',[0,0,0.1,1]);
        hp2 = uipanel(hf,'position',[0.1,0,0.9,1]);
        uicontrol(hp1,'string','SAVE','units','normalized','position',[0.05,0.85,0.9,0.1],...
            'fontsize',12,'callback',@save_new_mask);
        ha = axes(hp2,'position',[0.05,0.05,0.9,0.9]);
        
        imshow(I,'parent',ha);
        roi = drawpolygon(ha);
        
            function save_new_mask(~,~)
            %Callback function for save button to save current polygon to
            %a binary mask file.
            roi = roi.Position;
            M = poly2mask(roi(:,1),roi(:,2),size(I,1),size(I,2));
            if M(1)~=1
                %Invert so that the mask is true for the top hemisphere
                M = M~=1;
            end
            fn_out = fullfile(self.dir_split_mask,fname_tif);
            imwrite(uint8(mat2gray(M)*255),fn_out);
            
            %Save new mask-segmentation overlay image for visualization
            fn_in = fullfile(self.dir_ilastik_out,fname_tif);
            B_ = self.rm_tissue_background(fn_in,self.BKGRD_label,self.area_thresh);
            B = double(B_);
            B(B~=0) = 255;
            fn_out_vis = fullfile(self.dir_split_mask_vis,fname_tif);
            imwrite(uint8(mat2gray(M+double(B_))*255),fn_out_vis)
            
            fprintf('\tManually selected mask saved: %s.\n',fname_tif);
            f = true;
            end
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [C,header,f] = run_area_analysis(self,fname_tif)
        %[C,header,f] = RUN_AREA_ANALYSIS(SELF,FNAME_TIF)
        
        %Check for hemisphere splitting masks
        is_split_mask = true;
        if ~isfolder(self.dir_split_mask)
            is_split_mask = false;
        else
            D = dir(fullfile(self.dir_split_mask,'*.tif'));
            if isempty(D)
                is_split_mask = false;
            end
        end
        
        fprintf('\tAnalyzing spot and segmented area: %s.\n',fname_tif);
    
        B = self.rm_tissue_background(fullfile(self.dir_ilastik_out,fname_tif),self.BKGRD_label,self.area_thresh); %whole sample/tissue segmentation
        S = imread(fullfile(self.dir_final_seg_out,fname_tif));
        S = S~=0;
        
        if is_split_mask
            M = imread(fullfile(self.dir_split_mask,fname_tif)); %split mask
            M = M~=0;
            split_name = {'top_','bottom_'};
            C = cell(2,5);
            for jx = 0:1
                M_ = M==jx;
                r_B = regionprops(B&M_,'Area');
                r_S = regionprops(S&M_,'Area');
                a_B = [r_B.Area];
                a_S = [r_S.Area];
                
                C(jx+1,:) = [[split_name{jx+1} fname_tif], {numel(a_S),sum(a_S),sum(a_B),sum(a_S)/sum(a_B)}];
            end
        else
            r_B = regionprops(B,'Area'); %calculate total segmented area
            r_S = regionprops(S,'Area'); %calculate area for each DAB spot

            a_B = [r_B.Area]; %vector of all segmented tissue areas
            a_S = [r_S.Area]; %vector of all spot areas 
            
            C = [fname_tif, {numel(a_S),sum(a_S),sum(a_B),sum(a_S)/sum(a_B)}]; %store results
        end
        
        %Save total segmented area to file
        fn_out = fullfile(self.dir_total_seg_area,fname_tif);
        imwrite(uint8(mat2gray(B)*255), fn_out);
        
        %Column titles corresponding to C
        header = {'filename','n_spots','total_spot_area','total_seg_area','area_ratio'};
        f = true;
        
        end
        
    end
    
    
    %%
    methods (Static = true)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [n_spots,area_spots,area_total,a_S,a_B] = area_analysis(S,B)
        %[n_spots,area_spots,area_total] = area_analysis(S,B)
        %From the spot segmented image S and segmented tissue sample B,
        %calculates the number of spots, the area that they
        %comprise, and the total area segmented.
        %
        %Additional outputs
        %
        % a_S : vector containing individual areas of segmented spots
        %
        % a_B : vector containing individual areas of segmented tissue
        %       sample, e.g. in case of hemisphere/lobe separation
        
        r_B = regionprops(B,'Area'); %calculate total segmented area
        r_S = regionprops(S,'Area'); %calculate area for each DAB spot

        a_B = [r_B.Area]; %vector of all segmented tissue areas
        a_S = [r_S.Area]; %vector of all spot areas
        
        n_spots = numel(a_S);
        area_spots = sum(a_S);
        area_total = sum(a_B);
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function clean_tif_filenames(dirname)
        %f = CLEAN_TIF_FILENAMES(dirname)
        %If TIF files in the directory dirname have extra '.tif' appended
        %to the end, remove the additional file extension.
        
        fnames = dir(fullfile(dirname,'*.tif.tif'));
        fnames = {fnames(:).name};
        if ~isempty(fnames)
            fprintf('Cleaning up TIF filenames in: %s\n',dirname)
            fprintf('--------------------\n')
            for k = 1:numel(fnames)
                movefile(fullfile(dirname,fnames{k}), fullfile(dirname,fnames{k}(1:end-4)));
                if mod(k,25) || k==numel(fnames)
                    fprintf('\t%d of %d files renamed.\n',k,numel(fnames))
                end
            end
            fprintf('\n')
        end
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fnames = gen_zen_filenames(dir_src)
        %FNAMES = GEN_ZEN_FILENAMES(DIR_SRC)
        %Generates the relative path list of images converted by Zen
        %contained in the folder dir_src.
        % fnames : struct with field 'name'
        
        D = dir(dir_src);    
        fnames = struct;
        fnames.name = '';
        for ix = 3:numel(D)
            d = dir(fullfile(dir_src,D(ix).name,'*.tif'));
            fn = strcat([D(ix).name filesep],{d.name});
            fn = cell2struct(fn,'name',1);
            fnames = [fnames; fn]; %#ok<*AGROW>
        end
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function idx = get_unique_roi(omepyr)
        %idx = GET_UNIQUE_ROI(OMEPYR)
        %Retrieves a list of indices for unique ROIs within a given CZI at the 
        %highest resolution.
        %
        % omepyr : OMEPyramid metadata, i.e. any cell in the fourth column of a
        %           bfopen read-in.
        %
        %2021-03-26 Mia Panlilio. Solution to retrieving multiple ROI from CZI 
        %series, omitting lower resolution duplicates that might also be stored by
        %looking up (x,y,z) stage motor positions.
        %
        %2021-09-27 Update by MP. Fixed unique ROI finding. Encountered dataset
        %where NaNs were yielded only for z-stage at lower magnification level, 
        %resulting in artificially unique positions.

        %Retrieve (x,y,z) stage positions and image area in pixels
        n = omepyr.getImageCount();
        pos = nan(n,3);
        area = nan(n,1);
        roi_found = true(n,1);

        for k = 1:n
            try
                pos(k,:) = [stagepos(omepyr.getStageLabelX(k-1)), ...
                    stagepos(omepyr.getStageLabelY(k-1)), ...
                    stagepos(omepyr.getStageLabelZ(k-1))];

                area(k) = omepyr.getPixelsSizeX(k-1).getValue()*...
                    omepyr.getPixelsSizeY(k-1).getValue()*...
                    omepyr.getPixelsSizeZ(k-1).getValue();
            catch
                roi_found(k) = false;
            end
        end
        
        %Ensure that a stage position was retrieved for all dimensions
        missing_dim = any(isnan(pos),2);
        roi_found(missing_dim) = false;
        
        %Rearrange indices to ensure highest resolution images are found first
        [~,sortID] = sort(area,'descend');
        pos = pos(sortID);
        roi_found = roi_found(sortID);

        %Find set of unique stage positions
        [~,uniqueID] = unique(pos,'rows');
        idx = sortID(uniqueID);

        %Ensure that an ROI was actually obtained
        roi_found = roi_found(uniqueID);

        idx = idx(roi_found);
        idx = sort(idx,'ascend');
        
        %--------ADDITIONAL FUNCTIONS FOR ROI DETECTION--------------------

        function val = stagepos(str)
            r = regexp(char(str),'value\[[\-]*(\d*\.\d*)\]','tokens');
            if ~isempty(r)
                val = str2double(r{1});
            else
                val = nan;
            end
        end

        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function B = rm_tissue_background(tif_src,BKGRD_label,area_thresh)
        %B = RM_TISSUE_BACKGROUND(tif_src,BKGRD_label,area_thresh)
        %On the Ilastik labelled image tif_src, performs cleanup to isolate 
        %the tissue sample from background. Returns binary image where 1 is
        %sample and 0 the background.
        % BKGRD_label : binary label for the image background (1 for Ilastik)
        %
        % area_thresh : lower object area threshold to clean up background
        %               using bwareaopen
        
        B = imread(tif_src);
        B = B~=BKGRD_label; %isolate tissue from background
        B = bwareaopen(B,area_thresh,4); %remove small objects
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function segment_spots(tif_src,dir_dest,DAB_label,BKGRD_label,dilation_r)
        %SEGMENT_SPOTS(TIF_SRC,DIR_DEST,DAB_label,BKGRD_label,dilation_r)
        %Further segments Ilastik labelled image in tif_src to identify
        %plaques. Resulting image is saved to dir_dest folder.
        
        im = imread(tif_src);
        I = im;
        I(I==DAB_label) = 255;  % get DAB_label
        I(I<128) = 0;
        
        imb = im; % get BKGRD_label
        imb(imb==BKGRD_label) = 255;
        imb(imb<128) = 0;

        if dilation_r==0
            imd = imdilate(imb, strel('disk', dilation_r));
        else
            imd = imb
        end
        I(imd==255) = 0;
        
        [~,fname,~] = fileparts(tif_src);
        fn_out = fullfile(dir_dest,[fname,'.tif']);
        imwrite(uint8(mat2gray(I)*255), fn_out);
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

end
