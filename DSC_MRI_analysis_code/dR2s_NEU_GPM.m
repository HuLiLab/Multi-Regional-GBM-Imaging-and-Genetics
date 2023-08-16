% This code computes the mean dR2s for the perfusion data (glycolytic 
% and neuronal) on a 10x10 ROI mask and also the WBNE dR2s curves.

fprintf('-------------------------------------------------------------------------------------\n');
t = datetime('now');
display(['Code started at ' num2str(t.Hour) ':' num2str(t.Minute, '%02i')]);

%% DSC MRI data

directory = '/Users/stokeslab/Desktop/BNI/Leland/Code/DSC_MRI_analysis_code/Sample_data';
cd(directory)

%get into the patient data list
patients = dir('MCH*');

for p = 1:size(patients,1)
    %% Data intake
    patientFolder = fullfile(directory,patients(p).name);         
    cd(patientFolder);
    
    folder = dir(fullfile(patientFolder,'MCH_*'));
    files = {folder.name};
    
    % Read data
    epi = niftiread(fullfile(patientFolder,files{2}));
    for i = 1:size(epi,4)
        w = squeeze(epi(:,:,:,i));
        perf(:,:,:,i) = imrotate3(w,90,[0 0 1]);
    end

    clear epi w i
    
    % extract TR and TE
    fname = jsondecode(fileread(files{1}));   
    TR = fname.RepetitionTime;
    TE = fname.EchoTime;
    disp(files{1});
    
    clear files patientFolder fname

    %brain mask - skull stripping
    for i = 1:size(perf,3)
    
        J = cast(perf(:,:,i,1),'double');
        a = squeeze(J(:,:,1));

        normMriImage1 = a./max(a(:));

        % Creating binary mask
        normBW = im2bw(normMriImage1,0.15);

        % clean mask of small areas
        cleanBW = bwareaopen(normBW,500);

        % Label Regions
        labelBW = bwlabel(cleanBW);
        % Calculate number of voxels in each region
        regStats = regionprops(cleanBW,'area');
        allAreas = [regStats.Area];

        % Find region with largest area
        [brainArea brainInd] = max(allAreas);

        % Extract the largest region using ismember()
        brainRegion = ismember(labelBW, brainInd);

        % Convert from integer labeled image into binary image.
        brainBW = brainRegion > 0;

        % Fill holes in brain mask
        brain_mask(:,:,i) = imfill(brainBW,'holes');

    % figure();imagesc(double(brain_mask).*normMriImage1);colormap('gray');
    % figure();imagesc(a);colormap('gray')
    % 
        clear allAreas brainArea brainBW brainInd brainRegion
        clear cleanBW labelBW normBW regStats a
        clear I normMriImage1
    end

    perf = cast(perf,'double');
    vol = perf.*brain_mask;

    %% determine steady state, gradient and peak timepoints

    [nx,ny,nz,nt] = size(vol);
    mean_tc = mean(reshape(vol, [nx*ny*nz nt])); %mean time course
    slope = floor(diff(mean_tc)); %find slope of curve to first determine steady state location
    % %         if slope(1) < -20 %AMS commented out 7/21/2016
    % %             ss_tp = find(slope >= -1, 1, 'first') + round(3/TR); %steady state reached, add 3 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1 %AMS commented out 7/21/2016
    if slope(1) < -20 || slope(1) > 20 %include positive slope - important for SAGE
        ss_tp = find(slope >= -1 & slope <= 1, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
        if ss_tp > nt/2
            warning('Did not find steady-state in first half of time-course, trying again')
            ss_tp = find(slope >= -20 & slope <= 20, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
        end
    else
        ss_tp = 1; %assume no dummy scans included in data
    end
    if isempty(ss_tp) %catch empty ss_tp
        ss_tp = 1;
    end
    [pks,locs,w] = findpeaks(-mean_tc(:,ss_tp:end));
    [~, gd_index] = max(pks);
    gd_tp = locs(gd_index) + ss_tp + 1 - round(3/TR) - floor(w(gd_index)*1.5); %contrast arrival
    pk_tp = locs(gd_index) + ss_tp - round(3/TR) + 2; %time to peak

%     h(3) = figure('Color', [1 1 1]); 
%     title('Mean Time Course Across Volume'); hold on;
%     plot(mean_tc, 'k', 'LineWidth', 2'); hold on;

    clear h locs mean_tc nx ny nz ne nt pks w gd_index slope
    
    %% extract WBNE voxels - peak tp > -6SD & postcontrast tp +/- 2SD

    % vol = cast(perf(:,:,z,:),'double');
    [nx,ny,nz,nt] = size(vol);
    ne = 1;
    volVec1 = reshape(vol, [nx*ny*nz ne nt]);

    mnSI_end = mean(volVec1(:,1,end-10:end),3,"omitnan");

    pos_threshold = mean(volVec1(:,1,ss_tp:gd_tp),3,"omitnan") + (2*std(squeeze(volVec1(:,1,ss_tp:gd_tp)), [], 2,"omitnan"));
    neg_threshold = mean(volVec1(:,1,ss_tp:gd_tp),3,"omitnan") - (2*std(squeeze(volVec1(:,1,ss_tp:gd_tp)), [], 2,"omitnan"));    
    peak_threshold = mean(volVec1(:,1,ss_tp:gd_tp),3,"omitnan") - (6*std(squeeze(volVec1(:,1,ss_tp:gd_tp)), [], 2,"omitnan")); 
    min_value = squeeze(volVec1(:,1,pk_tp));

    nonenhance_map1 = zeros(nx,ny,nz);
    nonenhance_map2 = zeros(nx,ny,nz);
    nonenhance_map3 = zeros(nx,ny,nz);

    nonenhance_map1(mnSI_end < pos_threshold) = 1;
    nonenhance_map2(mnSI_end > neg_threshold) = 1;
    nonenhance_map3(min_value < peak_threshold) = 1;

    nonenhance_map = nonenhance_map1 + nonenhance_map2 + nonenhance_map3;
    nonenhance_map(nonenhance_map == 1) = 0;
    nonenhance_map(nonenhance_map == 2) = 0; 
    nonenhance_map(nonenhance_map == 3) = 1;


    clear pos_threshold neg_threshold mnSI_end peak_threshold min_value
    clear nonenhance_map1 nonenhance_map2 nonenhance_map3

    S0_TE = squeeze(mean(volVec1(:,ss_tp:gd_tp),2,"omitnan")); S0_TE = repmat(S0_TE, [1 nt]);

    % delta R2 star
    dR2_UC = -(1/(TE)).*(log(squeeze(volVec1(:,:))./S0_TE));
    dR2_UC(isinf(dR2_UC)) = 0;
    dR2_UC(isnan(dR2_UC)) = 0; 

    % Average delta R2 star based on non-enhancing mask
    dR2s_WBNE = dR2_UC .* repmat(nonenhance_map(:),[1 nt]);
    dR2s_WBNE(dR2s_WBNE == 0) = NaN;
    dR2s_WBNE = mean(dR2s_WBNE,"omitnan");
    baseline_dR2sWBNE = mean(dR2s_WBNE(ss_tp:gd_tp),"omitnan");
    dR2s_WBNE(1:ss_tp) = repmat(baseline_dR2sWBNE, [1 ss_tp]);

    clear baseline_dR2sWBNE J dR2_UC nonenhance_map S0_TE volVec1 nx ny nz ne nt brain_mask

    %% Tumor ROI Coordinates

    sheet = readtable('/Users/stokeslab/Desktop/BNI/Leland/Code/DSC_MRI_analysis_code/ROI_Coordinates.xlsx');
    
    sub_name = folder.name; 
    name = sub_name(1:8);

    [~,ind] = ismember(name,sheet.Subject);

    x = sheet.X(ind,1);
    y = sheet.Y(ind,1);
    z = sheet.Z(ind,1);
    
    fprintf("X: %d,Y: %d,Z: %d\n",x,y,z)

    a = (x-5);
    b = (y-5);

    roi = squeeze(perf(a:a+9,b:b+9,z,:));   %10x10 ROI

    % delta R2 star ROI
    [nx,ny,nt] = size(roi); 
    roi = cast(roi,'double');

    volVec1 = reshape(roi, [nx*ny nt]);

    S0_TE1 = squeeze(mean(volVec1(:,ss_tp:gd_tp),2,"omitnan")); S0_TE1 = repmat(S0_TE1, [1 nt]);

    dR2s_UC = -(1/(TE)).*(log(squeeze(volVec1(:,:))./S0_TE1));
    dR2s_UC(isinf(dR2s_UC)) = 0;
    dR2s_UC(isnan(dR2s_UC)) = 0;    

    delta_R2 = mean(dR2s_UC);

    clear x y z a b roi m nx ny nz volVec1 S0_TE1 dR2s_UC
    
    WBNE(p,:) = dR2s_WBNE;
    ROI(p,:) = delta_R2;
    
    clear dR2s_WBNE delta_R2 ss_tp gd_tp pk_tp TE TR i nt vol perf
    
    
end

%% Plotting signal for a single subject

figure()
plot(ROI,'Color','r','LineWidth', 3)
set(gca,'FontSize',17,'FontWeight','bold')
hold on
plot(WBNE,'Color','k','LineWidth', 3)
xlabel('Time (sec)','FontSize',20,'FontWeight','bold')
ylabel('Delta R2* (1/sec)','FontSize',20,'FontWeight','bold')
legend('ROI', 'WBNE','FontSize',12)


%% Plotting mean signals of all subjects

% avg_WBNE = mean(WBNE);
% avg_del = mean(ROI);
% 
% figure()
% plot(avg_del,'Color','r','LineWidth', 3)
% set(gca,'FontSize',17,'FontWeight','bold')
% hold on
% plot(avg_WBNE,'Color','k','LineWidth', 3)
% % title('Glycolytic','FontSize',18)
% xlabel('Time (sec)','FontSize',20,'FontWeight','bold')
% ylabel('Delta R2* (1/sec)','FontSize',20,'FontWeight','bold')
% legend('ROI', 'WBNE','FontSize',12)
% 
% 
t = datetime('now');
display(['Code ended at ' num2str(t.Hour) ':' num2str(t.Minute, '%02i')]);
fprintf('-------------------------------------------------------------------------------------\n');
