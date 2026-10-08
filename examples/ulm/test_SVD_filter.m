%% Set paths
local_path = fullfile(ustb_path(), 'data');

input_folder = fullfile(local_path, 'RF_channeldata');
output_folder = fullfile(local_path, 'beamformed_slidingSVD');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

files = dir(fullfile(input_folder, '*.uff'));

[~, sortOrder] = sort({files.name});
files = files(sortOrder);

uffFiles = fullfile({files.folder}, {files.name});

nFiles = numel(uffFiles);

fprintf('Fant %d UFF-filer.\n', nFiles);

for i = 1:nFiles
    fprintf('%2d: %s\n', i, files(i).name);
end

%% Scan
invivo_scan = uff.read_object(fullfile(local_path, 'InVivoRatBrain_scan.uff'), '/scan');
%%
invivo_ch_data = uff.read_object(fullfile(input_folder, 'InVivoRatBrain_001.uff'), '/channel_data');


%%

invivo_scan = tools.scan_integer_upscale(invivo_scan, 2);

% DAS
das = midprocess.das();
das.dimension = dimension.transmit;
das.scan = invivo_scan;
das.receive_apodization.f_number = 1.5; 
das.receive_apodization.window = uff.window.hamming;

% CF
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
% SVD
windowFiles = 3;
cutoff = 5;

for targetFileIndex = 1:nFiles
    [svd_ch_data_all, svdModel] = tools.slidingTemporalSVD(uffFiles, targetFileIndex, windowFiles, cutoff);

    total_frames = size(svd_ch_data_all.data, 4);
    batch_size = 40;
    invivo_b_cf_full = [];

    for start_frame = 1:batch_size:total_frames
        end_frame = min(start_frame + batch_size - 1, total_frames);
        fprintf('Processing frames %d to %d of %d...\n', start_frame, end_frame, total_frames);
        
        % Slice current batch of frames
        %invivo_ch_data.data = raw_frames(:,:,:, start_frame:end_frame);
        frame_idx = start_frame:end_frame;
        % SVD
    
        svd_ch_data_batch = uff.channel_data(svd_ch_data_all);
        svd_ch_data_batch.data = svd_ch_data_all.data(:,:,:,frame_idx);
    
    
        % DAS
        %das.channel_data = svd_ch_data;
        das.channel_data = svd_ch_data_batch;
        b_data_tx = das.go();
    
        % CF
        cf.input = b_data_tx;
        invivo_b_cf = cf.go();
        invivo_b_cf.frame_rate = 1000; % WHY from param.mat?
    
        % Coherent compounding (ADD)
        % das_rx.input = b_data_tx;
        % invivo_b_das = das_rx.go();
    
        if isempty(invivo_b_cf_full)
            invivo_b_cf_full = invivo_b_cf;
            %invivo_b_das_full = invivo_b_das;
        else
            invivo_b_cf_full.data = cat(4, invivo_b_cf_full.data, invivo_b_cf.data);
            %invivo_b_das_full.data = cat(4, invivo_b_das_full.data, invivo_b_das.data);
        end
    
        clear svd_ch_data b_data_tx invivo_b_cf;
    
    end
    [~, fileName] = fileparts(uffFiles{targetFileIndex});

    beamformed_path_CF = fullfile(output_folder, ...
        [fileName '_beamformed_CF_slidingSVD.uff']);

    uff.write_object(beamformed_path_CF, ...
        invivo_b_cf_full, ...
        'b_data');

    save(fullfile(output_folder, ...
        [fileName '_slidingSVD_model.mat']), ...
        'svdModel');

    clear svd_ch_data_all invivo_b_cf_full

end

%%
beamformed_file_1 = fullfile(output_folder, 'InVivoRatBrain_001_beamformed_CF_slidingSVD.uff' );
beamformed_file_2 = fullfile(output_folder, 'InVivoRatBrain_002_beamformed_CF_slidingSVD.uff' );

invivo_b_cf_full_1 = uff.read_object(beamformed_file_1, '/b_data');
invivo_b_cf_full_2 = uff.read_object(beamformed_file_2, '/b_data');
%%
figure;
invivo_b_cf_full.plot([], 'InVivo Ratbrain Sliding SVD CF', 60);

%%
invivo_b_cf_full = invivo_b_cf_full_1;
invivo_b_cf_full.data = cat(4, invivo_b_cf_full_1.data, invivo_b_cf_full_2.data);

%%
u = ulm.ULM();
u.framerate = 500;
u.fwhm = [5 5];
u.numberOfParticles = 40;
u.NLocalMax = 2;
u.max_linking_distance = 3; 
u.min_length = 15;
u.algorithm = ulm.algorithm.radial;
u.tracking = ulm.tracking.tracks;
u.lambda = invivo_ch_data.lambda;
u.input = invivo_b_cf_full;
u.scan = invivo_scan;

% Then simply execute the ULM process
tracks = u.go();


%%
figure_path = fullfile(ustb_path(), 'examples', 'ulm', 'Figures');

ulm_img = u.create_image(tracks, "tracks");
fig1 = figure('Visible', 'off');
imagesc(invivo_scan.x_axis * 1e3, invivo_scan.z_axis * 1e3, ulm_img);
xlabel("X [mm]");ylabel("z [mm]");
title('ULM InVivo Rat Brain')
colormap turbo;
exportgraphics(gca, fullfile(figure_path, 'ULM_InVivo_Rat_Brain_sliding_SVD.png'));

%%
u.tracking = ulm.tracking.velocity_interpolation;
ulm_img = u.create_image(u.go(), "tracks");
fig2 = figure('Visible', 'off');
imagesc(invivo_scan.x_axis * 1e3, invivo_scan.z_axis * 1e3, ulm_img);
xlabel("X [mm]");ylabel("z [mm]");
title('ULM InVivo Rat Brain with interpolation')
colormap turbo;
exportgraphics(gca, fullfile(figure_path, 'ULM_InVivo_Rat_Brain_interpolated_sliding_SVD.png'));


