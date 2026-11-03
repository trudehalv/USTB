%% 1) Convert to uff
clear all; close all;

% Define paths
local_file = fullfile(ustb_path(), '/sandbox/OpenHRF/Data/OpenH-RF/ulmshare/mouse_42_acquisition_1/buffer_050.hdf5' );
local_path = fullfile(ustb_path(), '/sandbox/OpenHRF/Data/OpenH-RF/ulmshare/mouse_42_acquisition_1/' );


% Convert to uff
channel_data = openhrf_to_uff(local_file);
% %% Download sequence.json file
% base_url = 'https://www.frdr-dfdr.ca/repo/files/33/published/publication_1817/submitted_data/mouse_18/acquisition_2/';
% seq_url = [base_url 'sequence.json'];
% seq_file = [local_path 'sequence.json'];
% websave(seq_file, seq_url)
% 
% %% Read parameters in json file
% fid = fopen([local_path filesep 'sequence.json'], 'r');
% seq = fscanf(fid, '%c'); fclose(fid); % use fread instead??
% seq = jsondecode(seq);

%%
% Looking at raw channeldata, inspecting singular value spectrum
total_frames = size(channel_data.data, 4); 

% Reshape raw channel data into a 2D Casorati matrix: [(Time * Channels * Emissions) x frames]
casorati = reshape(channel_data.data, [], total_frames);

% Compute SVD
[~, S, V] = svd(casorati, 'econ');
singular_vals = diag(S);

%%
figure;
plot(log10(singular_vals), 'b-', 'LineWidth', 1.5);
grid on;
%xlim([0 50]);
% ylim([-80 0]); % Standard dynamic range matching typical SVD plots
xlabel('Singular value orders');
ylabel('Singular value (Log10)');
title('Mouse ... - Singular Value Spectrum');

%% Undersøk frekvensinnholdet til tidsvektorene (V)
Fs = 1000; % Eksempel: Bildefrekvens (Hz)
N = size(V, 1);
frequencies = (0:N/10)/N * Fs; % Sjekk de første frekvensene

% Gjør FFT av de første 30 tidsvektorene i V
V_fft = abs(fft(V(:, 1:30)));

figure;
imagesc(1:30, frequencies, V_fft(1:length(frequencies), :));
colorbar;
xlabel('Singular value index');
ylabel('Frequency (Hz)');
title('Mouse ... - Frequency per singular value');


%% Make scan object
scan = uff.linear_scan();
scan.x_axis = (-5e-3:25e-6:5e-3-25e-6).';
scan.z_axis = (0.5e-3:25e-6:8.5e-3-25e-6).';
%scan = tools.scan_integer_upscale(scan, 2); % Not needed with the scan grid resolution


%% 2) Beamform
total_frames = size(channel_data.data, 4); 
lambda = channel_data.lambda;

% Preprocess SVD-filter
svd = preprocess.svd_filter();
svd.cutoff = 5; 

% DAS
das = midprocess.das();
das.dimension = dimension.transmit;
das.scan = scan;
das.receive_apodization.f_number = 1.4; 
das.receive_apodization.window = uff.window.hamming;

% CF
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;

% Coherent Compounding 
das_rx = postprocess.coherent_compounding();
das_rx.dimension = dimension.receive;

batch_size = 10;
frame_rate = 1000; 

% Actual beamforming
svd.input = channel_data;
svd_ch_data = svd.go();

for start_frame = 1:batch_size:total_frames
    end_frame = min(start_frame + batch_size - 1, total_frames);
    fprintf('Processing frames %d to %d of %d...\n', start_frame, end_frame, total_frames);
    
    % Slice current batch of frames
    frame_idx = start_frame:end_frame;

    % SVD
    svd_ch_data_batch = uff.channel_data(svd_ch_data); 
    svd_ch_data_batch.data = svd_ch_data.data(:,:,:,frame_idx);

    % DAS
    das.channel_data = svd_ch_data_batch;
    b_tx_batch = das.go();

    % CF
    cf.input = b_tx_batch;
    b_CF_batch = cf.go();
    b_CF_batch.frame_rate = frame_rate; 

    % Coherent compounding 
    das_rx.input = b_tx_batch;
    b_DAS_batch = das_rx.go();
    b_DAS_batch.frame_rate = frame_rate;

    if start_frame == 1
        b_CF_file = b_CF_batch;
        b_DAS_file = b_DAS_batch;
    else
        b_CF_file.data = cat(4, b_CF_file.data, b_CF_batch.data);
        b_DAS_file.data = cat(4, b_DAS_file.data, b_DAS_batch.data);
    end

    clear svd_ch_data_batch b_tx_batch b_CF_batch b_DAS_batch;

end

%% 3) Visualise Beamformed Images
% DAS
figure;
b_DAS_file.plot([], 'ULMShare Mouse 42 (acq.1) (DAS)', 60);

%%
% CF
figure;
b_CF_file.plot([], 'ULMShare Mouse 42 (acq.1) (CF)', 60);

%% ULM

u = ulm.ULM();
u.framerate = frame_rate;
u.fwhm = [5 5];
u.numberOfParticles = 90; 
u.intensity_threshold = 1.5;
u.NLocalMax = 7; % 7 for fwhm = [5 5]
u.max_linking_distance = 3; 
u.min_length = 15;
u.algorithm = ulm.algorithm.radial;
u.tracking = ulm.tracking.velocity_interpolation;
u.lambda = lambda;

scan = b_CF_file.scan;

env = abs(b_CF_file.data);
env = reshape(env, scan.N_z_axis, scan.N_x_axis, []);
env = env./imgaussfilt(sqrt(mean(env.^2, 3)), 9);
b_buffer = uff.beamformed_data();
b_buffer.scan = scan;
b_buffer.data = reshape(env, scan.N_pixels, 1, 1, []);

u.input = b_buffer;
u.scan = scan;
tracks = u.go();
tracks = tracks(:);

%%
% DAS
ulm_img = u.create_image(tracks, "tracks");
fig = figure('Visible', 'on');
imagesc(scan.x_axis * 1e3, scan.z_axis * 1e3, ulm_img);
xlabel("X [mm]");ylabel("z [mm]");
title('Mouse 42 (acq.1) ULM Interpolated - CF')
colormap turbo;
