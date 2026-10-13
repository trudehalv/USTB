%% Read OpenH-RF ULMShare plane-wave IQ data into USTB
%
% Converts one ULMShare acquisition (transcranial mouse-brain ULM, Verasonics
% Vantage 256 with an L22-14v probe, plane-wave IQ) from the OpenH-RF zea
% format to a uff.channel_data object, beamforms it with USTB DAS on the
% ULMShare reference grid, and shows a B-mode image and a power-Doppler image
% after SVD clutter filtering.
%
% Dataset: https://huggingface.co/datasets/nvidia/OpenH-RF/tree/main/ulmshare
% Explorer: https://open-h.github.io/OpenH-RF/#dataset/ulmshare
% Citation: Rauby et al. (2026), ULMShare, doi:10.20383/103.01550 (CC BY 4.0)
%
% The full acquisitions are 60-370 GB, so first stream one 400-frame buffer
% (~720 MB local) with the helper in this folder:
%
%   python openhrf_download_subset.py ^
%     --source hf://nvidia/OpenH-RF/ulmshare/data/mouse_10_acquisition_1_zea.hdf5 ^
%     --output C:/Data/OpenH-RF/ulmshare/mouse_10_acquisition_1_frames_20000_20399.hdf5 ^
%     --first-frame 20000 --n-frames 400
%
% (needs: pip install h5py hdf5plugin huggingface_hub numpy)
%
% Frame 20000 (buffer 50) is used because the first buffers of an acquisition
% are typically recorded before the microbubble bolus arrives (e.g. frames
% 8000-8399 of mouse 10 show tissue and noise only in power Doppler).

clear all; close all;
local_file = fullfile(ustb_path(), '/sandbox/OpenHRF/Data/OpenH-RF/ulmshare/mouse_10_acquisition_1/buffer_050.hdf5' );
%local_file = 'C:/Data/OpenH-RF/ulmshare/mouse_10_acquisition_1_frames_20000_20399.hdf5';

%% Convert to UFF
channel_data = openhrf_to_uff(local_file);
disp(channel_data.name.');

%% ULMShare reference beamforming grid (25 um pixels)
scan = uff.linear_scan();
scan.x_axis = (-5e-3:25e-6:5e-3-25e-6).';
scan.z_axis = (0.5e-3:25e-6:8.5e-3-25e-6).';


%%
% Pre-beamforming SVD on all frames in this buffer
svd = preprocess.svd_filter();
svd.input = channel_data;
svd.cutoff = 5; % How many?
svd_channel_data = svd.go();

% DAS: compound transmit plane waves, preserve receive dimension
das = midprocess.das();
das.scan = scan;
das.dimension = dimension.transmit;
das.receive_apodization.window = uff.window.hanning;
das.receive_apodization.f_number = 1.4;
das.transmit_apodization.window = uff.window.none;

% Coherence factor over receive aperture
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;

% Ordinary coherent receive compounding
das_rx = postprocess.coherent_compounding();
das_rx.dimension = dimension.receive;

n_frames = size(svd_channel_data.data, 4);
batch_size = 10;

% Store final, compact outputs.
% These will contain only image pixels × frames, not 128 receive channels.
b_data_das_all = [];
b_data_cf_all  = [];

batch_number = 0;

for start_frame = 1:batch_size:n_frames
    end_frame = min(start_frame + batch_size - 1, n_frames);
    frame_idx = start_frame:end_frame;
    batch_number = batch_number + 1;

    fprintf('Frames %d–%d of %d\n', ...
        start_frame, end_frame, n_frames);

    % Extract only the current frame batch
    ch_batch = uff.channel_data(svd_channel_data);
    ch_batch.data = svd_channel_data.data(:, :, :, frame_idx);

    % DAS over transmit, while retaining receive information
    das.channel_data = ch_batch;
    b_data_tx = das.go();

    % CF image
    cf.input = b_data_tx;
    b_data_cf_batch = cf.go();
    %b_data_cf_batch.frame_rate = 1000;

    % Standard coherent-compounded DAS image
    das_rx.input = b_data_tx;
    b_data_das_batch = das_rx.go();
    %b_data_das_batch.frame_rate = 1000;

    if isempty(b_data_cf_all)
        b_data_cf_all = uff.beamformed_data(b_data_cf_batch);
        b_data_das_all = uff.beamformed_data(b_data_das_batch);

        b_data_cf_all.data = zeros( ...
            [size(b_data_cf_batch.data, 1), ...
             size(b_data_cf_batch.data, 2), ...
             size(b_data_cf_batch.data, 3), ...
             n_frames], ...
            'like', b_data_cf_batch.data);

        b_data_das_all.data = zeros( ...
            [size(b_data_das_batch.data, 1), ...
             size(b_data_das_batch.data, 2), ...
             size(b_data_das_batch.data, 3), ...
             n_frames], ...
            'like', b_data_das_batch.data);
    end
    % Insert this batch at its original frame positions
    b_data_cf_all.data(:,:, :, frame_idx) = b_data_cf_batch.data;
    b_data_das_all.data(:, :,:, frame_idx) = b_data_das_batch.data;

    % Critical: release receive-resolved DAS data before next batch
    clear ch_batch b_data_tx b_data_cf_batch b_data_das_batch
end


%% B-mode DAS of the first frame
b_mode_das_all = uff.beamformed_data(b_data_das_all);
b_mode_das_all.data = b_data_das_all.data(:, :, :, 1);
b_mode_das_all.plot([], 'ULMShare mouse 10 (OpenH-RF) - DAS', 50);

%% B-mode CF of the first frame
b_mode_cf_all = uff.beamformed_data(b_data_cf_all);
b_mode_cf_all.data = b_data_cf_all.data(:, :, :, 1);
b_mode_cf_all.plot([], 'ULMShare mouse 10 (OpenH-RF) - CF', 50);
%% Power Doppler after SVD clutter filtering of the beamformed frames
n_frames = b_data_das.N_frames;
casorati = reshape(b_data_das.data, [], n_frames);              % [pixel x frame]
n_cut = max(round(0.05*n_frames) - 1, 0);                    % ULMShare: reject 5% of the components
[U, S, V] = svd(casorati, 'econ');
filtered = casorati - U(:, 1:n_cut)*S(1:n_cut, 1:n_cut)*V(:, 1:n_cut)';
power_doppler = reshape(mean(abs(filtered).^2, 2), scan.N_z_axis, scan.N_x_axis);

figure;
imagesc(scan.x_axis*1e3, scan.z_axis*1e3, 10*log10(power_doppler./max(power_doppler(:))));
axis image; colormap hot; colorbar; clim([-40 0]);
xlabel('x [mm]'); ylabel('z [mm]');
title(sprintf('ULMShare mouse 10 (OpenH-RF) - power Doppler, %d frames', n_frames));
