%% Step 1: Data Loading
% Start by loading the first chunk of channel data and the uff scan object 
local_path = fullfile(ustb_path(), 'data');
RF_channeldata_path  = fullfile(local_path, 'RF_channeldata');

% Creating file for generated figures
figure_path = fullfile(ustb_path(), 'examples', 'ulm', 'Figures');
if ~exist(figure_path, 'dir')
    mkdir(figure_path);
end

invivo_ch_data = uff.read_object(fullfile(RF_channeldata_path, 'InVivoRatBrain_001.uff'), '/channel_data');
invivo_scan = uff.read_object(fullfile(local_path, 'InVivoRatBrain_scan.uff'), '/scan');

fprintf("\nInVivo Data:\nSamples: %d\nReceive: %d\nTransmits: %d\nFrames: %d\n", deal(size(invivo_ch_data.data)));


%% Step 2: Set up Beamforming 
% For this example, the beamforming step uses DAS only on transmit and
% combines the receive with the coherence-factor beamformer. The parameters
% used follows the tuned parameter case for lambda/2 sampling from Simon A.
% Bjørn's master's thesis, Figure 5.1.
% http://hdl.handle.net.ezproxy.uio.no/10852/120821

% Since the scan is by default using lambda-by-lamture bda sampling, upscaling
% the scan resolution by 2, yields a lambda/2-by-lambda/2 spatial sampling.
% A convinience script for this was added to the +tools module
invivo_scan = tools.scan_integer_upscale(invivo_scan, 2);

% Preprocess SVD-filter
svd = preprocess.svd_filter();
svd.cutoff = 2;

% DAS
das = midprocess.das();
das.dimension = dimension.transmit;
das.scan = invivo_scan;
das.receive_apodization.f_number = 1.5; 
das.receive_apodization.window = uff.window.hamming;

% CF
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;

% Coherent Compounding (ADD)
% das_rx = postprocess.coherent_compounding();
% das_rx.dimension = dimension.receive;

%% Step 3: Beamforming batches
beamformed_path_CF = fullfile(local_path, 'beamformed_RatBrain_CF_001.uff');
%beamformed_path_DAS = fullfile(local_path, 'beamformed_RatBrain_DAS_001.uff');

raw_frames = invivo_ch_data.data;
total_frames = size(invivo_ch_data.data, 4);
batch_size = 40;

for start_frame = 1:batch_size:total_frames
    end_frame = min(start_frame + batch_size - 1, total_frames);
    fprintf('Processing frames %d to %d of %d...\n', start_frame, end_frame, total_frames);

    % Slice current batch of frames
    invivo_ch_data.data = raw_frames(:,:,:, start_frame:end_frame);
  
    % SVD
    svd.input = invivo_ch_data;
    svd_ch_data = svd.go();

    % DAS
    das.channel_data = svd_ch_data;
    b_data_tx = das.go();

    % CF
    cf.input = b_data_tx;
    invivo_b_cf = cf.go();
    invivo_b_cf.frame_rate = 100;

    % Coherent compounding (ADD)
    % das_rx.input = b_data_tx;
    % invivo_b_das = das_rx.go();

    if start_frame == 1
        invivo_b_cf_full = invivo_b_cf;
        %invivo_b_das_full = invivo_b_das;
    else
        invivo_b_cf_full.data = cat(4, invivo_b_cf_full.data, invivo_b_cf.data);
        %invivo_b_das_full.data = cat(4, invivo_b_das_full.data, invivo_b_das.data);
    end

    clear svd_ch_data b_data_tx invivo_b_cf;

end

uff.write_object(beamformed_path_CF, invivo_b_cf_full, 'b_data');
% uff.write_object(beamformed_path_das, invivo_b_das_full, 'b_data');


%% Visualization
% CF
figure;
invivo_b_cf_full.plot([], 'InVivo Rat Brain SVD CF', 60);

% Coherent Compounding (ADD)

%% Step 4: The actual ULM part
% With some data beamformed and configured, this step walks through a basic
% configurationg of a ULM pipeline.
% The basis for all pipelines is the ULM handle object. This creature
% behaves similar to postprocesses, where a uff.beamformed_data is fed into
% its 'input' property, and 'go()' executes its algorithm. However, there
% are quite a few parameters and options that goes along with it. Some of
% which will be covered here.

% Start by instantiating a simple ULM object:
u = ulm.ULM();

% The following parameters used in chapter 5 of Simon A.
% Bjørn's master's thesis, Table 5.1b: CF @ lambda/2:
% The framerate specifies the framerate at which the RF data is captured.
% This is required due its implications on linking particales across
% frames.
u.framerate = 500;

% The full-width half-maximum (fwhm) parameter tunes the kernel sizes of
% initial particle position guesses, and is configured in number of pixels.
% For this case, the fwhm is estimated to be 3x3 (width, height) pixels.
u.fwhm = [3 3];

% The numberOfParticles parameter sets the upper limit for how many
% particles the ULM process tries to localize. If more than this number is
% located, only the specified number of strongest points will be recorded.
% This is per frame without heuristics, so overshooting a bit is better, as
% this allows for a small amount of false positives in addition to true 
% positives, instead of potentially discarding true positives. For this
% pure example, however, we know there are a maximum of 41 particles in 
% any given frame, meaning its a good estimate to use here.
u.numberOfParticles = 40;

% The last 3 parameters are a bit more advanced and not well documented in
% the toolbox yet. For now, see Chapter 3.3 "ULM implementation in USTB" 
% of Simon A. Bjørn's master's thesis.
u.NLocalMax = 2;
u.max_linking_distance = 3; 
u.min_length = 15;


% The next two options are not parameters, but settings deciding how the
% ULM process is performed. 

% The first is the algorithm, which specifies
% which localization algorithm to use. The radial symmetry localization
% algorithm is implemented under the ulm.algorithm.radial, and is quite
% fast and efficient, so its used in this example. For other options, see
% enumeration('ulm.algorithm')
u.algorithm = ulm.algorithm.radial;

% The second option specifies how the linker operates. In this example, we
% want to compare the raw localized datapoints against ground truth, and so
% we do not want any interpolation of the trackes based on velocity or
% positions. This equates to using the "tracks" tracking algorithm, which
% simply outputs a cell array, containing arrays of particle tracks over
% with positions at each frame the track exists.
u.tracking = ulm.tracking.tracks;

% Lastly, the data from the previous step is supplied. Lambda must be
% supplied separately, as beamformed_data has no lambda property.
u.lambda = invivo_ch_data.lambda;
u.input = invivo_b_cf_full;
u.scan = invivo_scan;

% Then simply execute the ULM process
tracks = u.go();

%% Step 5a: ULM Image construction
% ULM is no fun without images.
% To synthesize an ULM image from tracks, the aptly named "create_image"
% method in the ULM process can be used. Simply feed it the tracks created
% from the main "go()" process.

% Synthesize a "track" image. Other image modes are available. See
% enumeration("ulm.image_mode") for more
ulm_img = u.create_image(tracks, "tracks");
fig1 = figure('Visible', 'off');
imagesc(invivo_scan.x_axis * 1e3, invivo_scan.z_axis * 1e3, ulm_img);
xlabel("X [mm]");ylabel("z [mm]");
title('ULM InVivo Rat Brain')
colormap turbo;
exportgraphics(gca, fullfile(figure_path, 'ULM_InVivo_Rat_Brain.png'));


%% Step 5b: ULM Image Construction with interpolation
% As is visible, the image is quite choppy, and the tracks are pixelated.
% This is because the tracking algorithm only used the particle positions
% at every "full frame". This is especially visable for fast moving 
% particles, as large gaps are formed as the particle moves multiple pixels
% between frames. No interpolation, no smooth paths.
% Changing the tracking algorithm to use velocity_interpolation instead
% yields a much better image.

u.tracking = ulm.tracking.velocity_interpolation;
ulm_img = u.create_image(u.go(), "tracks");
fig2 = figure('Visible', 'off');
imagesc(invivo_scan.x_axis * 1e3, invivo_scan.z_axis * 1e3, ulm_img);
xlabel("X [mm]");ylabel("z [mm]");
title('ULM InVivo Rat Brain with interpolation')
colormap turbo;
exportgraphics(gca, fullfile(figure_path, 'ULM_InVivo_Rat_Brain_interpolated.png'));


fprintf('\nSaved output images to:\n - %s\n - %s\n - %s\n - %s\n - %s\n', ...
    [figure_path 'Beamformed_DAS_SVD_filtered.png'], ...
    [figure_path 'Beamformed_CF_SVD_filtered.png'], ...
    [figure_path 'Beamformed_CF_SVD_filtered_MIP.png'], ...
    [figure_path 'ULM_InVivo_Rat_Brain.png'], ...
    [figure_path 'ULM_InVivo_Rat_Brain_interpolated.png']);

