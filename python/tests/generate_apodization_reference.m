% generate_apodization_reference.m
% Saves MATLAB uff.apodization matrices for the configurations compared in
% test_apodization_vs_matlab.py (transmit: plane and diverging waves;
% receive: linear and phased array, linear and sector scan).
% Run from the repository root:
%   matlab -batch "addpath('.'); run('python/tests/generate_apodization_reference.m');"

addpath(ustb_path());
out = fullfile(fileparts(mfilename('fullpath')), 'apodization_reference.h5');
if isfile(out), delete(out); end
d = [ustb_path() '/data/'];
url = tools.zenodo_dataset_files_base();
tools.download('PICMUS_experiment_resolution_distortion.uff', url, d);
tools.download('Verasonics_P2-4_parasternal_long_small.uff', url, d);
pic = uff.read_object([d 'PICMUS_experiment_resolution_distortion.uff'], '/channel_data');
pscan = uff.read_object([d 'PICMUS_experiment_resolution_distortion.uff'], '/scan');
w = @(name, v) h5save(out, name, v);
w('/linear/x', pscan.x); w('/linear/z', pscan.z);

% transmit: PICMUS plane waves, every window
wins = {'boxcar','hanning','hamming','tukey25','tukey50','tukey75','triangle'};
for k = 1:numel(wins)
    a = uff.apodization('sequence', pic.sequence, 'focus', pscan);
    a.window = uff.window.(wins{k}); a.f_number = 1.7;
    w(['/tx_plane/' wins{k}], a.data);
end
a = uff.apodization('sequence', pic.sequence, 'focus', pscan, 'window', uff.window.tukey50, 'f_number', 1.7, 'tilt', [0.05 0]);
w('/tx_plane/tukey50_tilt', a.data);

% transmit: diverging waves (virtual sources behind the probe) and an STA-like wave
seq = uff.wave();
xs = [-10 -3 0 4 12]*1e-3;
for n = 1:numel(xs)
    seq(n) = uff.wave(); seq(n).wavefront = uff.wavefront.spherical;
    seq(n).source.xyz = [xs(n) 0 -15e-3]; seq(n).origin.xyz = [xs(n) 0 0];
end
seq(end+1) = uff.wave(); seq(end).wavefront = uff.wavefront.spherical; seq(end).source.xyz = [2e-3 0 0]; seq(end).origin.xyz = [2e-3 0 0];
src = reshape([seq.source], 1, []); org = reshape([seq.origin], 1, []);
w('/tx_div/source', [src.x; src.y; src.z]); w('/tx_div/origin', [org.x; org.y; org.z]);
for k = {'hanning','tukey50','boxcar'}
    a = uff.apodization('sequence', seq, 'focus', pscan, 'window', uff.window.(k{1}), 'f_number', 1.2);
    w(['/tx_div/' k{1}], a.data);
end
a = uff.apodization('sequence', seq, 'focus', pscan, 'window', uff.window.tukey50, 'f_number', 1.2, 'grating_lobe_angle', 0.6, 'minimum_aperture', 3e-3, 'maximum_aperture', 20e-3);
w('/tx_div/tukey50_limits', a.data);

% receive: PICMUS linear probe and linear scan
for k = 1:numel(wins)
    a = uff.apodization('probe', pic.probe, 'focus', pscan, 'window', uff.window.(wins{k}), 'f_number', 1.7);
    w(['/rx_linear/' wins{k}], a.data);
end
a = uff.apodization('probe', pic.probe, 'focus', pscan, 'window', uff.window.hamming, 'f_number', [1.4 1.4], 'tilt', [0.1 0], 'minimum_aperture', 4e-3, 'maximum_aperture', 15e-3);
w('/rx_linear/hamming_tilt_limits', a.data);

% receive: phased array with a sector scan (origin at the apex)
ver = uff.read_object([d 'Verasonics_P2-4_parasternal_long_small.uff'], '/channel_data');
az = linspace(ver.sequence(1).source.azimuth, ver.sequence(end).source.azimuth, ver.N_waves).';
sscan = uff.sector_scan('azimuth_axis', az, 'depth_axis', linspace(1e-3, 110e-3, 64).');
w('/sector/x', sscan.x); w('/sector/y', sscan.y); w('/sector/z', sscan.z);
w('/sector/azimuth_axis', az); w('/sector/depth_axis', sscan.depth_axis);
for k = {'hanning','tukey25'}
    a = uff.apodization('probe', ver.probe, 'focus', sscan, 'window', uff.window.(k{1}), 'f_number', 1.5);
    w(['/sector/' k{1}], a.data);
end
% transmit: scanline (MLA) apodization, sector scans starting at depth 0
s1 = uff.sector_scan('azimuth_axis', az, 'depth_axis', linspace(0, 110e-3, 32).');
s2 = uff.sector_scan('azimuth_axis', linspace(az(1), az(end), 2*numel(az)).', 'depth_axis', s1.depth_axis);
w('/scanline/sector1/x', s1.x); w('/scanline/sector1/z', s1.z); w('/scanline/sector1/azimuth_axis', s1.azimuth_axis);
w('/scanline/sector2/x', s2.x); w('/scanline/sector2/z', s2.z); w('/scanline/sector2/azimuth_axis', s2.azimuth_axis);
w('/scanline/depth_axis', s1.depth_axis);
a = uff.apodization('sequence', ver.sequence, 'focus', s1, 'window', uff.window.scanline);
w('/scanline/sector_mla1', a.data);
a = uff.apodization('sequence', ver.sequence, 'focus', s2, 'window', uff.window.scanline, 'MLA', [2 1]);
w('/scanline/sector_mla2', a.data);
a = uff.apodization('sequence', ver.sequence, 'focus', s2, 'window', uff.window.scanline, 'MLA', [2 1], 'MLA_overlap', [1 0]);
w('/scanline/sector_mla2_overlap1', a.data);

% transmit: scanline apodization on a linear scan (only the number of waves matters)
lscan = uff.linear_scan('x_axis', linspace(-8e-3, 8e-3, 32).', 'z_axis', linspace(0, 30e-3, 20).');
w('/scanline/linear/x', lscan.x); w('/scanline/linear/z', lscan.z);
w('/scanline/linear/x_axis', lscan.x_axis); w('/scanline/linear/z_axis', lscan.z_axis);
lseq = repmat(uff.wave(), 1, 32);
a = uff.apodization('sequence', lseq, 'focus', lscan, 'window', uff.window.scanline);
w('/scanline/linear_mla1', a.data);
for ov = [0 1 2]
    a = uff.apodization('sequence', lseq(1:8), 'focus', lscan, 'window', uff.window.scanline, 'MLA', [4 1], 'MLA_overlap', [ov 0]);
    w(sprintf('/scanline/linear_mla4_overlap%d', ov), a.data);
end

fprintf('saved %s\n', out);

function h5save(f, name, v)
    h5create(f, name, size(v)); h5write(f, name, double(v));
end
