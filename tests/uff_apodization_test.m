classdef uff_apodization_test < matlab.unittest.TestCase

    methods (Test)
        function test_default_construction(testCase)
            apo = uff.apodization();
            testCase.verifyClass(apo, ?uff.apodization);
        end

        function test_default_window_is_none(testCase)
            apo = uff.apodization();
            testCase.verifyEqual(apo.window, uff.window.none);
        end

        function test_default_f_number(testCase)
            apo = uff.apodization();
            testCase.verifyEqual(apo.f_number, [1, 1]);
        end

        function test_set_window_types(testCase)
            apo = uff.apodization();

            apo.window = uff.window.hamming;
            testCase.verifyEqual(apo.window, uff.window.hamming);

            apo.window = uff.window.hanning;
            testCase.verifyEqual(apo.window, uff.window.hanning);

            apo.window = uff.window.boxcar;
            testCase.verifyEqual(apo.window, uff.window.boxcar);

            apo.window = uff.window.tukey25;
            testCase.verifyEqual(apo.window, uff.window.tukey25);
        end

        function test_set_f_number(testCase)
            apo = uff.apodization();
            apo.f_number = [1.75, 1.75];
            testCase.verifyEqual(apo.f_number, [1.75, 1.75]);
        end

        function test_probe_assignment(testCase)
            prb = uff.linear_array();
            prb.N = 64;
            prb.pitch = 300e-6;

            apo = uff.apodization();
            apo.probe = prb;
            testCase.verifyClass(apo.probe, ?uff.linear_array);
        end

        function test_receive_data_for_every_window(testCase)
            % Computes the data (not just the property) for every window
            % in the enum, so a window without a case in apply_window fails
            prb = uff.linear_array('N', 65, 'pitch', 300e-6);  % odd: one element at x = 0
            scn = uff.linear_scan('x_axis', linspace(-10e-3, 10e-3, 21).', ...
                'z_axis', linspace(1e-3, 40e-3, 40).');
            for win = windows_with_f_number()
                apo = uff.apodization('probe', prb, 'focus', scn, 'window', win, 'f_number', 1.7);
                data = apo.data;
                testCase.verifySize(data, [scn.N_pixels, prb.N_elements], char(win));
                testCase.verifyTrue(all(isfinite(data(:))), char(win));
                testCase.verifyGreaterThanOrEqual(min(data(:)), 0, char(win));
                testCase.verifyLessThanOrEqual(max(data(:)), 1, char(win));
                % The element right above a pixel gets full weight
                [~, el] = min(abs(prb.x - 0));
                px = find(abs(scn.x) < eps & abs(scn.z - 20e-3) < 1e-4, 1);
                testCase.verifyEqual(data(px, el), 1, 'AbsTol', 1e-2, char(win));
            end
        end

        function test_transmit_data_for_every_window(testCase)
            scn = uff.linear_scan('x_axis', linspace(-10e-3, 10e-3, 21).', ...
                'z_axis', linspace(1e-3, 40e-3, 40).');
            seq = uff.wave();
            angles = [-0.2, 0, 0.2];
            for n = 1:numel(angles)
                seq(n) = uff.wave();
                seq(n).wavefront = uff.wavefront.plane;
                seq(n).source.azimuth = angles(n);
                seq(n).source.distance = Inf;
            end
            for win = windows_with_f_number()
                apo = uff.apodization('sequence', seq, 'focus', scn, 'window', win, 'f_number', 1.7);
                data = apo.data;
                testCase.verifySize(data, [scn.N_pixels, numel(seq)], char(win));
                testCase.verifyTrue(all(isfinite(data(:))), char(win));
                % The unsteered plane wave gets full weight everywhere
                testCase.verifyEqual(data(:, 2), ones(scn.N_pixels, 1), 'AbsTol', 1e-12, char(win));
            end
        end
    end

end

function wins = windows_with_f_number()
% All uff.window values computed by apply_window (none, sta and scanline
% have their own code paths)
wins = enumeration('uff.window').';
wins = wins(~ismember(wins, [uff.window.none, uff.window.sta, uff.window.scanline]));
end
