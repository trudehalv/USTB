classdef uff_probe_test < matlab.unittest.TestCase

    methods (Test)
        function test_linear_array_construction(testCase)
            prb = uff.linear_array();
            prb.N = 128;
            prb.pitch = 300e-6;

            testCase.verifyEqual(prb.N_elements, 128);
        end

        function test_elements_centered_at_origin(testCase)
            prb = uff.linear_array();
            prb.N = 64;
            prb.pitch = 300e-6;

            testCase.verifyEqual(mean(prb.x), 0, 'AbsTol', 1e-12);
        end

        function test_element_spacing_matches_pitch(testCase)
            pitch = 250e-6;
            prb = uff.linear_array();
            prb.N = 32;
            prb.pitch = pitch;

            dx = diff(prb.x);
            testCase.verifyEqual(dx, repmat(pitch, 31, 1), 'AbsTol', 1e-12);
        end

        function test_element_width_defaults_to_pitch(testCase)
            pitch = 300e-6;
            prb = uff.linear_array();
            prb.N = 16;
            prb.pitch = pitch;

            testCase.verifyEqual(prb.element_width, pitch);
        end

        function test_geometry_shape(testCase)
            prb = uff.linear_array();
            prb.N = 64;
            prb.pitch = 300e-6;

            testCase.verifyEqual(size(prb.geometry), [64, 7]);
        end

        function test_elements_on_z_equals_zero(testCase)
            prb = uff.linear_array();
            prb.N = 32;
            prb.pitch = 300e-6;

            testCase.verifyEqual(prb.z, zeros(32, 1));
            testCase.verifyEqual(prb.y, zeros(32, 1));
        end

        function test_linear_array_manual_geometry_locked(testCase)
            prb = uff.linear_array();
            custom_geom = rand(10, 7);
            prb.geometry = custom_geom;

            % Modifying properties that would trigger update() should not overwrite geometry
            prb.pitch = 500e-6;
            prb.N = 32;
            prb.element_width = 400e-6;
            prb.element_height = 5e-3;

            testCase.verifyEqual(prb.N_elements, 10);
            testCase.verifyEqual(prb.geometry, custom_geom);
        end

        function test_linear_array_normal_workflow_updates_geometry(testCase)
            prb = uff.linear_array();
            prb.N = 16;
            prb.pitch = 300e-6;
            testCase.verifyEqual(prb.N_elements, 16);

            % Subsequent modification in normal workflow should recompute geometry
            prb.N = 32;
            testCase.verifyEqual(prb.N_elements, 32);
            testCase.verifyEqual(size(prb.geometry, 1), 32);

            pitch = 200e-6;
            prb.pitch = pitch;
            dx = diff(prb.x);
            testCase.verifyEqual(dx, repmat(pitch, 31, 1), 'AbsTol', 1e-12);
        end

        function test_curvilinear_array_manual_geometry_locked(testCase)
            prb = uff.curvilinear_array();
            custom_geom = rand(8, 7);
            prb.geometry = custom_geom;

            prb.pitch = 300e-6;
            prb.N = 64;
            prb.radius = 40e-3;
            prb.element_width = 250e-6;
            prb.element_height = 4e-3;

            testCase.verifyEqual(prb.N_elements, 8);
            testCase.verifyEqual(prb.geometry, custom_geom);
        end

        function test_matrix_array_manual_geometry_locked(testCase)
            prb = uff.matrix_array();
            custom_geom = rand(12, 7);
            prb.geometry = custom_geom;

            prb.pitch_x = 300e-6;
            prb.pitch_y = 300e-6;
            prb.N_x = 8;
            prb.N_y = 8;
            prb.element_width = 250e-6;
            prb.element_height = 250e-6;

            testCase.verifyEqual(prb.N_elements, 12);
            testCase.verifyEqual(prb.geometry, custom_geom);
        end

        function test_curvilinear_matrix_array_manual_geometry_locked(testCase)
            prb = uff.curvilinear_matrix_array();
            custom_geom = rand(15, 7);
            prb.geometry = custom_geom;

            prb.pitch_x = 300e-6;
            prb.pitch_y = 300e-6;
            prb.N_x = 8;
            prb.N_y = 8;
            prb.radius_x = 40e-3;

            testCase.verifyEqual(prb.N_elements, 15);
            testCase.verifyEqual(prb.geometry, custom_geom);
        end

        function test_manual_geometry_auto_padding(testCase)
            prb = uff.linear_array();
            xyz = rand(6, 3);
            prb.geometry = xyz;

            testCase.verifyEqual(size(prb.geometry), [6, 7]);
            testCase.verifyEqual(prb.geometry(:, 1:3), xyz);
            testCase.verifyEqual(prb.geometry(:, 4:7), zeros(6, 4));

            % Check lock also holds when initialized with <7 columns
            prb.pitch = 300e-6;
            prb.N = 128;
            testCase.verifyEqual(prb.N_elements, 6);
            testCase.verifyEqual(prb.geometry(:, 1:3), xyz);
        end
    end

end
