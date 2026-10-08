function [filtered_ch_data, model] = slidingTemporalSVD( ...
    uffFiles, targetFileIndex, windowFiles, cutoff)
%SLIDINGTEMPORALSVD Sliding-window temporal SVD for UFF channel data.
%
% Example:
%   [chFiltered, model] = tools.slidingTemporalSVD( ...
%       uffFiles, 5, 3, 5);
%
% targetFileIndex = 5 and windowFiles = 3:
%
%   files [4 5 6] --> one local temporal SVD --> filtered file 5
%
% The function reads /channel_data from each file in the local window once.

    %% Input checks

    uffFiles = cellstr(uffFiles);
    nFiles = numel(uffFiles);

    assert(nFiles >= 1, 'No UFF files were provided.');

    assert(targetFileIndex >= 1 && targetFileIndex <= nFiles, ...
        'targetFileIndex is outside the valid file range.');

    assert(windowFiles >= 1 && windowFiles <= nFiles && ...
        mod(windowFiles, 1) == 0, ...
        'windowFiles must be an integer between 1 and nFiles.');

    assert(cutoff >= 1 && mod(cutoff, 1) == 0, ...
        'cutoff must be a positive integer.');

    %% Select the local sliding window

    halfWindow = floor(windowFiles / 2);

    firstFile = targetFileIndex - halfWindow;
    firstFile = max(firstFile, 1);
    firstFile = min(firstFile, nFiles - windowFiles + 1);

    lastFile = firstFile + windowFiles - 1;

    windowIndices = firstFile:lastFile;
    targetPosition = targetFileIndex - firstFile + 1;

    windowUffFiles = uffFiles(windowIndices);

    fprintf('\n--- Sliding-window temporal SVD ---\n');
    fprintf('Target file index: %d of %d\n', targetFileIndex, nFiles);
    fprintf('Window file indices: %s\n', mat2str(windowIndices));

    %% Read each channel_data object ONCE

    channelDataWindow = cell(windowFiles, 1);

    for i = 1:windowFiles

        fprintf('Reading channel_data: %s\n', windowUffFiles{i});

        channelDataWindow{i} = uff.read_object( ...
            windowUffFiles{i}, '/channel_data');
    end

    %% Determine data dimensions

    target_ch_data = channelDataWindow{targetPosition};
    targetData = target_ch_data.data;

    originalSize = size(targetData);

    % Expected layout:
    % [samples x channels x waves x frames]
    nFrames = size(targetData, 4);
    nFeatures = numel(targetData) / nFrames;
    nWindowFrames = windowFiles * nFrames;

    assert(cutoff < nWindowFrames, ...
        'cutoff must be smaller than number of frames in the full window.');

    % Verify equal dimensions for all files in the window.
    for i = 1:windowFiles

        assert(isequal(size(channelDataWindow{i}.data), originalSize), ...
            'All files in one SVD window must have equal dimensions.');
    end

    fprintf('Features per frame: %d\n', nFeatures);
    fprintf('Frames per file:    %d\n', nFrames);
    fprintf('Frames in window:   %d\n', nWindowFrames);
    fprintf('Removed ranks:      1:%d\n', cutoff);

    %% Reshape each file into [features x frames]

    Xwindow = cell(windowFiles, 1);

    for i = 1:windowFiles

        Xwindow{i} = reshape( ...
            channelDataWindow{i}.data, ...
            nFeatures, ...
            nFrames);
    end

    %% Build temporal Gram matrix
    %
    % Conceptually:
    %
    % Xall = [X1 X2 ... XW]
    % Gtime = Xall' * Xall
    %
    % Here, Gtime is built from smaller blocks.

    Gtime = zeros(nWindowFrames, nWindowFrames, 'like', targetData);

    for i = 1:windowFiles

        rowIndices = (i - 1) * nFrames + (1:nFrames);

        for j = i:windowFiles

            colIndices = (j - 1) * nFrames + (1:nFrames);

            gramBlock = Xwindow{i}' * Xwindow{j};

            Gtime(rowIndices, colIndices) = gramBlock;

            if i ~= j
                Gtime(colIndices, rowIndices) = gramBlock';
            end
        end
    end

    % Avoid small numerical asymmetries.
    Gtime = (Gtime + Gtime') / 2;

    %% Eigendecomposition of the temporal Gram matrix

    [V, lambda] = eig(Gtime, 'vector');

    [lambda, order] = sort(real(lambda), 'descend');

    lambda = max(lambda, 0);
    V = V(:, order);

    removeRanks = 1:cutoff;
    Vclutter = V(:, removeRanks);

    %% Calculate clutter basis coefficients
    %
    % A = Xall * Vclutter
    %
    % calculated blockwise:
    %
    % A = X1*V1 + X2*V2 + ... + XW*VW

    clutterCoefficients = zeros( ...
        nFeatures, cutoff, 'like', targetData);

    for i = 1:windowFiles

        temporalRows = (i - 1) * nFrames + (1:nFrames);

        V_i = Vclutter(temporalRows, :);

        clutterCoefficients = clutterCoefficients + ...
            Xwindow{i} * V_i;
    end

    %% Filter only the requested target file

    Xtarget = Xwindow{targetPosition};

    targetRows = (targetPosition - 1) * nFrames + (1:nFrames);

    Vtarget = Vclutter(targetRows, :);

    % Remove dominant temporal modes, assumed to be clutter.
    Xfiltered = Xtarget - clutterCoefficients * Vtarget';

    %% Restore the original UFF channel_data layout

    filtered_ch_data = uff.channel_data(target_ch_data);

    filtered_ch_data.data = reshape(Xfiltered, originalSize);

    %% Diagnostic output

    model.targetFileIndex = targetFileIndex;
    model.windowIndices = windowIndices;
    model.windowUffFiles = windowUffFiles;
    model.targetPositionInWindow = targetPosition;
    model.windowFiles = windowFiles;
    model.cutoff = cutoff;
    model.removeRanks = removeRanks;
    model.nFeatures = nFeatures;
    model.nFramesPerFile = nFrames;
    model.nWindowFrames = nWindowFrames;
    model.singularValues = sqrt(lambda);

end
