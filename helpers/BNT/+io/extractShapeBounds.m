% cameraParams is optional
function [points, width, height] = extractShapeBounds(recording, arenaShape, cameraParams)
    if nargin < 3
        cameraParams = [];
    end
    pos = io.neuralynx.readVideoData(fullfile(recording, 'VT1.nvt'));
    goodInd = ~isnan(pos(:, 2));
    X = pos(goodInd, 2:3);
    if ~isempty(cameraParams)
        X = helpers.undistortPoints(X, cameraParams);
    end

    switch arenaShape
        case bntConstants.ArenaShape.Circle
            error('Not supported');
        case bntConstants.ArenaShape.Track
            error('Not supported');
        case bntConstants.ArenaShape.Box
%             % principal point
%             cx = 357.531774;
%             cy = 274.247960;
%
%             [~, c] = kmeans(X, 4);
%             c = round(c);
%             d = pdist2(c, [0 0]);
%             [~, idx] = min(d); % minimum distance from the origin
%             c_rect(1, :) = c(idx, :);
%
%             selection = true(size(c, 1), 1);
%             selection(idx) = false;
%
%             [~, sIdx] = sort(c(:, 1), 'descend');
%             [~, pointerToPointer] = min(c(sIdx(1:2), 2));
%             idx = sIdx(pointerToPointer);
%             selection(idx) = false;
%             c_rect(2, :) = c(idx, :);
%
%             [~, sIdx] = sort(c(:, 2), 'descend'); % get upper-right corner
%             [~, pointerToPointer] = max(c(sIdx(1:2), 1));
%             idx = sIdx(pointerToPointer);
%             selection(idx) = false;
%             c_rect(3, :) = c(idx, :);
%
%             c_rect(4, :) = c(selection, :);
%
%             width(1) = c_rect(2, 1) - c_rect(1, 1);
%             width(2) = c_rect(3, 1) - c_rect(4, 1);
%             width(3) = abs(width(1) - width(2));
%
%             height(1) = c_rect(4, 2) - c_rect(1, 2);
%             height(2) = c_rect(3, 2) - c_rect(2, 2);
%             height(3) = abs(height(1) - height(2));
%             if width(3) > 0 || height(3) > 0
%                 % adjust rectangle
%             end
            c_rect = general.fitSquare(X);
            height = max(c_rect(:, 2)) - min(c_rect(:, 2));
            width = max(c_rect(:, 1)) - min(c_rect(:, 1));
            if height < 50 || width < 50
                error('BNT:badArenaBound', 'Failed to extract arena bounds from provided information (%s).\nCheck that all LEDs were visible.', recording);
            end

        otherwise
            error('Unknown arena shape');
    end
    points = c_rect;
end