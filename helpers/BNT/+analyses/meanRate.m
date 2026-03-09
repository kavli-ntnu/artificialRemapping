% Calculate spike mean firing rate
%
% A mean firing rate is obtained by dividing the number of spikes by the duration of a trial/position samples.
% Duration is calculated as the number of tracked position samples multipled by sampling rate. Sampling rate
% is estimated from positions as mean(diff(t)).
%
%  USAGE
%   rate = analyses.meanRate(spikes, positions);
%   spikes      Vector. Spike timestamps.
%   positions   Nx3 or Nx5 ([t x y] or [t x1 y1 x2 y2]) matrix of position samples.
%   rate        Mean firing rate. Units of rate are based on units of time from positions matrix.
%               If time is in seconds, then rate is in Hz.
%
%  EXAMPLE
%
%    data.loadSessions('C:\home\workspace\Li\inputfile-short.cfg');
%    pos = data.getPositions();
%    spikes = data.getSpikeTimes([1 1]);
%    rate = analyses.meanRate(spikes, pos);
%
function rate = meanRate(spikes, positions)
    if isempty(spikes) || isempty(positions)
        rate = 0;
        return;
    end

    sampleTime = mean(diff(positions(:, 1)));
    duration = size(positions, 1) * sampleTime;

    rate = length(spikes) / duration;
end
