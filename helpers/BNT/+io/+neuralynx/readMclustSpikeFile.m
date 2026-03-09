% readMclustSpikeFile - Read Mclust spike files (t-files)
%
function ts = readMclustSpikeFile(spikeFile)
    % Open file for binary read, big-endian ordering
    fid = data.safefopen(spikeFile, 'rb', 'b');

    beginHeader = '%%BEGINHEADER';
    endHeader = '%%ENDHEADER';

    % Current file position
    curFilePos = ftell(fid);

    str = strtrim(fgets(fid));

    if strcmpi(str, beginHeader)
        while ~feof(fid) && ~strcmp(str, endHeader)
            str = strtrim(fgets(fid));
        end
    else
        % No header
        fseek(fid, curFilePos, 'bof');
    end

    curFilePos = ftell(fid);

    [~, ~, ext] = helpers.fileparts(spikeFile);
    if ext(1) == '.'
        ext(1) = [];
    end
    if length(ext) > 1
        % determine number of bits from the extension
        bits = regexp(ext, '\d*', 'match');
        source = sprintf('uint%s', bits{1});
    else
        % guess number of bits
        source = 'uint64';

        % we can have files with values either 32 or 64 bits.
        ts32 = fread(fid, 6, 'uint32');
        if length(ts32) == 6 && (ts32(1) == ts32(3) && ts32(3) == ts32(5))
            source = 'uint64';
        else
            source = 'uint32';
        end
    end
    
    fseek(fid, curFilePos, 'bof');

    % Read timestamps, 32 bit long
    ts = fread(fid, inf, source);

    % Transform the timestamps to seconds
    ts = ts / 10000;
end
