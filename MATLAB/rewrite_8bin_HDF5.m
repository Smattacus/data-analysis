function rewrite_8bin_HDF5(filename, bool_8bits, compression)
%This is a utility function in order to rewrite a file from the old binary 
%format to the new HDF5 format. 
%
% rewrite_8bin_HDF5(filename, bool_8bits, compression)
%
% This program will take the old file and create a new one with the .h5
% extension, but with the same filename. Additionally, all the header
% values will be written to the new .h5 file as attributes in the same that
%the data are placed in.
%
%
%INPUTS
%filename       - Filename of the old data file.
%bool_8bits     - Boolean whether the input data file is in 32bit or 8bit
%                   format.
%compression    - Set how aggressive compression is. 9 - smallest, but most
%                   time. 1 is largest, but least time. Default is 6.
%
if nargin < 3
    compression = 6;
end
[data, H] = getPMTData(filename, bool_8bits);
newfn = strrep(filename, '.DAT', '.h5');
if bool_8bits == true
    dtype = 'uint8';
else
    dtype = 'uint32';
end
h5create(newfn, '/PMT_DATA_8BIT', [size(data, 1), size(data,2)], ...
    'Datatype', dtype, 'ChunkSize', [32, 4096], 'Deflate', compression);
h5write(newfn, '/PMT_DATA_8BIT', data);

for i=2:(size(H,1)-1)
    [temp, b] = textscan(H(i,:), '%s = %s');
    name = temp{1}; name = name{1};
    val = temp{2}; val = val{1};
    h5writeatt(newfn, '/PMT_DATA_8BIT', name, val);
end
end