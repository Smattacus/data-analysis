function [A H] = read_GAGE(filename)
%Reads a GAGE Scope file into MATLAB.
%Outputs:
%A   - Data from GAGE file
%H   - Array of header lines
f = fopen(filename);
H = char([]);
while 1
    t = deblank(fread(f, 128, '*char').');
    H = strvcat(H, t);
    if strcmp(t, 'BEGIN DATA')
        break
    end
end
A = fread(f, '*float');
H = cellstr(H);
fclose(f);