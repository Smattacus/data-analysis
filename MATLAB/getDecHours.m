function hrs = getDecHours(filename)
%Function which reads a filename and gives the decimal hour the file was
%acquired.
%
%hrs = getDecHours(filename)
%
%To be used on PMT h5 data.
time = filename(1:20);
hrs = str2double(time(12:13)) + str2double(time(15:16)) / 60 + str2double(time(18:19)) / 3600;