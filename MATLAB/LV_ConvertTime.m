function ts = LV_ConvertTime(Date1, time1, Date2, time2)
%Helper function to convert a LABVIEW file format date and time to seconds.
%
%In a labview text file, the date and time are formatted like:
%
%    'Date	YYYY/MM/DD'
%    'Time	17:59:18.0940682100351102877'
%
%This function takes these two strings, parses them, and returns the number
%of seconds elapsed since the original time.
%
%INPUTS:
%   Date1   - Date string.
%   Time1   - Time string.
%   Date2   - Reference Date String.
%   Time2   - Reference Time String.
    elaps_days = convertDateToDecDays(Date2) - convertDateToDecDays(Date1);
    elaps_seconds = convertTimeToDecSeconds(time2) - convertTimeToDecSeconds(time1);
    if elaps_seconds < 0
        ts = elaps_days * 24 * 3600 - elaps_seconds;
    else
        ts = elaps_days * 24 * 3600 + elaps_seconds;
    end
end

function decs = convertTimeToDecSeconds(time)
%This function converts the labview time string to decimal seconds.    
% The string is formatted like
%'Time	17:59:18.0940682100351102877'
    decs = str2num(t(6:7)) / 3600 + str2num(t(9:10)) / 60 + str2num(t(12:33));
end

function decd = convertDateToDecDays(date)
%This function gives the number of days since Jan 1, 2010, using the
%labview interface string.
%
%   INPUT:
%       date    - Date string.
%   OUTPUT:
%       decd    - Decimal days since Jan 1 1970.
    y = str2num(date(6:9));
    m = str2num(date(11:12));
    d = str2num(date(14:15));
    eyears = y - 2010;
    leaps = floor(eyears/4);
    switch(m)
        case 1
            edays = d;
        case 2
            edays = 31 + d;
        case 3
            edays = 59 + d;
        case 4
            edays = 90 + d;
        case 5
            edays = 120 + d;
        case 6
            edays = 151 + d;
        case 7
            edays = 181 + d;
        case 8 
            edays = 212 + d;
        case 9
            edays = 243 + d;
        case 10
            edays = 273 + d;
        case 11
            edays = 304 + d;
        case 12
            edays = 334 + d;
        otherwise
            display('Not a month between 1 and 12. Check input?');
    end
    decd = eyears * 365 + edays;
end