function  Tracks = importXMLfile(filename)
%IMPORTFILE Import data from a text XML file
%  IMAGE000016TRACKS = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  Example:
%  Tracks = importXMLfile("V:\Inyoung\IJ2021_Flow_measurement\20211119ij_flow_measurement\fish3\Image 000016_Tracks.xml", [1, Inf]);
%
%  See also READTABLE.

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = '"';

% Specify column names and types
opts.VariableNames = ["detection time", "time in frame", "x position", "x in pixel", "y position", "y in pixel"];
opts.SelectedVariableNames = ["detection time", "time in frame", "x position", "x in pixel", "y position", "y in pixel"];
opts.VariableTypes = ["categorical", "double", "categorical", "double", "categorical", "double"];

% Specify file level properties
opts.PreserveVariableNames = true;
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Import the data
Tracks = readtable(filename, opts);

end