%% nSTAT_Install.m
% Author: Iahn Cajigas
% This script adds all of the appropriate folders to the matlab path and
% makes the nSTAT help files searchable within the MATLAB help browser.

%
% nSTAT v1 Copyright (C) 2012 Masschusetts Institute of Technology
% Cajigas, I, Malik, WQ, Brown, EN
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as published 
% by the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU General Public License for more details.
%  
% You should have received a copy of the GNU General Public License 
% along with this program; if not, write to the Free Software Foundation, 
% Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

fileLocation = which('nSTAT_Install'); 
index = strfind(fileLocation,'nSTAT_Install.m')-1;
rootDir =fileLocation(1:index);
helpDir = strcat(rootDir,'helpfiles');
display('Adding nSTAT help files to the matlab search database');
builddocsearchdb(helpDir); %make the toolbox files searchable

display('Adding nSTAT to the top of the search path');
addpath(genpath(rootDir),'-begin');
display('Saving path');
savepath;
clear fileLocation index rootDir helpDir dataDir libraries;
