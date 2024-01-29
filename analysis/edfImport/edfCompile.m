%%  edfImport library v1.0 
%  Alexander Pastukhov 
%  kobi.nat.uni-magdeburg.de/sasha
%  email: pastukhov.alexander@gmail.com
%
%  edfCompile
%  OS-sensitive script for library compilation. Plese, modify paths  to
%  suite your environment. 

%% Compiles edfImport library
if (ispc())
  %% Windows % MGL DEC 2021: edited to x64 paths
  %edfapiIncludesFolder= 'C:/Program Files/SR Research/EyeLink/EDF_Access_API/Example';
  edfapiIncludesFolder = 'C:/Program Files (x86)/SR Research/EyeLink/Includes/eyelink';
  %edfapiLibraryFolder= 'C:/Program Files/SR Research/EyeLink/EDF_Access_API/lib/win32';
  edfapiLibraryFolder= 'C:/Program Files (x86)/SR Research/EyeLink/libs/x64';
  eval(sprintf('mex -I''%s'' edfMexImport.cpp EDFFILE2.cpp ''%s/edfapi64.lib''', edfapiIncludesFolder, edfapiLibraryFolder));
%   eval(sprintf('mex -g -I''%s'' edfMexImport.cpp EDFFILE2.cpp ''%s/edfapi.lib''', edfapiIncludesFolder, edfapiLibraryFolder));
elseif (ismac())
  %% Mac OS
  edfapiIncludesFolder= '/Library/Frameworks/edfapi.framework/Headers/';
  edfapiLibraryFolder= '/Library/Frameworks/edfapi.framework/';
  eval(sprintf('mex -I''%s'' edfMexImport.cpp EDFFILE2.cpp -L''%s''  LDFLAGS=''\\$LDFLAGS -framework edfapi'' ', edfapiIncludesFolder, edfapiLibraryFolder));  
elseif (isunix())
  %% Linux/Unix
  mex edfMexImport.cpp EDFFile2.cpp -ledfapi -lm -lz -lrt; 
else
  fprintf('Sorry, have no clue what OS you''ve got!\n');
end;
  
fprintf('Library was compiled successfully!\n');