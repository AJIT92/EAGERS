%RUN: runs the EAGERS software
global Intitialized Model_dir
if isempty(Intitialized)
    Model_dir=strrep(which('RUN.m'),fullfile('main','RUN.m'),'');
    AddPaths()
    Intitialized=1;
end
K = menu('Select Option', 'Simple Dispatch Test', 'Launch EAGERS Interface');
if K==1
    TestDispatch
else
    EAGERS
end