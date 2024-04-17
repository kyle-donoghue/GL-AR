poolobj = gcp('nocreate'); % If pool, delete
if ~isempty(poolobj)
    delete(poolobj);
else
    parpool('local',32);
end
tic
