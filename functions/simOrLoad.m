function [Model,tit] = simOrLoad(name, runSimulation, simFcn, keepIdx, keepVars)
% simOrLoad  Run a model simulation and save the part a figure needs, or load the saved copy.
%
%   [Model,tit] = simOrLoad(name, runSimulation, simFcn, keepIdx, keepVars)
%
%   name           file name in simulation_results/ (without .mat)
%   runSimulation  true  -> call simFcn() (getScore), keep only the stimulation protocols keepIdx and the
%                           variables keepVars, and save them to simulation_results/<name>.mat
%                  false -> load simulation_results/<name>.mat (falls back to simulating if the file is missing)
%   simFcn         function handle returning [Model,tit], e.g. @() getScore(Ms,pp,tt,nj,dt,tspan,yin,Par,M)
%   keepIdx        indices of the protocols (Model{i}) used by the figure; other cells are left empty
%   keepVars       variable names kept from each protocol table (e.g. {'NfkB','A20'})
%
% A full 1000-cell result holds every protocol and variable (~1.6 GB); the saved files keep only what the
% figure plots, so the figures can be reproduced in seconds without re-running the simulations.

resdir = 'simulation_results';   % relative to the main code folder (the current folder)
f = fullfile(resdir, [name '.mat']);
if runSimulation || ~isfile(f)
    [Model, tit] = simFcn();
    for i = 1:numel(Model)
        if ismember(i, keepIdx)
            Model{i} = Model{i}(:, keepVars);
        else
            Model{i} = [];
        end
    end
    if ~isfolder(resdir), mkdir(resdir); end
    save(f, 'Model', 'tit', 'keepIdx', 'keepVars', '-v7');
    fprintf('simOrLoad: simulated and saved %s\n', f);
else
    S = load(f, 'Model', 'tit');
    Model = S.Model; tit = S.tit;
    fprintf('simOrLoad: loaded %s\n', f);
end
end
