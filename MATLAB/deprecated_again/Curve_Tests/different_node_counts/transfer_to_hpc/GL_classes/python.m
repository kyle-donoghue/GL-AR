classdef python
    
    methods (Static=true)
        function imports
            pathname = strcat(fileparts(which('python')),'\python');
            if count(py.sys.path,pathname) == 0
                insert(py.sys.path,int32(0),pathname);
            end
            lib = py.importlib.import_module('HH_brian2');
            py.importlib.reload(lib);
            py.importlib.import_module('brian2');
            py.importlib.import_module('numpy');
        end
        function [t, v, I_ext] = simulate_HH_network(A,stim,duration)
            N = size(A,1);
            output = py.HH_brian2.simulate_HH_network(py.numpy.array(A),py.numpy.array(stim),duration);
            output = double(output);
            t = output(1,:);
            v = output(2:N+1,:);
            I_ext = output(N+2:end,:);
        end
    end
end

