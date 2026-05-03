%% Setup
ops = sdpsettings('verbose', 0);  % <-- suppress output
basePath = "Petz recovery map\final data\";
nFiles=4;
s = sdpvar(4,1);
rho = [s(1), s(2)+1i*s(3); s(2)-1i*s(3), s(4)];
F_constraints = [rho>=0, trace(rho)==1];
% Clear results file at start of each run
fid = fopen("results.txt", "w"); fclose(fid);

%% Channel base info: {subfolder, state sub folder channel tag in filename}
channelInfo = {'AD', "\ad-petz(","ad"; 
    'PD', "\pd-petz(","pd" };

%% States per channel: {density matrix, label}
stateMap = struct();

stateMap.AD = { [1 0;0 0], "0", "|0>" ,"zero" ;[0 0;0 1], "1", "|1>", "one"; [0.5 0.5;0.5 0.5], "+x", "|+>", "plus";
    [0.859063 -1i*0.347956;1i*0.347956 0.140937], "0.93", "0.93|0> + 0.37|1>", "extra"};

stateMap.PD = {[0.5 0.5;0.5 0.5], "+x" , "|+>", "plus"; [0.5 -0.5;-0.5 0.5], "-x", "|->","minus"; [1 0;0 0], "0", "|0>", "zero"; 
    [0.859063 -1i*0.347956;1i*0.347956 0.140937], "0.93" , "0.93|0> + 0.37|1>", "extra"};

%% Count total iterations
totalIter = 0;
for c = 1:size(channelInfo, 1)
    totalIter = totalIter + size(stateMap.(channelInfo{c,1}), 1) * 44 * nFiles;
end
doneIter = 0;
 %% Waitbar
 wb = waitbar(0, 'Starting...', 'Name', 'Petz Recovery Fidelity');


%% Main loop
for c = 1:size(channelInfo, 1)
    chanName  = channelInfo{c, 1};
    chanPre   = channelInfo{c, 2};
    chanTag   = channelInfo{c, 3};
    chanStates = stateMap.(chanName);

    for st = 1:size(chanStates, 1)
        the       = chanStates{st, 1};
        statelabel = chanStates{st, 2};
        stateName = chanStates{st, 3};
        sub_folder = chanStates{st, 4};

        % Build dynamic filename: e.g. sub_folder\ad-petz(0,ad,0.2,0.5,0.8)_real1
        fileBase = basePath + chanName + "\" + sub_folder + chanPre  + statelabel + "," + chanTag + ",0.2,0.5,0.8)_";
        fidelity_all   = zeros(44, nFiles);

        for f = 1:nFiles
            % Import and sum integrals
            [r57,r13,r68,r24] = import_integral_3qubit(fileBase+"real"+f, 16, 37+22*4);
            [i57,i13,i68,i24] = import_integral_3qubit(fileBase+"imag"+f, 16, 37+22*4);
            Freal = r57 + r13 + r68 + r24;
            Fimag = i57 + i13 + i68 + i24;

            % Build B matrix
            idx = 2*(1:44)-1;
            B = [Freal(idx), Freal(idx+1), Fimag(idx), Fimag(idx+1)];

            for i = 1:44
                % Pick A based on channel, state and epsilon column
                eps_col = floor((i-1)/11) + 1;  % 1=p, 2=0.2, 3=0.5, 4=0.8

                if strcmp(chanName, 'AD') && statelabel == "0" && eps_col == 2
                    A = [0 1 0 0;0 1 0 0;0 0 1 0;0.5 0 0 -0.5];         %x(pi/2) tomo pulse
                else
                    A =  [0 1 0 0; 0.5 0 0 -0.5; 0 0 1 0; 0 0 1 0];     %y(pi/2) tomo pulse
                end

                optimize([rho>=0, trace(rho)==1], norm(A*s - B(i,:)')^2, ops);
                sol = value(s);
                output = [sol(1), sol(2)+1i*sol(3); sol(2)-1i*sol(3), sol(4)];
                fidelity_all(i,f) = real(trace(sqrtm(sqrtm(the)*output*sqrtm(the)))^2);
                % Update waitbar
                doneIter = doneIter + 1;
                waitbar(doneIter/totalIter, wb,sprintf('Channel: %s | State: %s | File: %d/%d | Point: %d/44', chanName, stateName, f, nFiles, i));
                
            end
        end
        % Mean and std across files
        fidelity_mean = mean(fidelity_all, 2);
        fidelity_std  = std(fidelity_all,  0, 2);

        % Reshape to 11x4
        Fmat_mean = reshape(fidelity_mean, 11, 4);
        Fmat_std  = reshape(fidelity_std,  11, 4);

        fid = fopen("results.txt", "a");
        fprintf(fid, '\nChannel: %s | State: %s\n', chanName, stateName);
        fprintf(fid, '--------------------------------------------------------------------------------\n');
        fprintf(fid, 'damping strength |          Petz recovery map with epsilon             \n');
        fprintf(fid, '     p           |       0.2              0.5              0.8        \n');
        fprintf(fid, '--------------------------------------------------------------------------------\n');
        for i = 1:11
            fprintf(fid, ' %7.4f±%6.4f |  %6.4f±%6.4f   %6.4f±%6.4f   %6.4f±%6.4f\n', ...
                Fmat_mean(i,1), Fmat_std(i,2),Fmat_mean(i,2), Fmat_std(i,2), Fmat_mean(i,3), Fmat_std(i,3), ...
                Fmat_mean(i,4), Fmat_std(i,4));
        end
        fprintf(fid, '--------------------------------------------------------------------------------\n');
        fclose(fid);
    end
end

