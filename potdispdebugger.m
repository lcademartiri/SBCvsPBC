[~, cmdout] = system('wmic cpu get L2CacheSize, L3CacheSize /format:list');
pattern = 'L3CacheSize=(\d+)';

% --- 2. Run the Extraction ---
% 'tokens': Returns the content captured within parentheses
% 'once': Stops after the first match
tokens = regexp(cmdout, pattern, 'tokens', 'once');

if isempty(tokens)
    warning('L3 Cache Size pattern not found in the command output.');
    L3_Cache_KB = NaN;
else
    % 3. Convert the captured string to a number
    L3_Cache_KB = str2double(tokens{1,1});
    
    fprintf('Extracted L3 Cache Size: %d KB\n', L3_Cache_KB);
end
L3_Cache_MB=L3_Cache_KB/1024;
for nm=1:14
    p=rand(2^nm,3).*(2*S.br)-S.br;
    N(nm,1)=size(p,1);
    for i0=1:5
        if nm<15
            tic
            [disppot2,coll1,coll2,dists]=potential_displacements_v2(p, S, H, H_interpolant, 0);
            times2(nm,i0)=toc;
    
            tic
            [disppot3,coll1,coll2,dists]=potential_displacements_v3(p, S, H, H_interpolant, 0);
            times3(nm,i0)=toc;
    
            tic
            [disppot4,coll1,coll2,dists]=potential_displacements_v4(p, S, H, H_interpolant, 0);
            times4(nm,i0)=toc;
        end

        % tic
        % [disppot5,coll1,coll2,dists]=potential_displacements_v5(p, S, H, H_interpolant, 0);
        % times5(nm,i0)=toc;
        % 
        % tic
        % [disppot6,coll1,coll2,dists]=potential_displacements_v6(p, S, H, H_interpolant);
        % times6(nm,i0)=toc;

        tic
        [disppot7,coll1,coll2,dists]=potential_displacements_v7(p, S, H, H_interpolant,0,L3_Cache_MB);
        times7(nm,i0)=toc;
        disp([nm,i0]);
    end
end