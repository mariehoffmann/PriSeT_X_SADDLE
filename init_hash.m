function H = init_hash(S, mask, dimer_range)
%Accumulate 1/(1+d) per pattern for all (s_i, s_j) for i <= j
    arguments
        S (:,1) {mustBeVector}
        mask (:, 1) {mustBeNumericOrLogical}
        dimer_range (1, :) {mustBePositive}
    end
    if ~islogical(mask)
        mask = logical(mask);
    end
    H = containers.Map('KeyType', 'char', 'ValueType', 'double');
    idcs = 1:numel(S);
    idcs = idcs(mask);

    for i_pos = 1:numel(idcs)
        s1 = S{idcs(i_pos)};
        for k = dimer_range
            if numel(k) > 1
                error('Iterate index-wise over dimer range')
            end
            for pos = 1:length(s1) - k + 1
                kmer = s1(pos:pos+k-1);   
                if isKey(H, kmer)
                    continue
                end
%                 d1 = length(s1) - (k + pos) + 1;
                H(kmer) = 0; % 1 / (d1 + 1);
                for j_pos = i_pos:numel(idcs)
%                     sprintf('j_pos = %i', j_pos);
                    s2 = S{idcs(j_pos)};
% %                     sprintf('seq = %s, pos = %d, k = %d, kmer = %s', s2, pos, k, kmer);
                    locs = strfind(s2, kmer);
                    for loc = locs
                        d = length(s2) - (k + loc) + 1;
%                         sprintf('', )
                        H(kmer) = H(kmer) + 1 / (d + 1);
%                         sprintf('add %d, result = %d', 1 / (d + 1), H(kmer));
                    end
                end
                
            end
        end
    end
end
