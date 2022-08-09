function H = update_hash(H, S, dimer_range, i_del, i_ins)
%Remove or add all patterns from indicated sequence from hash. 
%In case only one sequence is added or removed, set the other index to -1.
    arguments
        H {mustBeNonempty}
        S {mustBeVector}
        dimer_range (1, :) {mustBeNonempty, mustBeVector, mustBeInteger}
        i_del (1,1) {mustBeInteger}  % set to -1 to ignore
        i_ins (1,1) {mustBeInteger}  % set to -1 to ignore
    end

    if i_del == i_ins
        return
    end

    if i_del > 0 && H.Count == 0
        error('Cannot delete from empty map!');
    end
    
    nearZero = @(x) abs(x) < .00001;
    
    if i_del > 0
        s = S{i_del};
        for k = dimer_range
            for pos = 1:length(s) - k + 1
                kmer = s(pos:pos+k-1);   
                if ~isKey(H, kmer)
                    error('kmer not found in cash: %s from %s', kmer, s);
                end
                d = length(s) - (k + pos) + 1;
                H(kmer) = H(kmer) - 1 / (d + 1);
                if nearZero(H(kmer))
                    remove(H, kmer);
                end
            end
        end
    end
    if i_ins > 0
        s = S{i_ins};
        for k = dimer_range
            for pos = 1:length(s) - k + 1
                kmer = s(pos:pos+k-1);   
                if ~isKey(H, kmer)
                    H(kmer) = 0;
                end
                d = length(s) - (k + pos) + 1;
                H(kmer) = H(kmer) + 1 / (d + 1); 
            end
        end
    end
end