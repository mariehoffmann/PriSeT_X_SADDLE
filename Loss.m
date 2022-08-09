function loss = Loss(S, mask, dimer_range, H, B_id)
    arguments
        S (:, 1) cell
        mask (:, 1) {mustBeNumericOrLogical}
        dimer_range (1, :) {mustBeNumeric}
        H {mustBeNonempty}
        B_id {mustBeNonempty}
    end

    idcs = 1:numel(S);
    idcs = idcs(mask);

    % accumulate badness relating to self-dimerization
    l2 = .0;
    for i_pos = 1:numel(idcs)
        l2 = l2 + B_id(idcs(i_pos));
    end
 
    % iterate over all sequences of current set S_g, extract patterns penalize
    % when reverse complement present
    l1 = .0;
    for i_pos = 1:numel(idcs)
        s = S{idcs(i_pos)};
        GC_cs = cumsum(cat(2, [0], (s == 'C' | s == 'G'))); 
        for k = dimer_range
            for pos = 1:length(s) - k + 1
                kmer = s(pos:pos+k-1);
                kmer_rc = revcomp(kmer);
                if ~isKey(H, kmer_rc)
                    continue
                end
                GC = GC_cs(pos + k) - GC_cs(pos);
                d = length(s) - (k + pos) + 1;
                l1 = l1 + (2^k * 2^GC / (d + 1)) * H(kmer_rc);
            end
        end
    end
    loss = (l1 + l2) / 2;
end