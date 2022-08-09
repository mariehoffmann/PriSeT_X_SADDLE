function b = badness(S, dimer_range, i1, i2)
%Compute badness score of S_i1 vs S_i2, i.e., some likeliness of
%dimerization. Most correct implementation would start with longest common 
%substrings and cancelling of all shorter, contained substrings via an 
%updated bit mask (here: seen).
arguments
    S (:,1) {mustBeVector}
    dimer_range (1,:) {mustBeVector, mustBeInteger}
    i1 (1,1) {mustBeInteger}
    i2 (1,1) {mustBeInteger}
end
    b = .0;
    s1 = S{i1};
    s2 = revcomp(S{i2});
    if isstring(s1)
        s1 = convertStringsToChars(s1);
        s2 = convertStringsToChars(s1);
        disp('WARNING: got string, expected char array.')
    end
    % do not count shorter substrings - store flag for processed longer patterns
%     seen = zeros(length(s2), 1);
    % cumulative sum of GC content for fast computation
    GC_cs = cumsum(cat(2, [0], (s1 == 'C' | s1 == 'G'))); 
    for k = dimer_range %flip(dimer_range)
%         seen_new = zeros(length(s2), 1);
        
        for pos = 1:length(s1) - k + 1
            GC = GC_cs(pos + k) - GC_cs(pos);
            kmer = s1(pos:pos+k-1);
            sprintf('pos = %d, k = %d, GC = %d, kmer = %s', pos, k, GC, kmer);
            d1 = length(s1) - (pos + k) + 1;
            locs = strfind(s2, s1(pos:pos+k-1));
            for loc = locs
                % skip if longer substring previously seen
%                 if seen(loc) 
%                     continue
%                 end
                d2 = loc - 1 %length(s2) - (k + loc) + 1;
                sprintf('len = %i, GC = %i, d1 = %i, d2 = %i', k, GC, d1, d2);
                b = b + 2^k * 2^GC/((d1 + 1) * (d2 + 1));
%                 % mark all subsequent substrings as seen
%                 for offset = 0 : k-1
%                     seen_new(loc + offset) = 1;
%                 end
            end
        end
%         seen = seen | seen_new;
    end
end