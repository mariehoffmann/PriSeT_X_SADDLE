function [H, keys_ref, values_ref] = test()
    S = [{'ACGTACGT'}; {'TACGTACGT'}];
    mask = logical([1; 1]);
    dimer_range = 4:8;
    H = init_hash(S, mask, dimer_range);
    keys_ref = [{'ACGT'} {'ACGTA'} {'ACGTAC'} {'ACGTACG'} {'ACGTACGT'} ...
        {'CGTA'} {'CGTAC'} {'CGTACG'} {'CGTACGT'} {'GTAC'} {'GTACG'} ... 
        {'GTACGT'} {'TACG'} {'TACGT'} {'TACGTA'} {'TACGTAC'} {'TACGTACG'}];
    values_ref = [{[2.4]} {[0.5]} {[0.6667]} {[1]} {[2]} {[0.5]} ...
        {[0.6667]} {[1]} {[2]} {[0.6667]} {[1]} {[2]} {[1.1667]} ...
        {[2.2]} {[0.25]} {[0.3333]} {[0.5]}];
    
    if numel(keys_ref) ~= numel(H.keys)
        error('Unexpected number of keys in hash H')
    end

    index = cellfun(@strcmp, keys_ref, H.keys);
    if sum(index) ~= numel(keys_ref)
        error('Check generated keys in H')
    end

    nearEqual = @(x, y) abs(x - y) < .001;
    index = cellfun(nearEqual, values_ref, H.values);
    if sum(index) ~= numel(values_ref)
        error('Check computed values in H')
    end

    % badness B map
    B_id = zeros(numel(S), 1);
    ref1 = [115.2711 149.3333 597.3333 2048 4096];
    if ~nearEqual(ref1(1), badness(S, 4, 1, 1))
        error('expected 115.2711 for badness(S, [4], 1, 1)')
    end
    if ~nearEqual(ref1(2), badness(S, 5, 1, 1))
        error('expected 149.3333 for badness(S, [5], 1, 1)')
    end
    if ~nearEqual(ref1(3), badness(S, 6, 1, 1))
        error('expected 597.3333 for badness(S, [6], 1, 1)')
    end
    if ~nearEqual(ref1(4), badness(S, 7, 1, 1))
        error('expected 2048 for badness(S, [7], 1, 1)')
    end
    if ~nearEqual(ref1(5), badness(S, 8, 1, 1))
        error('expected 4096 for badness(S, [7], 1, 1)')
    end

    if ~nearEqual(sum(ref1), badness(S, dimer_range, 1, 1))
        error('expected %d for badness(S, [7], 1, 1)', sum(ref1))
    end
    

    B_id(1) = badness(S, dimer_range, 1, 1);
    B_id(2) = badness(S, dimer_range, 2, 2);
    
    

end

