function idx = get_rand_index(mask, bit)
%Generate random index where mask(idx)==0 (bit=0) or mask(idx)==1 (bit=1)

    arguments
        mask (:, 1) {mustBeNumericOrLogical, mustBeNonempty}
        bit (1,1) {mustBeNumericOrLogical} 
    end

    n = length(mask);
    if bit ~= 0 && bit ~= 1
        error('Expected bit to be in [0, 1]')
    end
    if bit && sum(mask) == 0
        error('Expected at least one bit set in mask')
    end
    if ~bit && sum(mask) == n
        error('Expected at least one unset bit in mask')
    end

    idx = randi(length(mask), 1);
    if mask(idx) && ~bit
        while mask(idx)
            idx = max(1, mod(idx + 1, n+1));
        end
    elseif ~mask(idx) && bit
        while ~mask(idx)
            idx = max(1, mod(idx + 1, n+1));
        end
    end
end