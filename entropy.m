function e = entropy(s, dS_table)
%Compute free entropy [cal/mol]
    if size(s, 1) < size(s, 2)
        s = s';
    end
    if iscell(s)
        s = s{1, 1};
    end

    if ischar(s)
        s = chars2idcs(s);
    end
    e = 15.1;  % initialization value as in [2]
    for i = 1:numel(s)-1
        e = e + dS_table{s(i, 1), s(i+1, 1)};
    end
    e = -e;
end