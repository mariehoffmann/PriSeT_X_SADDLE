function h = enthalpy(s, dH_table)
%Compute free enthalpy [cal/˚K/mol]
%   s      cell or array of integers or chars
%   dH_table    table of pairwise enthalpy scores
    if size(s, 1) < size(s, 2)
        s = s';
    end
    if iscell(s)
        s = s{1, 1};
    end
    if ischar(s)
        s = chars2idcs(s);
    end
    h = .0;
    for i = 1:numel(s)-1
        h = h - dH_table{s(i, 1), s(i+1, 1)};
    end
end