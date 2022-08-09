function t = revcomp(s)
    arguments
        s {mustBeText}
    end
    t = strrep(strrep(strrep(s,'A','X'),'T','A'),'X','T');
    t = strrep(strrep(strrep(t,'C','X'),'G','C'),'X','G');
    t = reverse(t);
end