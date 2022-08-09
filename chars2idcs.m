function idcs = chars2idcs(cs)
%Encode char array or cell of char array to index array 
% encoding scheme {A: 1, C: 2, G: 3, T: 4}
    c2i = zeros(256,1);
    c2i(double('A')) = 1;
    c2i(double('C')) = 2;
    c2i(double('G')) = 3;
    c2i(double('T')) = 4;
    sprintf('chars2idcs, type input = %s', class(cs));
    if iscell(cs)
        cs = cs{1,1};
    end
    idcs = c2i(cs); 
end