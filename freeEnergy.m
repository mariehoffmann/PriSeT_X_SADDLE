% see http://www.premierbiosoft.com/netprimer/netprlaunch/Help/Theories_and_Formulas.htm

function dG = freeEnergy(s, dH_table, dS_table, argv)
%freeEnergy compute free energy of an oligonucleotide using NN method [2].
%unit of return value is kcal/mol
%   s   index or char array representing oligonucleotide
%   dH_table  enthalpy table loaded from enthalpy.csv
%   dS_table  entropy table loaded from entropy.csv
    if ischar(s)
        s = chars2idcs(s);
    end

    % enthalpy for helix formation
    dH = enthalpy(s, dH_table);
    
    % entropy for helix formation
    dS = entropy(s, dS_table);
    
    % Tm based on nearest neighbour method using thermodynamics [3]
   
%     T = Tm(s, "nearest_neighbour", "Kelvin", argv)
    dG = dH - argv('DG_temp') * dS;
    dG = dG / 1000; %  
end

function t = Tm(s, formula, unit, option)
%melting temperature in degree "Celsius" or "Kelvin" given formula
arguments 
    s (:, 1)
    formula (1,1) string
    unit (1,1) string
    % option {mustBeMember(option,["linear","cubic"])}
    option {mustBeNonempty} = containers.Map()
end
    if formula == "wallace"
        t = wallace(s);
    elseif formula == "nearest_neighbour"
        argv = option;
        K_con = argv('mvi') + 4 * sqrt(argv('Mg') * 1000); % total Na+ equivalent
        t = argv('dH') / (argv('dS') + argv('R') * log(argv('C') / 4)) + 16.6 * log10(K_con) / (1 + .7 * K_con);
    else 
        error("Unknown Tm method");
    end
    if unit == "Kelvin"
        t = t - 273.15;
    end
end

function t = wallace(s)
%wallace melting temperature according to Wallace rule
    AT = 0;
    GC = 0;
    if iscell(s)
        error('Expected input to be array of chars|ints, but not cell.')
    end
    if iscell(s) && ischar(s{1,1})
        s = cell2mat(s);
        AT = length(find(s == 'A' | s == 'T'));
        GC = length(find(s == 'C' | s == 'G'));
    elseif isnumeric(s)
        AT = length(find(s == 1 | s == 4));
        GC = length(find(s == 2 | s == 3));
    else
        error('Unexpected input type: %s', class(s));
    end
    t = AT * 2 + GC * 4;
end
