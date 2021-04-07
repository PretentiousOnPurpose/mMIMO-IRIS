function lssb = GetIndexOfSSBurst(cases, fc)

switch upper(cases)
    case 'A'
        if fc <= 3e9
            n = [0,1];
            lssb = [2; 8] + 14*n;
            lssb = lssb(:).' + 1;
        elseif fc > 3e9 && fc <= 6e9
            n = [0,1,2,3];
            lssb = [2; 8] + 14*n;
            lssb = lssb(:).' + 1;
        end
        
    case 'B'
        if fc <= 3e9
            n = 0;
            lssb = [4; 8; 16; 20] + 28*n;
            lssb = lssb(:).' + 1;
        elseif fc > 3e9 && fc <= 6e9
            n = [0, 1];
            lssb = [4; 8; 16; 20] + 28*n;
            lssb = lssb(:).' + 1;
        end
        
    case 'C'
        if fc <= 3e9
            n = [0, 1];
            lssb = [2; 8] + 14*n;
            lssb = lssb(:).' + 1;
        elseif fc > 3e9 && fc <= 6e9
            n = [0, 1, 2, 3];
            lssb = [2; 8] + 14*n;
            lssb = lssb(:).' + 1;
        end
        
    case 'D'
        if fc > 6e9
            n = [0:3, 5:8, 10:13, 15:18];
            lssb = [4; 8; 16; 20] + 28*n;
            lssb = lssb(:).' + 1;
        end
        
    case 'E'
        if fc > 6e9
            n = [0:3, 5:8];
            lssb = [8; 12; 16; 20; 32; 36; 40; 44] + 56*n;
            lssb = lssb(:).' + 1;
        end
end


