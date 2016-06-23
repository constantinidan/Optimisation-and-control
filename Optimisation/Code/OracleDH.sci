function [F,G,H] = OracleDH(lambda,ind)
    if ind == 2 | ind == 3 | ind == 4 then
        [F,G] = OracleDG(lambda,ind)
        H = 0
    elseif ind == 5 then
        temp = Ar'*pr + Ad'*lambda
        temp2 = 2*sqrt(r .* abs(temp))
        temp3 = ones(temp) ./ temp2
        H = Ad * diag(temp3) * Ad'
        F = 0
        G = 0
    elseif ind == 6 then
        [F,G] = OracleDG(lambda,3)
        temp = Ar'*pr + Ad'*lambda
        temp2 = 2*sqrt(r .* abs(temp))
        temp3 = ones(temp) ./ temp2
        H = Ad * diag(temp3) * Ad'
        
    elseif ind == 7 then
        [F,G] = OracleDG(lambda,4)
        temp = Ar'*pr + Ad'*lambda
        temp2 = 2*sqrt(r .* abs(temp))
        temp3 = ones(temp) ./ temp2
        H = Ad * diag(temp3) * Ad'
    end
endfunction
