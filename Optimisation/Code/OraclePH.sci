function [F,G,H] = OraclePH(qc,ind)
    if ind == 2 | ind == 3 | ind == 4 then
        [F,G] = OraclePG(qc,ind)
        H = 0
    elseif ind == 5 then
        H = 2*B'*diag(r.*(abs(q0+B*qc)))*B
        F = 0
        G = 0
    elseif ind == 6 then
        [F,G] = OraclePG(qc,3)

        H = 2*B'*diag(r.*(abs(q0+B*qc)))*B
    elseif ind == 7 then
        [F,G] = OraclePG(qc,4)

        H = 2*B'*diag(r.*(abs(q0+B*qc)))*B
    end
endfunction
