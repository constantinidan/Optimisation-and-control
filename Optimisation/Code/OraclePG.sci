function [F,G] = OraclePG(qc,ind)
    if ind == 2 then
        temp = (q0 + B*qc)'*((r.*(q0 + B*qc)).*abs(q0 + B*qc))
        temp2 = pr'*(Ar*(q0 + B*qc))
        F = (1/3)*(temp) + temp2
        G = 0
    elseif ind == 3 then
        temp = ((r.*(q0 + B*qc)).*abs(q0 + B*qc))
        temp3 = (B')*temp
        temp2 = ((Ar*B)')*pr
        //disp(size(temp2))
        //disp(size(temp3))
        G = temp3 + temp2
        F = 0
    elseif ind == 4 then
        temp = (q0 + B*qc)'*((r.*(q0 + B*qc)).*abs(q0 + B*qc))
        temp2 = pr'*(Ar*(q0 + B*qc))
        F = (1/3)*(temp) + temp2
        temp3 = ((r.*(q0 + B*qc)).*abs(q0 + B*qc))
        temp4 = (B')*temp3
        temp5 = ((Ar*B)')*pr
        G = temp4 + temp5
    end
endfunction
