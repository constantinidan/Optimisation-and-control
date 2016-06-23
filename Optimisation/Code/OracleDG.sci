function [F,G] = OracleDG(lambda,ind)
    if ind == 2 then
        temp = Ar'*pr + Ad'*lambda
        ql = - sign(temp) .* sqrt(abs((1 ./ r) .*(temp)))
        F = -(-(1/3)*ql'*(temp) + pr'*Ar*ql + lambda'*(Ad*ql - fd))
        G = 0
    elseif ind == 3 then
        temp = -((Ar'*pr)+(Ad'*lambda))./r
        q = sqrt(abs(temp)).*sign(temp)
        G = Ad*q - fd
        F = 0
    elseif ind == 4 then
        temp = Ar'*pr + Ad'*lambda
        ql =   - sign(temp) .* sqrt(abs((1 ./ r) .*(temp)))
        F = -(-(1/3)*ql'*(temp) + pr'*Ar*ql + lambda'*(Ad*ql - fd))
        G =-(Ad*ql - fd)
    end
endfunction
