function v = Tp(T,x,p)
    %%% Input: tensor T, vector x, paramter p
    %%% Output: action of the operator T_p on x
    
    Mux = meanp(x,x,p);
    spMux = sptensor(Mux);
    v = ttt(T,spMux,[2 3], [1 2]);
    v = double(v);
           
end

function m = meanp(a,b,p)
    if p == 2
        u = abs(a).^2;
        v = abs(b).^2;
        m = sqrt((u + v') ./ 2);
    elseif p == 1
        m = (a./2) + (b./2)' ;
    elseif p == 0
        m = sqrt(a * b');
    else
        u = abs(a).^p;
        v = abs(b).^p;
        m = ((u + v') ./ 2).^(1/p);
    end
end


         