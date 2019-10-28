function [T] = build_triangles_tensor(A,varargin)
    addpath tensor_toolbox/
    
    par = inputParser;
    par.addParameter('type','standard');
    par.parse(varargin{:});
    normalization_type = par.Results.type;
    
    n = size(A,1);
    T = sptensor([n,n,n]); 
    v = sum(A,2); idv = 1:n; idv = idv(v>0); %B = A(idv,:);
    d = sum(A,1);
    
    if strcmp(normalization_type, 'random_walk')
        M = A .* A^2;
    elseif strcmp(normalization_type, 'local_closure')
        Ad = A*d(:);
    end
    
    
    for i = 1 : length(idv)
       r = A(idv(i),:);
       idr = 1:n;
       idr = idr(r>0);
             
       for j = 1 : length(idr)
           r2 = A(idr(j),:);
           
           idr2 = 1:n;
           idr2 = idr2(r2>0);
                      
           matches = intersect(idr,idr2);
           
           for m = 1 : length(matches)
               switch normalization_type
                   case 'standard'
                        T(idv(i), idr(j), matches(m)) = 1;
                   case 'pagerank'
                        T(idv(i), idr(j), matches(m)) = 1;
                   case 'watts_strogatz'
                        T(idv(i), idr(j), matches(m)) = 1./(d(idv(i)).*(d(idv(i))-1));
                   case 'local_closure'
                       T(idv(i), idr(j), matches(m)) = 1./(Ad(idv(i))-d(idv(i)));
                   case 'random_walk'
                        T(idv(i), idr(j), matches(m)) = 1 ./ M(idr(j), matches(m));
               end
               
           end
       end
    end
    
    if strcmp(normalization_type, 'pagerank')
        M = ttv(T,{ones(n,1)},[1]);
        vecM = tenmat(M,[1 2]);
        v = 1./double(vecM); v(isinf(v))=1;
        matT = tenmat(T, [1], [2 3]);
        matT2 = bsxfun(@times,double(matT)',v);
        T2 = tensor(double(matT2), [n n n]);
        T = sptensor(T2);
    end
end

