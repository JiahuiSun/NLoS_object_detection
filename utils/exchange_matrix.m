function E = exchange_matrix(dim)
E=zeros(dim,dim);
for i=1:dim
    E(i,dim-i+1)=1;
end
end