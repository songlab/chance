function d=dist_fun(Xi,Xj)

for i=1:size(Xj,1)
    d(i)=max(abs(Xi/Xi(end)-Xj(i,:)/Xj(i,end)));
end