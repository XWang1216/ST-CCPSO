function [x,v] = updatexv(x,v,w,c1,c2,fi,pBest,gBest,func_dim,xrange,subswarm_index,Obj_fitness,iter,Max_iter)
if iter<0.8*Max_iter
    [x,v] = CSO(x,v,w,c1,c2,fi,pBest,gBest,func_dim,xrange,subswarm_index,Obj_fitness);
else
    [x,v] =  pso(x,v,w,c1,c2,pBest,gBest,func_dim,xrange);
end