function Score = get1_IGD(PopObj,PF)
Distance = min(pdist2(PF,PopObj),[],2);
Score    = mean(Distance);
end