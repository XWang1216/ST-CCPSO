function [Score] = get_SP(PopObj)
% [M,~] = size(A);
% min_di = [];
% A1 = [];
% d = [];
% for i = 1:M
%     A1 = A(i,:)';
%     d = dist(A,repmat(A1,1,M));
%     d_raw = d(:,1);
%     [~,dmin_index] =  sort(d_raw);
%     min_di(i) = d_raw(dmin_index(2));
% end
% dmean = mean(min_di);
% for i = 1:M
%     min_di(i) = (min_di(i)-dmean)^2;
% end
% sp = sqrt(sum(min_di)/(M-1));
    Distance = pdist2(PopObj,PopObj,'cityblock');
    Distance(logical(eye(size(Distance,1)))) = inf;
    Score    = std(min(Distance,[],2));
end