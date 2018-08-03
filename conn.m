function [conn_cube,conncompscell]=conn(adj_cube)

for k=1:size(adj_cube,4)
for j=1:size(adj_cube,3)
    %[comps,comp_sizes] = get_components(adj_cube(:,:,j)); %labels the groups of all the components and figures out sizes of each group
    [~,comps]=graphconncomp(sparse(adj_cube(:,:,j,k)),'Directed',1','Weak',0);
    [connind,connsize]=mode(comps);
    %[connsize,connind]=max(comp_sizes);
    if connsize>1
        conncomps=(comps==connind);
        conn_mat=adj_cube(:,:,j,k);
        conn_mat(:,~conncomps)=0;
        conn_mat(~conncomps,:)=0;
        conn_mat=conn_mat(conncomps,conncomps);
        conn_cube{j,k}=conn_mat;
        conncompscell{j,k}=conncomps;
    else
        conn_cube{j,k}=[];
        conncompcell{j,k}=[];
    end
end
end



end