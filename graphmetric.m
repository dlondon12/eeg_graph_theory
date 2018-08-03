function metric=graphmetric(adj_cube)
num_nodes=size(adj_cube,1);
num_files=size(adj_cube,3);
metric=zeros(num_nodes,5,num_files);
metric(:,1,:)=degrees_und(adj_cube);
for i=1:num_files
metric(:,2,i)=betweenness_wei(adj_cube(:,:,i))/((num_nodes-1)*(num_nodes-2));
metric(:,3,i)=clustering_coef_bu(adj_cube(:,:,i));
end
metric(:,4,:)=metric(:,1,:).*metric(:,2,:);
metric(:,5,:)=metric(:,1,:).*metric(:,2,:)./metric(:,3,:);