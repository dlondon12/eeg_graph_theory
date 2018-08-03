function r=synch(conn_cube)

if iscell(conn_cube)==0
    conn_cube=bsxfun(@times,conn_cube,~eye(size(conn_cube,1))); %remove 0s from diagonal
    eigen_val=zeros(2,size(conn_cube,3));
    dg=squeeze(sum(conn_cube,1));
    if size(conn_cube,3)==1
        dg=dg';
    end
    dgmat=zeros(size(conn_cube));
    for i=1:size(conn_cube,3)
        dgmat(:,:,i)=diag(dg(:,i));
    end
    conn_cube=dgmat-conn_cube;
    for i=1:size(conn_cube,3)
        lamd=sort(eig(conn_cube(:,:,i)),'ascend');
        eigen_val(1,i)=lamd(2);
        eigen_val(2,i)=lamd(end);
    end
    
else
    eigen_val=zeros(2,length(conn_cube));
    for j=1:length(conn_cube)
        conn_mat=conn_cube{j};
        if ~isempty(conn_mat)
        
        conn_mat(logical(eye(size(conn_mat,1))))=0; %remove 0s from diagonal
        
        dg=squeeze(sum(conn_mat,1));
        dg=dg';
        dgmat=diag(dg);
        conn_mat=dgmat-conn_mat;
        lamd=sort(eig(conn_mat),'ascend');
        eigen_val(1,j)=lamd(2);
        eigen_val(2,j)=lamd(end);
        end
    end
    
end

r=eigen_val(2,:)./eigen_val(1,:);

end