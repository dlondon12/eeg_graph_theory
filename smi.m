%this program calculates the two different small world indices - omega (Lrand/L - C/Clatt) and
%sigma (C/Crand / L/Lrand) for a network with "n_nodes" nodes, "n_edg"
%edges
function [smlwrld,data_CL] = smi(conn_cube,n)


if iscell(conn_cube)==0
    Crand=zeros(size(conn_cube,3),1);
    Lrand=zeros(size(conn_cube,3),1);
    Clatt=zeros(size(conn_cube,3),1);
    data_CL=zeros(size(conn_cube,3),2);
    conn_cube=bsxfun(@times,conn_cube,~eye(size(conn_cube,1))); %remove 0s from diagonal
    
    for j=1:size(conn_cube,3)
        
        [comps,comp_sizes] = get_components(conn_cube(:,:,j));
        %n_nodes=max(comp_sizes);
        n_nodes=size(conn_cube,1);
        n_edg=nnz(conn_cube(:,:,j))/2;
        %n_edg=floor(n_edg);
        
        %calculate C & L for 100 random networks and average it to make
        %Crand and Lrand
        randnet_CL=ones(100,2);
        
        
        for i=1:n
            [CIJ] = makerandCIJ_und(n_nodes,n_edg);
            C = clustering_coef_bu(CIJ);
            randnet_CL(i,1)=mean(C);
            D=distance_bin(CIJ);
            randnet_CL(i,2)=charpath(D);
        end
        Crand(j)=mean(randnet_CL(:,1));
        Lrand(j)=mean(randnet_CL(:,2));
        
        %calculate C for 100 lattice networks
        latnet_C=ones(100,1);
        for i=1:n
            [CIJ]=makelatticeCIJ(n_nodes,n_edg);
            C=clustering_coef_bu(CIJ);
            latnet_C(i)=mean(C);
        end
        Clatt(j)=mean(latnet_C);
        
        C=clustering_coef_bu(conn_cube(:,:,j)); %calculate clustering coefficients for all the nodes
        data_CL(j,1)=sum(C)/n_nodes; %record the average clustering coefficient
        D=distance_bin(conn_cube(:,:,j));
        data_CL(j,2)=charpath(D);
        
    end
    
else
    Crand=zeros(length(conn_cube),1);
    Lrand=zeros(length(conn_cube),1);
    Clatt=zeros(length(conn_cube),1);
    data_CL=zeros(length(conn_cube),2);
    for j=1:length(conn_cube)
        conn_mat=conn_cube{j};
        conn_mat(logical(eye(size(conn_mat,1))))=0; %remove 0s from diagonal
        
        [comps,comp_sizes] = get_components(conn_mat);
        %n_nodes=max(comp_sizes);
        n_nodes=size(conn_mat,1);
        n_edg=nnz(conn_mat)/2;
        %n_edg=floor(n_edg);
        
        %calculate C & L for 100 random networks and average it to make
        %Crand and Lrand
        randnet_CL=ones(100,2);
        
        
        for i=1:n
            [CIJ] = makerandCIJ_und(n_nodes,n_edg);
            C = clustering_coef_bu(CIJ);
            randnet_CL(i,1)=mean(C);
            D=distance_bin(CIJ);
            randnet_CL(i,2)=charpath(D);
        end
        Crand(j)=mean(randnet_CL(:,1));
        Lrand(j)=mean(randnet_CL(:,2));
        
        %calculate C for 100 lattice networks
        latnet_C=ones(100,1);
        for i=1:n
            [CIJ]=makelatticeCIJ(n_nodes,n_edg);
            C=clustering_coef_bu(CIJ);
            latnet_C(i)=mean(C);
        end
        Clatt(j)=mean(latnet_C);
        
        C=clustering_coef_bu(conn_mat); %calculate clustering coefficients for all the nodes
        data_CL(j,1)=sum(C)/n_nodes; %record the average clustering coefficient
        D=distance_bin(conn_mat);
        data_CL(j,2)=charpath(D);
        
    end
end
smlwrld(:,1)=(data_CL(:,1)./Crand)./(data_CL(:,2)./Lrand); %sigma
smlwrld(:,2)=(Lrand./data_CL(:,2))-(data_CL(:,1)./Clatt); %omega

data_CL(:,1)=(data_CL(:,1)./Crand);
data_CL(:,2)=(data_CL(:,2)./Lrand);
end