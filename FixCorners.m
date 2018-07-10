function u = FixCorners(uc,invV,nN,nE)
% ut: modal values
u=zeros(nN,nE); u(:,:)=uc; ut=invV*u; %ut=V\u;

%% Fix corner elements to avoid weakness of BC
nEy=sqrt(nE); CE=[1,nEy,nE-nEy+1,nE]; u(:,CE)=repmat(0.5*ut(1,CE),nN,1);