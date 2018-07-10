function u = PPlimiter(uc,invV,nN,nE)
% ut: modal values
u=zeros(nN,nE); u(:,:)=uc; ut=invV*u; %ut=V\u;

% Find the smallest value
m=min(u); ut1=0.5*ut(1,:); 

%% CPUtime = 12secs
% Theta
%theta=zeros(size(ut1));
%for e=1:nE; theta(e)=min(1,ut1(e)./(ut1(e)-m(e))); end

% Apply positivity limiter
%for e=1:nE; for i=1:nN; u(i,e)=ut1(e)+theta(e).*(u(i,e)-ut1(e)); end; end

%% CPUtime = 4secs
% Theta (no for loop)
theta=min([ones(1,nE);ut1./(ut1-m)]);

% Apply positivity limiter
u=repmat(ut1,nN,1)+repmat(theta,nN,1).*(u-repmat(ut1,nN,1));

%% Fix corner elements to avoid weakness of BC
%nEy=sqrt(nE); CE=[1,nEy,nE-nEy+1,nE]; u(:,CE)=repmat(0.5*ut(1,CE),nN,1);