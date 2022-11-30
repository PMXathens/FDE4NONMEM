                                %----------------------------------------------%----------------------------------------------%  
                                                                    % Grunwald-Letnikov Scheme %
                                %----------------------------------------------%----------------------------------------------%                                                                           %%%%%  %%%%%



%----------------------------------------------%
         % Setting the Parameters %
%----------------------------------------------%

                  
alpha =1 ;                  % Fractional Order
T = 18;                     % Total Time
h = 0.05;                   % Step
N = T/h;                    % Number of Steps
k = 5;                      % Number of Corrections8
k01 = 0.3;                  % Parameter of E:uation
k2 = 0.6;                   % Parameter of E:uation
y0 = [5,0];                 % Initial Condition8
V = 1;
dim = 2;                    % Number of Equations
numd = 3;                   % Number of Doses
params = [k01,k2,V];
%----------------------------------------------%
          % Initializing Matrices %
%----------------------------------------------%

y = zeros(N,dim);           % Solution Matrix
t = zeros(1,N);             % Time Matrix
bi = zeros(1,N);            % Fractional Coefficient matrix
doses = [0,0,0];            % Doses Matrix
tdoses = [0,0,0];           % Dosing Times
temp = zeros(k,dim);
sum = zeros(1,dim);
psi = zeros(1,dim);
gf = zeros(k,dim);
yf = zeros(N+1,dim);
tf = zeros(N+1,1);
%----------------------------------------------%
          % Fractional Coefficients %
%----------------------------------------------%
omega0=1;
for i=1:N
   if i==1
    omega(i)=omega0*(1-(1-alpha)/(i));
   else
    omega(i)=omega(i-1)*(1-(1-alpha)/(i));
   end
end
for j=1:N
    t(1,j)=h*j;
end 

%----------------------------------------------%
       % Numerical Algorithm for Solver %
%----------------------------------------------%

for n=1:N
   
    sum=zeros(1,dim);
    if n>1
     for q=1:(numel(tdoses)-1) 
        if tdoses(q)<=t(n) && t(n)<=tdoses(q+1)
            temp(1,:)=y(n-1,:)+doses(q);                    
      elseif tdoses(numd)<=t(n)
            temp(1,:)=y(n-1,:)+doses(end);    
        end
     end
    y(n,:)= temp(1,:);                                     
    for j=1:n-1
    ffv =  ff(t(j),y(j,:),params);
    sum(:)=sum(:)+omega(n-j).*ffv(:);                       
    end
    else
    temp(1,:) = y0(:);
    y(n,:)=y0(:); 
    end
    for q=1:(numel(tdoses)-1) 
        if tdoses(q)<=t(n) && t(n)<=tdoses(q+1)
            psi(:)=y0(:)+h^alpha.*sum(:)+doses(q);          
        elseif tdoses(numd)<=t(n)
            psi(:)=y0(:)+h^alpha.*sum(:)+doses(end);    
        end
    end    

    for i=2:k
        
    for q=1:(numel(tdoses)-1) 
        if tdoses(q)<=t(n) && t(n)<=tdoses(q+1)
            temp(i,:) = y(n,:)+doses(q);                    
        elseif tdoses(numd)<=t(n)
            temp(i,:) = y(n,:)+doses(end);    
        end
    end 
        temp(i,:) = y(n,:);
        ffv(:) = ff(t(n),temp(i,:),params); 
        gf(i,:) = temp(i,:)-h^(alpha).*ffv(1,:)-psi(1,:);  
        ffvd = ffdot(t(n),temp(i,:),params);              
        gdf = eye(dim)-h^(alpha).*ffvd;
        y(n,:)=temp(i,:)-gf(i,:)/(gdf);                      
    end
end

tf(1) = 0;
tf(2:N+1) = t(1:N);
yf(2:N+1,:) = y(1:N,:);
yf(1,:) = y0 (:);
yf(2:end,:) = y(:,:);

%----------------------------------------------%
          % ODE15s for Comparison %
%----------------------------------------------%

sol=ode15s(@(t,y)rate(t,y,k01,k2),[0 T],y0)
ycl=deval(sol,tf);


%----------------------------------------------%
                   % Plots %
%----------------------------------------------%

 figure(1)
 plot(tf,yf,'o','linewidth',2)
 hold on;
 plot(tf,ycl/V,'linewidth',2)
 hold on;
 xlabel('t','Interpreter','latex')
 ylabel('y','Interpreter','latex')
 string_title={'\textbf{Two Compartment Model for  $a=$}'+string(a)};
 title(string_title,'Interpreter','latex');
 legend('\textbf{GL Scheme y1}','\textbf{GL Scheme y2}','\textbf{ode15s y1}','\textbf{ode15s y2}','Interpreter','latex')
 grid on;
 hold on;

                                             %----------------------------------------------%
                                                             % Defining the Model %
                                             %----------------------------------------------%
 


function f=ff(t,y,params)
%----------------------------------------------%
             % Model Parameters %
%----------------------------------------------%
 k1 = params(1);
 kd = params(2);
 
 %----------------------------------------------%
             % Model Equations %
%----------------------------------------------%
 f(1) = -k1*y(1);
 f(2) = k1*y(1)-kd*y(2);
end
 
function fdot=ffdot(t,y,params)
%----------------------------------------------%
             % Model Parameters %
%----------------------------------------------%
k1 = params(1);
k2 = params(2);

%----------------------------------------------%
             % Model Jackobian %
%----------------------------------------------%
fdot(1,1)=-k1;
fdot(1,2)=0;
fdot(2,1)=k1;
fdot(2,2)=-k2;
end  
                                             %----------------------------------------------%
                                                             % ODE15s Function %
                                             %----------------------------------------------%
                                             
function dydt = rate(~,y,k1,k2)
dydt(1) = -k1*y(1);
dydt(2)=k1*y(1)-k2*y(2);
dydt=[dydt(1),dydt(2)]';
end