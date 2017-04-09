%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%    Example of ADI Method for 2D heat equation                        %
%                                                                      %
%          u_t = u_{xx} + u_{yy} + f(x,t)                              %
%                                                                      %
%    Test problme:                                                     %
%      Exact solution: u(t,x,y) = exp(-t) sin(pi*x) sin(pi*y)          %
%      Source term:  f(t,x,y) = exp(-t) sin(pi*x) sin(pi*y) (2pi^2-1)  %
%                                                                      %
%    Files needed for the test:                                        %
%                                                                      %
%     adi.m:      This file, the main calling code.                    %
%     f.m:        The file defines the f(t,x,y)                        %
%     uexact.m:    The exact solution.                                 %
%                                                                      %
%     Results:         n              e            ratio               %
%                     10           0.0041                              %
%     t_final=0.5     20           0.0010           4.1                %
%                     40           2.5192e-04       3.97               %
%                     80           6.3069e-05       3.9944             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
tic;

%====================================
%INITIAL VALUES
%====================================
a=0; 
b=2;  
c=0; 
d=2;
tfinal = 0.5;
N = 80;  
M = N;
Mp=M+1;
Np=N+1;
h = (b-a)/N;
h1 = h*h;


x=a:h:b;
y=c:h:d;
% define the mesh in time
dt=0.01;%h;
k_t=fix(tfinal/dt)

r=dt/(h1);
%-- Initial condition:
    V1 = zeros(Np,Np);
    V2 = zeros(Np,Np);
    F = zeros(Np,Np);
   t = 0;
   for i=1:Np,
      for j=1:Np,
         V1(i,j) = uexact(t,x(i),y(j));
      end
   end
   W=V1;

    Ax = zeros(Np,Np);
    Ay = Ax;
           
    for i=1:Np
        if i==1
            Ax(i,i+1) = -r;
        elseif i==Np             
            Ax(i,i-1) =  -r;
        else
            Ax(i,i+1) = -r/2;          
            Ax(i,i-1) = -r/2;
        end
        Ax(i,i) = 1+r;
    end
    
    for j=1:Np         
        if j==1            
            Ay(j,j+1) = -r;
        elseif j==Np   
            Ay(j,j-1) =  -r;
        else
            Ay(j,j+1) = -r/2;      
            Ay(j,j-1) = -r/2;
        end
        Ay(j,j) = 1+r;
    end


%---------- Big loop for time t --------------------------------------


for k=1:k_t   
    %--- sweep in x-direction --------------------------------------    
        for j = 1:Np,                             % Look for fixed y(j)      
            b=zeros(Np,1);     
            for i=1:Np
                if j==1
                    b(i) = V1(i,j) + (r)*( -V1(i,j) + V1(i,j+1));
                elseif j==Np
                    b(i) = V1(i,j) + (r)*( V1(i,j-1) -V1(i,j));
                else
                    b(i) = V1(i,j) + (r/2)*(V1(i,j-1) -2*V1(i,j) + V1(i,j+1));
                end
%                   b(i-1) = (V1(i,j-1) -2*V1(i,j) + V1(i,j+1))/h1 + 2*C_m/dt*V1(i,j);
%                 if i == 2  
%                     b(i-1) = b(i-1); %+ u_exact(t,x(i-1),y(j))/h1; 
%                 else i==N
%                     b(i-1) = b(i-1); %+ u_exact(t,x(i+1),y(j))/h1;
%                 end
            end
            Vt = Ax\b;                          % Solve the diagonal matrix. 
            V2(1:Np,j) = Vt;
        end
    % Strang Splitting RK2 time-integration% Finish x-sweep.
    %==== solve ODEs ====
%F(:,:)=0;%(2*pi^2-2)*V2;%-2*W;

    %-------- RK2 for ODE ---------------------------
    for i=1:Np
        V_tmp(i,1:Np)=V2(i,1:Np) + 1/2*dt*F(i,1:Np);
        %---- compute m H J d f X Ca -----------------
        W_tmp(i,1:Np)  = W(i,1:Np) + 1/2*dt*(-2*V2(i,1:Np));    
    end % finish RK2 stage1
    
%F=0;%(2*pi^2-2)*V_tmp;%-2*W_tmp;
   
    for i=1:Np
        V2(i,1:Np)=V2(i,1:Np) + dt*F(i,1:Np);
        %---- compute m H J d f X Ca -----------------
        W(i,1:Np) = W(i,1:Np)  + dt*(-2*V_tmp(i,1:Np)); 
    end % finish RK2

 
    %-------------- loop in y -direction -------------------------------- 
    
    for i = 1:Np
        b=zeros(Np,1);
        for j=1:Np
            if i==1
                b(j) = V2(i,j) + (r/2)*( -2*V2(i,j) + 2*V2(i+1,j));
            elseif i==Np
                b(j) = V2(i,j) + (r/2)*( 2*V2(i-1,j) -2*V2(i,j) );
            else
                b(j) = V2(i,j) + (r/2)*(V2(i-1,j) -2*V2(i,j) + V2(i+1,j));
            end
%               b(j-1) = (V2(i-1,j) -2*V2(i,j) + V2(i+1,j))/h1 + 2*C_m/dt*V2(i,j);
%             if j == 2
%                 b(j-1) = b(j-1);% + u_exact(t1,x(i),y(j-1))/h1;     
%             else j==Np        
%                 b(j-1) = b(j-1);% + u_exact(t1,x(i),y(j+1))/h1;
%             end
        end
        ut = Ay\b;
         V1(i,1:Np) = ut';
    end

    %--- finish ADI method at this time level, go to the next time level.
                  
%     if mod(k,100)==0
%     %      set(wave_handle,'ZData',V1);
%     set(wave_handle2,'YData',y1);
%     title(sprintf('Time %f',t))               
%     drawnow;
%     end


end  % Finished with the loop in time
  for i=1:Np,
    for j=1:Np,
       ue(i,j) = uexact(tfinal,x(i),y(j));
    end
  end

  e = max(max(abs(V1-ue)))        % The infinity error.
toc

    
    
    
   