%% defining A,B,E,H,epsilon. Matrices should satisfy disturbance decoupling condition
eigenvalue = 1/10*[-5 0 0 0 0 0; 0 -50.1 0 0 0 0; 0 0 -10.3 0 0 0; 0 0 0 -30.4 0 0; 0 0 0 0 -40.5 0; 0 0 0 0 0 -60.6];
eigenvector = [0.1 -20 5 2 20 4; 0 9 -25 20 5 -10;2 -17 15 4 1 0; 5 16 -19 0 2 1;-6 11 -2 5 6 0;-8 9 -1 5 7 -16];
%eigenvector = [0.1 -20 5 2 20 -0.1667; 0 9 -25 20 5 -0.0417;2 -17 15 4 1 -0.0083; 5 16 -19 0 2 -0.0167;-6 11 -2 5 6 -0.05;-8 9 -1 5 7 -0.0583];
A = eigenvector*eigenvalue*(eigenvector^(-1));
B = [0 0 0;0 0 0;1 0 0;0 1 0;0 0 1;0 0 0];
%B = [0 0; 0 0; 1 0; 0 1; 0 0; 0 0];
E = [1;0;0;0;0;0];
%  H = [0 0 0 1 0 0; 0 0 0 0 1 0];
%H = [-1 3 5 1/2 9 8; -3 6 8 1 5 0] ;
H = [0 -1 0 1 0 0;0 0 0 0 -1 0];
epsilon = 0.1;

%% Find feedback matrix for disturbance decoupling
Vstar=null(H);
for i = 1:(size(H,2)-rank(H))  %rank-nullity theorem, giving dimension of kernel of H
  Z = null([Vstar';B']);    %nu-star algorithm
  Vstar = null([H;Z'*A]);
end
X = [];
for i = 1:size(Vstar,2)
    X =  [X linsolve([Vstar B],(A*Vstar(:,i)))];    %A*Vi = [V* B][X;Y]
end
 p = size(X,1);
 U_part = X((p-size(B,2)+1):p,:);
 U_ker_f = null([Vstar B]);
 U_ker = U_ker_f((p-size(B,2)+1):p,:);
 F_part = -(eye(size(U_ker,1))-(U_ker)*gen_inv(U_ker))*(-U_part)*gen_inv(-Vstar);
 transposeF = -linsolve(transpose(Vstar),transpose(X(p-1:p,:)));
 %Vstar_gi = gen_inv(Vstar);
 F = transpose(transposeF);
 
 %% proposing random initial conditions
%     x0 = [0;0;0;0;0;0];
   x0 =1/10*[0;0.1;4;0.15/sqrt(2);0.25/sqrt(2);0.2];
%    x0 = [0.595197704931052;-1.15722386019289;1.39115876943583;-1.21896742934121;0.0756117477677569;-0.203947480368569];
%   x0 = [0.182802211161469;0.101763171101706;-0.0951342391596149;0.0701977367720205;0.0446852844068269;0.0554638532199223];
%    x0 = 1/100*[10;56;-24;-30;80;1];

 %% finding second term of the expression
 T = 0.95;
 step = 0.001;
 n = (T/step);
 z_0 = H*(expm(A*T) + (expm(A*T)-eye(size(A)))*inv(A)*B*F_part)*x0;
 g_0 =  H*((expm(A*T)-eye(size(A)))*inv(A)*B);
 AA = g_0;
 BB = x0;
 CC  = g_0*(-(eye(size(-U_ker,1))-(U_ker)*gen_inv(U_ker)));
 DD = Vstar*gen_inv(Vstar)*x0;
 G = kron(BB',AA)+kron(DD',CC);
 ZZZ = linsolve(G,vec(-z_0));
 ZZ = reshape(ZZZ,[size(B,2),size(x0)]);
 F_ker = ZZ + (U_ker*gen_inv(U_ker)-eye(size(U_ker,1)))*ZZ*Vstar*gen_inv(Vstar);
 F_des = F_part+F_ker;
 norm(F_des,2)
 u_i = [0];
 for i = 1:(2^(size(g_0,2))-1)
    u_i = [u_i (u_i(end)+1)];
 end
 u_f = int2bit(u_i,size(g_0,2));
 u_max = F_ker*x0;

 u_f(u_f == 0) = -max(abs(u_max));
 u_f(u_f == 1) = max(abs(u_max));
 z_des = z_0 + g_0*u_max;
 z_f = z_0 + g_0*u_f;
 z_f_t = z_f';
 k_f = convhull(z_f_t);
 traj_f =[];
 norm_traj_f=[];
 derivative_norm_traj_f=[];
 traj_total_f =[];
 norm_traj_total_f =[];
 state_f =[];
 for t = 0:step:2*T
    traj_f =    [traj_f H*(expm(A.*t)*x0 +(expm(A.*t)-eye(size(A)))*inv(A)*B*F_part*x0)+H*((expm(A*t)-eye(size(A)))*inv(A)*B)*u_max];
    norm_traj_f = [norm_traj_f norm(traj_f(:,end))];
    derivative_norm_traj_f = [derivative_norm_traj_f (H*(expm(A.*t)*x0 +(expm(A.*t)-eye(size(A)))*inv(A)*B*F_part*x0)+H*((expm(A*t)-eye(size(A)))*inv(A)*B)*u_max)'*H*expm(A*t)*(A+B*F_des)*x0];
    traj_total_f = [traj_total_f H*(expm(A.*t)*x0 +(expm(A.*t)-eye(size(A)))*inv(A)*B*F_part*x0)+H*((expm(A*t)-eye(size(A)))*inv(A)*B)*u_max+H*(expm(A*t)-eye(size(A)))*inv(A)*E*(-1)];
    norm_traj_total_f =[norm_traj_total_f norm(traj_total_f(:,end))]; 
    state_f = [state_f (expm(A.*t)*x0 +(expm(A.*t)-eye(size(A)))*inv(A)*B*F_part*x0)+((expm(A*t)-eye(size(A)))*inv(A)*B)*u_max+(expm(A*t)-eye(size(A)))*inv(A)*E*(1)];
 end
 t = 0:step:2*T;
%  plot(t,norm_traj_f)
 grid on
  hold on
%  plot(t,derivative_norm_traj_f)
%  hold on
%   plot(t,norm_traj_total_f)
 hold on
 eps = zeros(floor(T/step)+1,1);
 eps(eps==0) = epsilon;
%   plot(t,eps)
 hold on
%  plot(z_0(1),z_0(2),'o-','MarkerFaceColor','white','MarkerEdgeColor','red')
%  plot(z_des(1),z_des(2),'x-','MarkerFaceColor','white','MarkerEdgeColor','green')
 hold on
 grid on
% req_vec = [0;0]-z_0;
% u_0 = linsolve(g_0,req_vec);
% f_16_ker = u_0(1)/x0(6);
% f_26_ker = u_0(2)/x0(6);

%% finding zonotope for disturbance
step_sampling_time = 0.1;
A_discrete = expm(A*step_sampling_time);
E_discrete = (expm(A*step_sampling_time)-eye(6))*inv(A)*E;
g_d = [];
for i = 0:1:((T/step_sampling_time)-1)
    g_d = [g_d (A_discrete^i)*E_discrete];
end
n = size(g_d,2);
u = [0];
for i = 1:(2^n-1)
    u = [u (u(end)+1)];
end
u_d = int2bit(u,n);
u_d(u_d == 0) = -1;
% u = [1 0 -1 0 1/sqrt(2) -1/sqrt(2) 1/sqrt(2) -1/sqrt(2);0 1 0 -1 1/sqrt(2) -1/sqrt(2) -1/sqrt(2) 1/sqrt(2)];
% u_d = repmat(u,(n/2),1);
z_d = H*g_d*u_d;
x_d = g_d*u_d;
%  plot(z_d(1,:),z_d(2,:))
 hold on
z_d_t = z_d';
k_d = convhull(z_d_t);
% plot(z_d_t(:,1),z_d_t(:,2))
hold on
plot(z_d_t(k_d,1),z_d_t(k_d,2))
grid on
hold on
%theta = 0:pi/50:2*pi;
%plot(epsilon*cos(theta),epsilon*sin(theta),'m')
hold on
grid on

