clear;

DOA_C=[-30,0,20];       %center DOAs of  scattering clusters
N_c=length(DOA_C);      %number of scattering clusters
N_s=20;                 %number of sub-paths per scattering
L=N_c*N_s;              %number of DoAs in totoal

DOA_all=[];
for cc=1:N_c
    DOA_all=[DOA_all, DOA_C(cc)+rand(1,N_s)*10];
end
      
N=150;                 % number of antennas
T=60;                  % Snapshot
N_point=N;             % number of grid points
SNR=10;
Ps=sqrt((   10.^(SNR/10)   )/2);
X=Ps*( randn(T,N)+1i*randn(T,N)  );
A = exp(-1i*pi*(0:N-1)'*sind(DOA_all));
H=A*(randn(L,1)+1i*randn(L,1));
noise=sqrt(1/2)*(randn(T,1)+1i*randn(T,1));
Y=X*H + noise;
etc=T;                  % or etc=length(DOA) if the number of DoAs is known


H_burst=Burst_CE(Y,X,N,N_point,etc);
norm(H-H_burst,'fro')^2/norm(H)^2