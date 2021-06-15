function [C_h] = H_memory_capacity_for_linear_memory_tasks(lag,a,b,c,tstep,gamma,xp,L,x,u_z,lamda)
%C 此处显示有关此函数的摘要
%   此处显示详细说明
Num_of_VitualNodes=fix(lag/tstep);
A=DxF(tstep,gamma,a,b,c,xp,Num_of_VitualNodes);
eta=log(1+tstep*b);
dfx=a*(1+(1-c)*xp^c)/((1+xp^c)^2);
dffx=a*(c*(1-c)*xp^(c-1)*(1+xp^c)-2*(1+(1-c)*xp^c))/((1+xp^c)^3);
q_r=f_q_r(u_z,x,dfx,dffx,eta,Num_of_VitualNodes);
cov_yx=zeros(1,Num_of_VitualNodes);
for i=1:Num_of_VitualNodes
   cov_yx(i)=1-exp(-eta);
   for j=1:length(L)
       tmp=0;
       tmpA=A^(j-1);
       for r=1:Num_of_VitualNodes
           tmp=tmp+L(j)*tmpA(i,r)*f_p_r(u_z,r,x,dfx,dffx,eta);
       end
   end
   cov_yx(i)=cov_yx(i)*tmp;
end
sigama_e=zeros(Num_of_VitualNodes,Num_of_VitualNodes);
for i=1:Num_of_VitualNodes
    for j=1:Num_of_VitualNodes
        sigama_e(i,j)=(1-exp(-eta))^2*(f_q_r_i_j(u_z,i,j,x,dfx,dffx,eta)-q_r(i)*q_r(j));
    end
end
V_sigama_e=Matrix2Vector(sigama_e);
tmpGa =pinv(eye(Num_of_VitualNodes^2,Num_of_VitualNodes^2)-kron(A,A));
V_sigama_0=tmpGa*V_sigama_e';
sigama_0=VectorToMatrix(V_sigama_0);
Vary=u_z(2)*norm(L,2);

tmp1=pinv(sigama_0+lamda*eye(Num_of_VitualNodes,Num_of_VitualNodes));
tmp2=sigama_0+2*lamda*eye(Num_of_VitualNodes,Num_of_VitualNodes);
C_h=(cov_yx*tmp1*tmp2*tmp1*cov_yx')/Vary;
end

function y=f_q_r(z,x,dfx,dffx,eta,N)
y=zeros(1,N);
tmpq_r_1=z(1)*dfx;
tmpq_r_2=(z(2)^2/2)*dffx;
for i=1:N
    tmp=zeros(1,i);
    for j=1:i
        tmp(j)=exp(-(i-j)*eta);
    end
    y(i)=tmpq_r_1*(tmp*x(1:i)')+tmpq_r_2*(tmp*(x(1:i).^2)');
end
end

function y=f_p_r(z,i,x,dfx,dffx,eta)
tmp_i=zeros(1,i);
for k=1:i
     tmp_i(k)=exp(-(i-k)*eta);
end
y=z(2)*dfx*(tmp_i*x(1:i)')+z(3)*dffx*(tmp_i*(x(1:i).^2)');
end

function y=f_q_r_i_j(z,i,j,x,dfx,dffx,eta)
% tmpq_r_1=z(1)*dfx;
% tmpq_r_2=(z(2)^2/2)*dffx;
tmp_i=zeros(1,i);
tmp_j=zeros(1,j);
for k=1:i
     tmp_i(k)=exp(-(i-k)*eta);
end
for k=1:j
     tmp_j(k)=exp(-(j-k)*eta);
end
ai_1=dfx*(tmp_i*x(1:i)');
ai_2=dffx*(tmp_i*(x(1:i).^2)');
aj_1=dfx*(tmp_j*x(1:j)');
aj_2=dffx*(tmp_j*(x(1:j).^2)');
y=z(2)*ai_1*aj_1+z(3)*(ai_1*aj_2+ai_2*aj_1)+z(4)*(ai_2*aj_2);
end



function [A]=DxF(d,gamma,a,b,c,x,N)
xp=(ones(N,1).*x).^c;
g=a.*(1+(1-c).*xp)./((1+xp).^2);
% g=cos(a.*xp).*a.*c.*x.^(c-1);
% g=cosh(a.*xp).^-2.*a.*c.*x.^(c-1);
% g=a.*c.*x.^(c-1);
eta=log(1+d*b);
tmp1=exp(-eta);
tmp2=b^(-1)*(1-tmp1).*g;
% tmp2=b^(-1)*(1-exp(-eta)).*g;
A=zeros(N,N);
for j=1:N
    if j~=N
        for i=j:N
            A(i,j)=tmp2(i)*tmp1^(i-j);
        end
    elseif j==N
        for i=1:N
            A(i,j)=tmp1^i;
        end
        A(N,N)=A(N,N)+tmp2(i);
    end
end
end

