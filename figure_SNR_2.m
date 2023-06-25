%这个代码中的模型接收端预编码矩阵不是时变的
%横坐标为信噪比
%prop
%5个真实相位
%全部真实相位
% clc;
%% 预设值
close all;
clear;
% par.runId = 1;
% rng(par.runId);%这两行代码用于固定程序中的随机量
Nit1=100;%实验重复次数
% rhos_db=-50:10:-10;
rhos_db=-10;
rhos=power(10,rhos_db/10);
L_SNR=length(rhos_db);
v=5000;%振荡器线宽,单位Hz
Nt=1;%接收天线数
Nr=16;%发射天线数
Nrf=2;%射频链数
L=4;%信道抽头
ofdm_num=L;
Nt_beam=Nt;%假设的值
Nr_beam=Nr;%假设的值
Nx=Nr_beam/Nrf;
N=Nt_beam*Nr_beam/Nrf;%载波数量
Nc=4;
Nsc=4;%子路径数
path_num=Nc*Nsc;%总路径数量
velocity_light=3e8;%光速
fc=1e9;
W=20*1e6;%系统带宽 w=0.02*fc
lambda=velocity_light/fc;%波长
d=lambda/2;%天线间距离
Ts=0.05*(1e-6);%符号时间
delta_x_square=10;
I_Nrf=eye(Nrf);
theta_r=zeros(path_num,1);
theta_t=zeros(path_num,1);%Nsc个为一簇
sub_theta_r=zeros(Nc,1);
sub_theta_t=zeros(Nc,1);%Nc个簇
%% 预编码矩阵
%%%时域
F_RF_time=zeros(Nt,N);%模拟预编码矩阵
sub_F_RF=zeros(Nt,Nt_beam);
for i=1:Nt_beam
    for j=1:Nt
        sub_F_RF(j,i)=(1/sqrt(Nt))*exp(-sqrt(-1)*(2*pi/lambda)*(j-1)*d*((2/Nt_beam)*(i-1)-1));
    end
end
for i=1:N
    xx=ceil(i/Nx);
    if (xx<=Nt)
        xx=ceil(i/Nx);
    else
        xx=xx-(floor((xx-0.1)/Nt_beam))*Nt_beam;
    end
    %     F_RF_time_index(:,i)=xx;
    F_RF_time(:,i)=sub_F_RF(:,xx);
end
W_RF_time=zeros(Nr,N*Nrf);%模拟预编码矩阵
sub_W_RF=zeros(Nr,Nr_beam);
for i=1:Nr_beam
    for j=1:Nr
        sub_W_RF(j,i)=(1/sqrt(Nr))*exp(-sqrt(-1)*(2*pi/lambda)*(j-1)*d*((2/Nr_beam)*(i-1)-1));
    end
end
for i=1:N*Nrf/Nr_beam
    W_RF_time(:,(i-1)*Nr_beam+1:i*Nr_beam)=sub_W_RF;
end
%%%频域
F_RF=zeros(Nt,N);%模拟预编码矩阵
for i=1:N
    for j=1:N
        F_RF(:,i)=F_RF(:,i)+F_RF_time(:,j)*exp(-sqrt(-1)*2*pi*(1/N)*(i-1)*(j-1));
    end
end
% test_F_RF=F_RF'*F_RF;
W_RF=zeros(Nr,N*Nrf);%模拟预编码矩阵
for i=1:N
    for j=1:N
        W_RF(:,(i-1)*Nrf+1:i*Nrf)=W_RF(:,(i-1)*Nrf+1:i*Nrf)+W_RF_time(:,(j-1)*Nrf+1:j*Nrf)*exp(-sqrt(-1)*2*pi*(1/N)*(i-1)*(j-1));
    end
end
%%%
%% 生成AR和AT
M=32;%把0-2pi角度划分成M个
AR=zeros(Nr,M);
for i=1:M
    for j=1:Nr
        AR(j,i)=(1/sqrt(Nr))*exp(-sqrt(-1)*(2*pi/lambda)*(j-1)*d*((2/M)*(i-1)-1));
    end
end
AT=zeros(Nt,M);
for i=1:M
    for j=1:Nt
        AT(j,i)=(1/sqrt(Nt))*exp(-sqrt(-1)*(2*pi/lambda)*(j-1)*d*((2/M)*(i-1)-1));
    end
end
%% 生成矩阵G
I_L=eye(L);
submatrix_ARAT=kron(AR,conj(AT));%维度：NrNt*MM
G=kron(I_L,submatrix_ARAT);
%% 生成P为N*M的矩阵(相位噪声功率谱密度未知)
M1=8;
P_LL=zeros(N,M1);
for n=1:N
    for m=1:M1
        if (n-1)>=(m-1)*(N-1)/(M1-1)&&(n-1)<m*(N-1)/(M1-1)
            P_LL(n,m)=m-((n-1)*(M1-1)/(N-1));
        elseif (n-1)>=(m-2)*(N-1)/(M1-1)&&(n-1)<(m-1)*(N-1)/(M1-1)
            P_LL(n,m)=((n-1)*(M1-1)/(N-1))-(m-2);
        else
            P_LL(n,m)=0;
        end
    end
end
P=P_LL;
%% 生成F_g为N*N的矩阵
Fg=zeros(N,N);
for i=1:N
    for j=1:N
        Fg(i,j)=exp(-sqrt(-1)*2*pi*(1/N)*(i-1)*(j-1));
    end
end
%% 生成MM为N*L的矩阵(单根天线的信道傅里叶变换矩阵)
MM=ones(N,L);
for i=1:N
    for j=1:L
        MM(i,j)=exp(-sqrt(-1)*2*pi*(1/N)*(i-1)*(j-1));
    end
end
%% 设计MIMO的信道傅里叶变换矩阵
I_NrNt=eye(Nr*Nt);
F_h=kron(MM,I_NrNt);
%% 实验重复开始
tic
for kkk=1:Nit1
    kkk
    %% 生成随机角度
    sigma_AS=1/180*pi;
    a_A=rand(Nc,1)*2*pi;
    b_A=sigma_AS/sqrt(2);
    a_D=rand(Nc,1)*2*pi;
    b_D=sigma_AS/sqrt(2);
    for kk=1:Nsc
        AoA(:,kk)=laplace(a_A,b_A,Nc,1);
        AoD(:,kk)=laplace(a_D,b_D,Nc,1);
    end
    for i=1:Nc
        theta_r(Nsc*(i-1)+1:Nsc*i)=AoA(i,:).';
    end
    for i=1:Nc
        theta_t(Nsc*(i-1)+1:Nsc*i)=AoD(i,:).';
    end
    matrix_a_theta_r=zeros(Nr,path_num);
    for i=1:path_num
        for j=1:Nr
            matrix_a_theta_r(j,i)=(1/sqrt(Nr))*exp(-sqrt(-1)*(2*pi/lambda)*(j-1)*d*cos(theta_r(i)));
        end
    end
    matrix_a_theta_t=zeros(Nt,path_num);
    for i=1:path_num
        for j=1:Nt
            matrix_a_theta_t(j,i)=(1/sqrt(Nt))*exp(-sqrt(-1)*(2*pi/lambda)*(j-1)*d*cos(theta_t(i)));
        end
    end
    alpha_path_gain=(randn(path_num,1)+randn(path_num,1)*sqrt(-1))/sqrt(2);
    %% 时域滤波器建模
    roll_down_factor=1;
    beta=roll_down_factor;
    filter_t=zeros(L,path_num);
    t_times=zeros(L,path_num);
    path_delay=unifrnd(0,L*Ts,path_num,1);% 均匀分布：unifrnd (a, b, m, n)； 产生m*n阶[a, b]均匀分布,每个簇内时延一样
    for i=1:L
        for j=1:path_num
            t_times(i,j)=(i-1)*Ts-path_delay(j,1);%行表示不同抽样时刻，列表示每一抽样时刻不同簇内的不同时延
        end
    end
    for i=1:L
        for j=1:path_num
            if  t_times(i,j)==Ts/(2*beta)||t_times(i,j)==-Ts/(2*beta)
                filter_t(i,j)=(pi/4)*((sin((1/(2*beta))*pi))/((1/(2*beta))*pi));
            else
                filter_t(i,j)=((sin((t_times(i,j)/(Ts))*pi))/((t_times(i,j)/(Ts))*pi))*((cos((pi*beta*t_times(i,j))/Ts))/(1-(2*beta*t_times(i,j)/Ts)^2));
            end
        end
    end
    %% 时域信道
    HHt=zeros(Nr,Nt,L);
    for n=1:L
        for i=1:path_num
            HHt(:,:,n)=HHt(:,:,n)+sqrt((Nt*Nr)/(path_num))*alpha_path_gain(i)*filter_t(n,i)*matrix_a_theta_r(:,i)*matrix_a_theta_t(:,i)';
        end
    end
    Ht=cell(L,L);
    for i=1:L
        for j=1:L
            if(i==j)
                Ht{i,i}=HHt(:,:,i);
            else
                Ht{i,j}=zeros(Nr,Nt);
            end
        end
    end
    Ht=cell2mat(Ht);
    %%%%%%%%%%%%%%%%%%%%%%%%%%将时域信道向量化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hht=zeros(Nr*Nt,1,L);%分页
    for l=1:L
        for i=1:Nr
            for j=1:Nt
                hht((i-1)*Nt+j,1,l)=HHt(i,j,l);
            end
        end
    end
    ht=zeros(Nr*Nt*L,1);
    for i=1:L
        for j=1:Nr*Nt
            ht((i-1)*Nr*Nt+j,1)=hht(j,1,i);
        end
    end
    %% 频域滤波器建模
    filter_k=zeros(N,path_num);%生成滤波器频域值，傅里叶变换,行表示不同子载波，列表示每一子载波不同簇内的不同时延
    for k=1:N
        for i=1:L
            filter_k(k,:)=filter_k(k,:)+filter_t(i,:)*exp(-sqrt(-1)*2*pi*(k-1)*(i-1)/N);
        end
    end
    %% 频域信道
    HH=zeros(Nr,Nt,N);%i*j*k   %频域
    for k=1:N
        for i=1:path_num
            HH(:,:,k)=HH(:,:,k)+sqrt((Nt*Nr)/(path_num))*alpha_path_gain(i)*filter_k(k,i)*matrix_a_theta_r(:,i)*matrix_a_theta_t(:,i)';
        end
    end
    H=cell(N,N);
    for i=1:N
        for j=1:N
            if(i==j)
                H{i,i}=HH(:,:,i);
            else
                H{i,j}=zeros(Nr,Nt);
            end
        end
    end
    H=cell2mat(H);
    power_H_sum=0;
    for i=1:N*Nr
        for j=1:N*Nt
            power_H_sum=power_H_sum+H(i,j)'*H(i,j);
        end
    end
    %%%%%%%%%%%%%%%%%%%%将频域信道向量化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hhk=zeros(Nr*Nt,1,N);%分页
    for l=1:N
        for i=1:Nr
            for j=1:Nt
                hhk((i-1)*Nt+j,1,l)=HH(i,j,l);
            end
        end
    end
    hk=zeros(Nr*Nt*N,1);
    for i=1:N
        for j=1:Nr*Nt
            hk((i-1)*Nr*Nt+j,1)=hhk(j,1,i);
        end
    end
    %信道的时频转换核对
    %     hk_test=F_h*ht;
    %     hk_test_error=hk_test-hk;
    %% 生成不同的OFDM符号
    ofdm_num=N;
    for ofdmi=1:ofdm_num
        ofdmi;
        %% MIMO发射信号建模
        nsymbol=N;%一共有N个符号
        MM_QAM=16;%阶数，表示16QAM
        graycode=[0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10];%格雷映射，十进制表示
        msg=randi([0,15],1,nsymbol);%随机产生0-15的符号，
        msg1=graycode(msg+1);%格雷映射
        msgmod=qammod(msg1,MM_QAM);%调用qammod函数，得到调制后的符号
        s1=msgmod';%得到发送信号x(列向量)
        s(:,1)=s1;
        s_array(:,ofdmi)=s;
        %% MIMO接收端相位噪声建模,GeR建模
        N0_Gr=v/pi;%单边带功率谱密度
        ud_Gr=Ts*normrnd(0,sqrt(N0_Gr*W),1,N);%N0*W为功率，功率为方差
        F_Gr=ones(1,N);%此矩阵储存相位噪声
        for i=1:N
            u_sum_Gr=0;
            for j=1:i-1
                u_sum_Gr=u_sum_Gr+ud_Gr(j);
            end
            F_Gr(i)=2*pi*(u_sum_Gr);
        end
        a_Gr=zeros(N,1);%此矩阵储存频域相位噪声
        for k=1:N
            sum_A_Gr=0;
            for n=1:N
                sum_A_Gr=sum_A_Gr+(1/N)*exp(sqrt(-1)*F_Gr(n))*exp(-sqrt(-1)*2*pi*(k-1)*(n-1)/N);
            end
            a_Gr(k)=sum_A_Gr;
        end
        %         a_Gr=zeros(N,1);%此矩阵储存接收端频域相位噪声
        %         a_Gr(1)=1;%去除相位噪声
        a_Gr_array(:,ofdmi)=a_Gr;
        R=kron(a_Gr,I_Nrf);
        Re=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
        for i=1:N
            Re((i-1)*Nrf+1:i*Nrf,:)=R((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
        end
        G_eR=zeros(Nrf*N,Nr*N);
        G_eR(:,1:Nr)=Re;
        Ree=Re;
        for i=2:N
            Ree=circshift(Ree,[Nrf,0]);%由a循环移位得到A
            G_eR(:,(i-1)*Nr+1:i*Nr)=Ree;
        end
        %% MIMO发射端相位噪声建模,GeT建模
        N0_Gt=v/pi;%单边带功率谱密度
        ud_Gt=Ts*normrnd(0,sqrt(N0_Gt*W),1,N);%N0*W为功率，功率为方差
        F_Gt=ones(1,N);%此矩阵储存相位噪声
        for i=1:N
            u_sum_Gt=0;
            for j=1:i-1
                u_sum_Gt=u_sum_Gt+ud_Gt(j);
            end
            F_Gt(i)=2*pi*(u_sum_Gt);
        end
        a_Gt=zeros(N,1);%此矩阵储存频域相位噪声
        for k=1:N
            sum_A_Gt=0;
            for n=1:N
                sum_A_Gt=sum_A_Gt+(1/N)*exp(sqrt(-1)*F_Gt(n))*exp(-sqrt(-1)*2*pi*(k-1)*(n-1)/N);
            end
            a_Gt(k)=sum_A_Gt;
        end
        %     a_Gt=zeros(N,1);%此矩阵储存发射端频域相位噪声
        %     a_Gt(1)=1;%去除相位噪声
        a_Gt_array(:,ofdmi)=a_Gt;
        t=a_Gt;
        te=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
        for i=1:N
            te((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t(i);
        end
        G_eT=zeros(Nt*N,N);
        G_eT(:,1)=te;
        tee=te;
        for i=2:N;
            tee= circshift(tee,[Nt,0]);%由循环移位得到
            G_eT(:,i)=tee;
        end
        %% 接收信号
        sub_y=G_eR*H*G_eT*s;
        y((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,1)=sub_y;%生成总的y
    end
    y_store=y;
    %% 不同信噪比开始
    for iii=1:L_SNR%每次迭代使用不同信噪比
        iii
        %% 接收信号
        %生成高斯白噪声向量
        power_w=1;
        w=sqrt(power_w/2)*(randn(Nrf*N*N,1)+sqrt(-1)*randn(Nrf*N*N,1));%生成高斯白噪声
        y=y_store*sqrt(rhos(iii))+w;
        %% prop
        adjust_num=3;%调整次数
        adjust_pn_num=5;%调整相位噪声系数数目
        r=y;%残差矩阵初值
        K=60;
        %% 相位噪声初始化
        a_Gr_array_et0=zeros(N,ofdm_num);
        a_Gr_array_et0(1,:)=ones(1,ofdm_num);
        a_Gt_array_et0=zeros(N,ofdm_num);
        a_Gt_array_et0(1,:)=ones(1,ofdm_num);
        a_Gt_array_et=a_Gt_array_et0;
        a_Gr_array_et=a_Gr_array_et0;
        array_ht_error_record_prop4_7_4=zeros(1,1);
        array_value_error_y_record_prop4_7_4=zeros(1,1);
        %% 调整发射端相位噪声
        for j4_transmit=1:adjust_num
            j4_transmit;
            %% 确定调整第几个OFDM符号
            for j5_transmit=1:N
                %% 计算未调整相位噪声之前的初始y误差范数
                a_Gt_array_et=a_Gt_array_et;
                for ofdmi=1:ofdm_num%获得总的Q
                    %% 获得发射端相位噪声等效矩阵用于估计信道
                    t_et=a_Gt_array_et(:,ofdmi);
                    te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                    for i=1:N
                        te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                    end
                    G_eT_et=zeros(Nt*N,N);
                    G_eT_et(:,1)=te_et;
                    tee_et=te_et;
                    for i=2:N;
                        tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                        G_eT_et(:,i)=tee_et;
                    end
                    %% 获得接收端相位噪声等效矩阵用于估计信道
                    R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                    Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                    for i=1:N
                        Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                    end
                    G_eR_et=zeros(Nrf*N,Nr*N);
                    G_eR_et(:,1:Nr)=Re_et;
                    Ree_et=Re_et;
                    for i=2:N
                        Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                        G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                    end
                    %% 获得sub_Q和Q
                    matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                    T_h_k=zeros(Nr,Nr*Nt);
                    T_h=zeros(Nr*N,Nr*Nt*N);
                    for i=1:N
                        T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                    end
                    for i=1:N
                        T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                    end
                    sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                    Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                end
                vector_j_data=zeros(K,1);%储存最大值数据
                r=y;
                Q_j=zeros(Nrf*N*ofdm_num,1);
                for nn=1:K
                    nn;
                    z=Q'*r;%(M*M*L,1)
                    for i=1:M*M*L
                        if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                            vector_j_data(nn,1)=z(i,1);
                            vector_j(nn,1)=i;%储存最大值的位置
                        end
                    end
                    Q_j(:,nn)=Q(:,vector_j(nn));
                    h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                    r=y-Q_j*h_et_sparse_value;
                    array2(nn,1)=r'*r;
                end
                h_vitural_et=zeros(M*M*L,1);
                for i=1:nn
                    h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                end
                ht_et=G*h_vitural_et;
                error_h=ht_et-ht;
                error_h_norm=(error_h'*error_h)/(ht'*ht);
                %记录y残差范数
                y_et=Q*h_vitural_et;
                value_error_y_initial=(y-y_et)'*(y-y_et);
                a_Gt_array_et_initial=a_Gt_array_et;
                if j4_transmit==1&&j5_transmit==1
                    array_ht_error_record_prop4_7_4(1)=error_h_norm;
                    y_et=Q*h_vitural_et;
                    array_value_error_y_record_prop4_7_4(1)=(y-y_et)'*(y-y_et);
                end
                %% 确定10个值调整方向
                for j1=1:adjust_pn_num%共调整五个值（adjust_pn_num=5）
                    for j3=1:2%分实部虚部
                        %% 调整步长
                        if j1==1||j1==4
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.008;
                            end
                        end
                        if j1==2||j1==3
                            if j3==1
                                step_length=0.005;
                            end
                            if j3==2
                                step_length=0.004;
                            end
                        end
                        if j1==5
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.03;
                            end
                        end
                        %% 判断调整第几个值
                        if j1==1%调整a[2]
                            j2=2;
                        end
                        if j1==2%调整a[3]
                            j2=3;
                        end
                        if j1==3%调整a[29]
                            j2=N-1;
                        end
                        if j1==4%调整a[30]
                            j2=N;
                        end
                        if j1==5%调整a[1]
                            j2=1;
                        end
                        %% 确定方向
                        if j3==1%实部
                            %% 计算判断调整方向（向上）用到的对比值
                            a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)+step_length;
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                                Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge1=(y-y_et)'*(y-y_et);
                            a_Gt_array_et(:,j5_transmit)=a_Gt_array_et_initial(:,j5_transmit);
                            %% 计算判断调整方向（向下）用到的对比值
                            a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)-step_length;
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                                Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge2=(y-y_et)'*(y-y_et);
                            a_Gt_array_et(:,j5_transmit)=a_Gt_array_et_initial(:,j5_transmit);
                            %% 储存方向调整系数
                            array_direct_judge_value_error_y(1)=value_error_y_initial;
                            array_direct_judge_value_error_y(2)=value_error_y_direct_judge1;
                            array_direct_judge_value_error_y(3)=value_error_y_direct_judge2;
                            if min(array_direct_judge_value_error_y)==value_error_y_initial
                                array_direct_judge((j1-1)*2+j3,j5_transmit)=1;%确定不调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge1
                                array_direct_judge((j1-1)*2+j3,j5_transmit)=2;%确定向上调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge2
                                array_direct_judge((j1-1)*2+j3,j5_transmit)=3;%确定向下调整
                            end
                            if j4_transmit==1&&j1==5&&j3==1
                                if array_direct_judge((j1-1)*2+j3,j5_transmit)==2
                                    array_direct_judge((j1-1)*2+j3,j5_transmit)=1;
                                end
                            end
                        end
                        if j3==2%虚部
                            %% 计算判断调整方向（向上）用到的对比值
                            a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)+step_length*sqrt(-1);
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge1=(y-y_et)'*(y-y_et);
                            a_Gt_array_et(:,j5_transmit)=a_Gt_array_et_initial(:,j5_transmit);
                            %% 计算判断调整方向（向下）用到的对比值
                            a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)-step_length*sqrt(-1);
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                                Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge2=(y-y_et)'*(y-y_et);
                            a_Gt_array_et(:,j5_transmit)=a_Gt_array_et_initial(:,j5_transmit);
                            %% 储存方向调整系数
                            array_direct_judge_value_error_y(1)=value_error_y_initial;
                            array_direct_judge_value_error_y(2)=value_error_y_direct_judge1;
                            array_direct_judge_value_error_y(3)=value_error_y_direct_judge2;
                            if min(array_direct_judge_value_error_y)==value_error_y_initial
                                array_direct_judge((j1-1)*2+j3,j5_transmit)=1;%确定不调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge1
                                array_direct_judge((j1-1)*2+j3,j5_transmit)=2;%确定向上调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge2
                                array_direct_judge((j1-1)*2+j3,j5_transmit)=3;%确定向下调整
                            end
                        end
                    end
                end
                %% 对10个值统一进行调整
                phase_noise_result_stored_1=a_Gt_array_et;
                for j1=1:adjust_pn_num%共调整五个值
                    for j3=1:2%分实部虚部
                        %% 调整步长
                        if j1==1||j1==4
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.008;
                            end
                        end
                        if j1==2||j1==3
                            if j3==1
                                step_length=0.005;
                            end
                            if j3==2
                                step_length=0.004;
                            end
                        end
                        if j1==5
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.03;
                            end
                        end
                        %% 判断调整第几个值
                        if j1==1%调整a[2]
                            j2=2;
                        end
                        if j1==2%调整a[3]
                            j2=3;
                        end
                        if j1==3%调整a[29]
                            j2=N-1;
                        end
                        if j1==4%调整a[30]
                            j2=N;
                        end
                        if j1==5%调整a[1]
                            j2=1;
                        end
                        %% 统一调整
                        if j3==1%调整实部
                            if  array_direct_judge((j1-1)*2+j3,j5_transmit)==1;%确定不调整
                                a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit);
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_transmit)==2;%确定向上调整
                                a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)+step_length;
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_transmit)==3;%确定向下调整
                                a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)-step_length;
                            end
                        end
                        if j3==2%调整虚部
                            if  array_direct_judge((j1-1)*2+j3,j5_transmit)==1;%确定不调整
                                a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit);
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_transmit)==2;%确定向上调整
                                a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)+step_length*sqrt(-1);
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_transmit)==3;%确定向下调整
                                a_Gt_array_et(j2,j5_transmit)=a_Gt_array_et(j2,j5_transmit)-step_length*sqrt(-1);
                            end
                        end
                    end
                end
                %% 根据时域相位噪声调整频域相位噪声幅值
                a_Gt_et_time_domain=N*inv(Fg)*a_Gt_array_et(:,j5_transmit);%转到时域
                alpha=1/a_Gt_et_time_domain(1);%计算幅值系数
                a_Gt_array_et(:,j5_transmit)=a_Gt_array_et(:,j5_transmit)*alpha;%进一步调整相位噪声
                phase_noise_result_stored_2=a_Gt_array_et;
                %% 记录估计过程中的信道误差
                for ofdmi=1:ofdm_num%获得总的Q
                    %% 获得发射端相位噪声等效矩阵用于估计信道
                    t_et=a_Gt_array_et(:,ofdmi);
                    te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                    for i=1:N
                        te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                    end
                    G_eT_et=zeros(Nt*N,N);
                    G_eT_et(:,1)=te_et;
                    tee_et=te_et;
                    for i=2:N;
                        tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                        G_eT_et(:,i)=tee_et;
                    end
                    %% 获得接收端相位噪声等效矩阵用于估计信道
                    R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                    Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                    for i=1:N
                        Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                    end
                    G_eR_et=zeros(Nrf*N,Nr*N);
                    G_eR_et(:,1:Nr)=Re_et;
                    Ree_et=Re_et;
                    for i=2:N
                        Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                        G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                    end
                    %% 获得sub_Q和Q
                    matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                    T_h_k=zeros(Nr,Nr*Nt);
                    T_h=zeros(Nr*N,Nr*Nt*N);
                    for i=1:N
                        T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                    end
                    for i=1:N
                        T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                    end
                    sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                    Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                    Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                end
                vector_j_data=zeros(K,1);%储存最大值数据
                r=y;
                Q_j=zeros(Nrf*N*ofdm_num,1);
                for nn=1:K
                    nn;
                    z=Q'*r;%(M*M*L,1)
                    for i=1:M*M*L
                        if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                            vector_j_data(nn,1)=z(i,1);
                            vector_j(nn,1)=i;%储存最大值的位置
                        end
                    end
                    Q_j(:,nn)=Q(:,vector_j(nn));
                    h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                    r=y-Q_j*h_et_sparse_value;
                    array2(nn,1)=r'*r;
                end
                h_vitural_et=zeros(M*M*L,1);
                for i=1:nn
                    h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                end
                ht_et=G*h_vitural_et;
                error_h=ht_et-ht;
                error_h_norm=(error_h'*error_h)/(ht'*ht);
                array_ht_error_record_prop4_7_4((j4_transmit-1)*L+j5_transmit+1)=error_h_norm;
                %% 记录估计过程中的y误差范数
                y_et=Q*h_vitural_et;
                value_error_y_final=(y-y_et)'*(y-y_et);
                array_value_error_y_record_prop4_7_4((j4_transmit-1)*L+j5_transmit+1)=value_error_y_final;
                %% 最终判断是否调整
                if  array_value_error_y_record_prop4_7_4((j4_transmit-1)*N+j5_transmit+1)>=array_value_error_y_record_prop4_7_4((j4_transmit-1)*N+j5_transmit)
                    a_Gt_array_et=phase_noise_result_stored_1;
                    array_value_error_y_record_prop4_7_4((j4_transmit-1)*N+j5_transmit+1)=array_value_error_y_record_prop4_7_4((j4_transmit-1)*N+j5_transmit);
                    array_ht_error_record_prop4_7_4((j4_transmit-1)*N+j5_transmit+1)=array_ht_error_record_prop4_7_4((j4_transmit-1)*N+j5_transmit);
                end
            end
            if array_value_error_y_record_prop4_7_4((j4_transmit-1)*N+j5_transmit+1)==array_value_error_y_record_prop4_7_4((j4_transmit-1)*N+j5_transmit+1-N)
                break;
            end
        end
        %%%发射端相位噪声调整完毕
        %% 调整接收端相位噪声
        for j4_receive=1:adjust_num
            j4_receive;
            %% 确定调整第几个OFDM符号
            for j5_receive=1:N
                %% 计算未调整相位噪声之前的初始y误差范数
                a_Gr_array_et=a_Gr_array_et;
                for ofdmi=1:ofdm_num%获得总的Q
                    %% 获得发射端相位噪声等效矩阵用于估计信道
                    t_et=a_Gt_array_et(:,ofdmi);
                    te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                    for i=1:N
                        te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                    end
                    G_eT_et=zeros(Nt*N,N);
                    G_eT_et(:,1)=te_et;
                    tee_et=te_et;
                    for i=2:N;
                        tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                        G_eT_et(:,i)=tee_et;
                    end
                    %% 获得接收端相位噪声等效矩阵用于估计信道
                    R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                    Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                    for i=1:N
                        Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                    end
                    G_eR_et=zeros(Nrf*N,Nr*N);
                    G_eR_et(:,1:Nr)=Re_et;
                    Ree_et=Re_et;
                    for i=2:N
                        Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                        G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                    end
                    %% 获得sub_Q和Q
                    matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                    T_h_k=zeros(Nr,Nr*Nt);
                    T_h=zeros(Nr*N,Nr*Nt*N);
                    for i=1:N
                        T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                    end
                    for i=1:N
                        T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                    end
                    sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                    Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                end
                vector_j_data=zeros(K,1);%储存最大值数据
                r=y;
                Q_j=zeros(Nrf*N*ofdm_num,1);
                for nn=1:K
                    nn;
                    z=Q'*r;%(M*M*L,1)
                    for i=1:M*M*L
                        if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                            vector_j_data(nn,1)=z(i,1);
                            vector_j(nn,1)=i;%储存最大值的位置
                        end
                    end
                    Q_j(:,nn)=Q(:,vector_j(nn));
                    h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                    r=y-Q_j*h_et_sparse_value;
                    array2(nn,1)=r'*r;
                end
                h_vitural_et=zeros(M*M*L,1);
                for i=1:nn
                    h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                end
                ht_et=G*h_vitural_et;
                error_h=ht_et-ht;
                error_h_norm=(error_h'*error_h)/(ht'*ht);
                %记录y残差范数
                y_et=Q*h_vitural_et;
                value_error_y_initial=(y-y_et)'*(y-y_et);
                a_Gr_array_et_initial=a_Gr_array_et;
                %% 确定10个值调整方向
                for j1=1:adjust_pn_num%共调整五个值（adjust_pn_num=5）
                    for j3=1:2%分实部虚部
                        %% 调整步长
                        if j1==1||j1==4
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.008;
                            end
                        end
                        if j1==2||j1==3
                            if j3==1
                                step_length=0.005;
                            end
                            if j3==2
                                step_length=0.004;
                            end
                        end
                        if j1==5
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.03;
                            end
                        end
                        %% 判断调整第几个值
                        if j1==1%调整a[2]
                            j2=2;
                        end
                        if j1==2%调整a[3]
                            j2=3;
                        end
                        if j1==3%调整a[29]
                            j2=N-1;
                        end
                        if j1==4%调整a[30]
                            j2=N;
                        end
                        if j1==5%调整a[1]
                            j2=1;
                        end
                        %% 确定方向
                        if j3==1%实部
                            %% 计算判断调整方向（向上）用到的对比值
                            a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)+step_length;
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                                Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge1=(y-y_et)'*(y-y_et);
                            a_Gr_array_et(:,j5_receive)=a_Gr_array_et_initial(:,j5_receive);
                            %% 计算判断调整方向（向下）用到的对比值
                            a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)-step_length;
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                                Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge2=(y-y_et)'*(y-y_et);
                            a_Gr_array_et(:,j5_receive)=a_Gr_array_et_initial(:,j5_receive);
                            %% 储存方向调整系数
                            array_direct_judge_value_error_y(1)=value_error_y_initial;
                            array_direct_judge_value_error_y(2)=value_error_y_direct_judge1;
                            array_direct_judge_value_error_y(3)=value_error_y_direct_judge2;
                            if min(array_direct_judge_value_error_y)==value_error_y_initial
                                array_direct_judge((j1-1)*2+j3,j5_receive)=1;%确定不调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge1
                                array_direct_judge((j1-1)*2+j3,j5_receive)=2;%确定向上调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge2
                                array_direct_judge((j1-1)*2+j3,j5_receive)=3;%确定向下调整
                            end
                            if j4_receive==1&&j1==5&&j3==1
                                if array_direct_judge((j1-1)*2+j3,j5_receive)==2
                                    array_direct_judge((j1-1)*2+j3,j5_receive)=1;
                                end
                            end
                        end
                        if j3==2%虚部
                            %% 计算判断调整方向（向上）用到的对比值
                            a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)+step_length*sqrt(-1);
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge1=(y-y_et)'*(y-y_et);
                            a_Gr_array_et(:,j5_receive)=a_Gr_array_et_initial(:,j5_receive);
                            %% 计算判断调整方向（向下）用到的对比值
                            a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)-step_length*sqrt(-1);
                            for ofdmi=1:ofdm_num%获得总的Q
                                %% 获得发射端相位噪声等效矩阵用于估计信道
                                t_et=a_Gt_array_et(:,ofdmi);
                                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                                for i=1:N
                                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                                end
                                G_eT_et=zeros(Nt*N,N);
                                G_eT_et(:,1)=te_et;
                                tee_et=te_et;
                                for i=2:N;
                                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                                    G_eT_et(:,i)=tee_et;
                                end
                                %% 获得接收端相位噪声等效矩阵用于估计信道
                                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                                for i=1:N
                                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                                end
                                G_eR_et=zeros(Nrf*N,Nr*N);
                                G_eR_et(:,1:Nr)=Re_et;
                                Ree_et=Re_et;
                                for i=2:N
                                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                                end
                                %% 获得sub_Q和Q
                                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                                T_h_k=zeros(Nr,Nr*Nt);
                                T_h=zeros(Nr*N,Nr*Nt*N);
                                for i=1:N
                                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                                end
                                for i=1:N
                                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                                end
                                sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                                Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                            end
                            vector_j_data=zeros(K,1);%储存最大值数据
                            r=y;
                            Q_j=zeros(Nrf*N*ofdm_num,1);
                            for nn=1:K
                                nn;
                                z=Q'*r;%(M*M*L,1)
                                for i=1:M*M*L
                                    if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                                        vector_j_data(nn,1)=z(i,1);
                                        vector_j(nn,1)=i;%储存最大值的位置
                                    end
                                end
                                Q_j(:,nn)=Q(:,vector_j(nn));
                                h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                                r=y-Q_j*h_et_sparse_value;
                                array2(nn,1)=r'*r;
                            end
                            h_vitural_et=zeros(M*M*L,1);
                            for i=1:nn
                                h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                            end
                            ht_et=G*h_vitural_et;
                            error_h=ht_et-ht;
                            error_h_norm=(error_h'*error_h)/(ht'*ht);
                            %%记录y残差范数
                            y_et=Q*h_vitural_et;
                            value_error_y_direct_judge2=(y-y_et)'*(y-y_et);
                            a_Gr_array_et(:,j5_receive)=a_Gr_array_et_initial(:,j5_receive);
                            %% 储存方向调整系数
                            array_direct_judge_value_error_y(1)=value_error_y_initial;
                            array_direct_judge_value_error_y(2)=value_error_y_direct_judge1;
                            array_direct_judge_value_error_y(3)=value_error_y_direct_judge2;
                            if min(array_direct_judge_value_error_y)==value_error_y_initial
                                array_direct_judge((j1-1)*2+j3,j5_receive)=1;%确定不调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge1
                                array_direct_judge((j1-1)*2+j3,j5_receive)=2;%确定向上调整
                            end
                            if min(array_direct_judge_value_error_y)==value_error_y_direct_judge2
                                array_direct_judge((j1-1)*2+j3,j5_receive)=3;%确定向下调整
                            end
                        end
                    end
                end
                %% 对10个值统一进行调整
                phase_noise_result_stored_1=a_Gr_array_et;
                for j1=1:adjust_pn_num%共调整五个值
                    for j3=1:2%分实部虚部
                        %% 调整步长
                        if j1==1||j1==4
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.008;
                            end
                        end
                        if j1==2||j1==3
                            if j3==1
                                step_length=0.005;
                            end
                            if j3==2
                                step_length=0.004;
                            end
                        end
                        if j1==5
                            if j3==1
                                step_length=0.01;
                            end
                            if j3==2
                                step_length=0.03;
                            end
                        end
                        %% 判断调整第几个值
                        if j1==1%调整a[2]
                            j2=2;
                        end
                        if j1==2%调整a[3]
                            j2=3;
                        end
                        if j1==3%调整a[29]
                            j2=N-1;
                        end
                        if j1==4%调整a[30]
                            j2=N;
                        end
                        if j1==5%调整a[1]
                            j2=1;
                        end
                        %% 统一调整
                        if j3==1%调整实部
                            if  array_direct_judge((j1-1)*2+j3,j5_receive)==1;%确定不调整
                                a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive);
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_receive)==2;%确定向上调整
                                a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)+step_length;
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_receive)==3;%确定向下调整
                                a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)-step_length;
                            end
                        end
                        if j3==2%调整虚部
                            if  array_direct_judge((j1-1)*2+j3,j5_receive)==1;%确定不调整
                                a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive);
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_receive)==2;%确定向上调整
                                a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)+step_length*sqrt(-1);
                            end
                            if  array_direct_judge((j1-1)*2+j3,j5_receive)==3;%确定向下调整
                                a_Gr_array_et(j2,j5_receive)=a_Gr_array_et(j2,j5_receive)-step_length*sqrt(-1);
                            end
                        end
                    end
                end
                %% 根据时域相位噪声调整频域相位噪声幅值
                a_Gr_et_time_domain=N*inv(Fg)*a_Gr_array_et(:,j5_receive);%转到时域
                alpha=1/a_Gr_et_time_domain(1);%计算幅值系数
                a_Gr_array_et(:,j5_receive)=a_Gr_array_et(:,j5_receive)*alpha;%进一步调整相位噪声
                phase_noise_result_stored_2=a_Gr_array_et;
                %% 记录估计过程中的信道误差
                for ofdmi=1:ofdm_num%获得总的Q
                    %% 获得发射端相位噪声等效矩阵用于估计信道
                    t_et=a_Gt_array_et(:,ofdmi);
                    te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                    for i=1:N
                        te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                    end
                    G_eT_et=zeros(Nt*N,N);
                    G_eT_et(:,1)=te_et;
                    tee_et=te_et;
                    for i=2:N;
                        tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                        G_eT_et(:,i)=tee_et;
                    end
                    %% 获得接收端相位噪声等效矩阵用于估计信道
                    R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                    Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                    for i=1:N
                        Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                    end
                    G_eR_et=zeros(Nrf*N,Nr*N);
                    G_eR_et(:,1:Nr)=Re_et;
                    Ree_et=Re_et;
                    for i=2:N
                        Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                        G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                    end
                    %% 获得sub_Q和Q
                    matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                    T_h_k=zeros(Nr,Nr*Nt);
                    T_h=zeros(Nr*N,Nr*Nt*N);
                    for i=1:N
                        T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                    end
                    for i=1:N
                        T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                    end
                    sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
                    Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
                    Q1((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=G_eR_et*T_h*F_h;
                end
                vector_j_data=zeros(K,1);%储存最大值数据
                r=y;
                Q_j=zeros(Nrf*N*ofdm_num,1);
                for nn=1:K
                    nn;
                    z=Q'*r;%(M*M*L,1)
                    for i=1:M*M*L
                        if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                            vector_j_data(nn,1)=z(i,1);
                            vector_j(nn,1)=i;%储存最大值的位置
                        end
                    end
                    Q_j(:,nn)=Q(:,vector_j(nn));
                    h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
                    r=y-Q_j*h_et_sparse_value;
                    array2(nn,1)=r'*r;
                end
                h_vitural_et=zeros(M*M*L,1);
                for i=1:nn
                    h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
                end
                ht_et=G*h_vitural_et;
                error_h=ht_et-ht;
                error_h_norm=(error_h'*error_h)/(ht'*ht);
                array_ht_error_record_prop4_7_4(j4_transmit*L+(j4_receive-1)*L+j5_receive+1)=error_h_norm;
                %% 记录估计过程中的y误差范数
                y_et=Q*h_vitural_et;
                value_error_y_final=(y-y_et)'*(y-y_et);
                array_value_error_y_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive+1)=value_error_y_final;
                %% 最终判断是否调整
                if  array_value_error_y_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive+1)>=array_value_error_y_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive)
                    a_Gr_array_et=phase_noise_result_stored_1;
                    array_value_error_y_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive+1)=array_value_error_y_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive);
                    array_ht_error_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive+1)=array_ht_error_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive);
                end
            end
            if array_value_error_y_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive+1)==array_value_error_y_record_prop4_7_4(j4_transmit*N+(j4_receive-1)*N+j5_receive+1-N)
                break;
            end
        end
        %%%接收端相位噪声调整完毕
        %% 统计最终结果
        a_Gt_array_et=a_Gt_array_et;
        a_Gr_array_et=a_Gr_array_et;
        for ofdmi=1:ofdm_num%获得总的Q
            %% 获得发射端相位噪声等效矩阵用于估计信道
            t_et=a_Gt_array_et(:,ofdmi);
            te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
            for i=1:N
                te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
            end
            G_eT_et=zeros(Nt*N,N);
            G_eT_et(:,1)=te_et;
            tee_et=te_et;
            for i=2:N;
                tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                G_eT_et(:,i)=tee_et;
            end
            %% 获得接收端相位噪声等效矩阵用于估计信道
            R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
            Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
            for i=1:N
                Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
            end
            G_eR_et=zeros(Nrf*N,Nr*N);
            G_eR_et(:,1:Nr)=Re_et;
            Ree_et=Re_et;
            for i=2:N
                Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
            end
            %% 获得sub_Q和Q
            matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
            T_h_k=zeros(Nr,Nr*Nt);
            T_h=zeros(Nr*N,Nr*Nt*N);
            for i=1:N
                T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
            end
            for i=1:N
                T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
            end
            sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
            Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
        end
        vector_j_data=zeros(K,1);%储存最大值数据
        r=y;
        Q_j=zeros(Nrf*N*ofdm_num,1);
        for nn=1:K
            nn;
            z=Q'*r;%(M*M*L,1)
            for i=1:M*M*L
                if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                    vector_j_data(nn,1)=z(i,1);
                    vector_j(nn,1)=i;%储存最大值的位置
                end
            end
            Q_j(:,nn)=Q(:,vector_j(nn));
            h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
            r=y-Q_j*h_et_sparse_value;
            array2(nn,1)=r'*r;
        end
        h_vitural_et=zeros(M*M*L,1);
        for i=1:nn
            h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
        end
        ht_et=G*h_vitural_et;
        error_h=ht_et-ht;
        error_h_norm=(error_h'*error_h)/(ht'*ht);
        array_error_prop4_7_4(kkk,iii)=error_h_norm;
        %% 参考论文中的对照方法
        a_Gr_array_et0=zeros(N,ofdm_num);
        a_Gr_array_et0(1,:)=ones(1,ofdm_num);
        a_Gt_array_et0=zeros(N,ofdm_num);
        a_Gt_array_et0(1,:)=ones(1,ofdm_num);
        %     a_Gt_array_et=a_Gt_array;
        %     a_Gr_array_et=a_Gr_array;
        adjust_num=10;%算法迭代次数
        a_Gt_array_et=a_Gt_array_et0;
        a_Gr_array_et=a_Gr_array_et0;
        %% 迭代开始
        for jjj=1:adjust_num%联合估计算法开始
            jjj
            %% 估计信道
            for ofdmi=1:ofdm_num%获得总的Q
                %% 获得发射端相位噪声等效矩阵用于估计信道
                t_et=a_Gt_array_et(:,ofdmi);
                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                for i=1:N
                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                end
                G_eT_et=zeros(Nt*N,N);
                G_eT_et(:,1)=te_et;
                tee_et=te_et;
                for i=2:N;
                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                    G_eT_et(:,i)=tee_et;
                end
                %% 获得接收端相位噪声等效矩阵用于估计信道
                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                for i=1:N
                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                end
                G_eR_et=zeros(Nrf*N,Nr*N);
                G_eR_et(:,1:Nr)=Re_et;
                Ree_et=Re_et;
                for i=2:N
                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                end
                %% 获得sub_Q和Q
                matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
                T_h_k=zeros(Nr,Nr*Nt);
                T_h=zeros(Nr*N,Nr*Nt*N);
                for i=1:N
                    T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
                end
                for i=1:N
                    T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
                end
                sub_Q=G_eR_et*T_h*F_h*sqrt(rhos(iii));
                Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
            end
            %用LS算法获得信道估计值
            ht_et=inv(Q'*Q)*Q'*y;%估计信道
            error_h=ht_et-ht;
            error_h_norm=(error_h'*error_h)/(ht'*ht);
            %% 统计一次实验中每次迭代后的信道误差
            array_ht_error_record_con1_1(jjj)=error_h_norm;
            %% 计算y残差并设置跳出条件
            y_et=Q*ht_et;
            array_value_error_y_con1_1(jjj)=(y-y_et)'*(y-y_et);
            if jjj>1
                if abs(array_value_error_y_con1_1(jjj)-array_value_error_y_con1_1(jjj-1))<0.01*abs(array_value_error_y_con1_1(2)-array_value_error_y_con1_1(1))||(array_value_error_y_con1_1(jjj)>array_value_error_y_con1_1(jjj-1))
                    break;
                end
            end
            %% 估计发射端相位噪声
            for ofdmi=1:N
                %%获得接收端相位噪声等效矩阵用于估计信道
                R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
                Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
                for i=1:N
                    Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
                end
                G_eR_et=zeros(Nrf*N,Nr*N);
                G_eR_et(:,1:Nr)=Re_et;
                Ree_et=Re_et;
                for i=2:N
                    Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                    G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
                end
                %估计发射端相位噪声
                hk_et=F_h*ht_et;
                hhk_et=zeros(Nr*Nt,1,N);%分页
                for i=1:N
                    for j=1:Nr*Nt
                        hhk_et(j,1,i)=hk_et((i-1)*Nr*Nt+j,1);
                    end
                end
                HH_et=zeros(Nr,Nt,N);%i*j*k   %频域
                for l=1:N
                    for i=1:Nr
                        for j=1:Nt
                            HH_et(i,j,l)=hhk_et((i-1)*Nt+j,1,l);
                        end
                    end
                end
                H_et=cell(N,N);
                for i=1:N
                    for j=1:N
                        if(i==j)
                            H_et{i,i}=HH_et(:,:,i);
                        else
                            H_et{i,j}=zeros(Nr,Nt);
                        end
                    end
                end
                H_et=cell2mat(H_et);
                HeT_et=zeros(Nr*N,Nrf_t*N);
                for i=1:N
                    HeT_et((i-1)*Nr+1:i*Nr,(i-1)*Nrf_t+1:i*Nrf_t)=H_et((i-1)*Nr+1:i*Nr,(i-1)*Nt+1:i*Nt)*F_RF(:,1);
                end
                T_TC=zeros(Nrf_t*N,N);
                x1=s_array(:,ofdmi);
                T_TC(:,1)=x1;
                for i=2:N
                    x1= circshift(x1,[Nrf_t,0]);%由a循环移位得到A
                    T_TC(:,i)=x1;
                end
                Q_TC=G_eR_et*HeT_et*T_TC*Fg*P*sqrt(rhos(iii));
                %             rank3=rank(T_TC);
                %             rank2=rank(Q_TC);
                a_Gt_p_et=N*inv(Q_TC'*Q_TC)*Q_TC'*y((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,1);
                a_Gt_array_et(:,ofdmi)=(1/N)*Fg*P*a_Gt_p_et;
            end
            %% 估计接收端相位噪声
            for ofdmi=1:N
                HeR_et=zeros(Nrf*N,Nt*N);
                for i=1:N
                    HeR_et((i-1)*Nrf+1:i*Nrf,(i-1)*Nt+1:i*Nt)=(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))'*H_et((i-1)*Nr+1:i*Nr,(i-1)*Nt+1:i*Nt);
                end
                %%获得发射端相位噪声等效矩阵用于估计信道
                t_et=a_Gt_array_et(:,ofdmi);
                te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
                for i=1:N
                    te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
                end
                G_eT_et=zeros(Nt*N,N);
                G_eT_et(:,1)=te_et;
                tee_et=te_et;
                for i=2:N;
                    tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                    G_eT_et(:,i)=tee_et;
                end
                T_RC=zeros(Nrf*N,N);
                x2=HeR_et*G_eT_et*s_array(:,ofdmi);
                T_RC(:,1)=x2;
                for i=2:N
                    x2= circshift(x2,[Nrf,0]);%由a循环移位得到A
                    T_RC(:,i)=x2;
                end
                Q_RC=T_RC*Fg*P*sqrt(rhos(iii));
                rank1=rank(Q_RC);
                a_Gr_p_et=N*inv(Q_RC'*Q_RC)*Q_RC'*y((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,1);
                a_Gr_array_et(:,ofdmi)=(1/N)*Fg*P*a_Gr_p_et;
            end
        end%联合估计结束
        %% 统计最终结果
        array_ht_et_error_final_con1_1(kkk,iii)=array_ht_error_record_con1_1(jjj);%ht_et是最终获得的信道估计值
        %% 方法五：全部初始化相位（用OMP计算信道）
        K=rank(Q);
        a_Gt_array_et=a_Gt_array_et0;
        a_Gr_array_et=a_Gr_array_et0;
        for ofdmi=1:ofdm_num%获得总的Q
            %% 获得发射端相位噪声等效矩阵用于估计信道
            t_et=a_Gt_array_et(:,ofdmi);
            te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
            for i=1:N
                te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
            end
            G_eT_et=zeros(Nt*N,N);
            G_eT_et(:,1)=te_et;
            tee_et=te_et;
            for i=2:N;
                tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                G_eT_et(:,i)=tee_et;
            end
            %% 获得接收端相位噪声等效矩阵用于估计信道
            R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
            Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
            for i=1:N
                Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
            end
            G_eR_et=zeros(Nrf*N,Nr*N);
            G_eR_et(:,1:Nr)=Re_et;
            Ree_et=Re_et;
            for i=2:N
                Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
            end
            %% 获得sub_Q和Q
            matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
            T_h_k=zeros(Nr,Nr*Nt);
            T_h=zeros(Nr*N,Nr*Nt*N);
            for i=1:N
                T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
            end
            for i=1:N
                T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
            end
            sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
            Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
        end
        vector_j_data=zeros(K,1);%储存最大值数据
        r=y;
        Q_j=zeros(Nrf*N*ofdm_num,1);
        for nn=1:K
            nn;
            z=Q'*r;%(M*M*L,1)
            for i=1:M*M*L
                if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                    vector_j_data(nn,1)=z(i,1);
                    vector_j(nn,1)=i;%储存最大值的位置
                end
            end
            Q_j(:,nn)=Q(:,vector_j(nn));
            h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
            r=y-Q_j*h_et_sparse_value;
            array2(nn,1)=r'*r;
        end
        h_vitural_et=zeros(M*M*L,1);
        for i=1:nn
            h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
        end
        ht_et=G*h_vitural_et;
        error_h=ht_et-ht;
        error_h_norm=(error_h'*error_h)/(ht'*ht);
        array_error_method5(kkk,iii)=error_h_norm;
        %方法六估计结束
        %% 方法六：全部真实相位（用OMP计算信道）
        a_Gt_array_et=a_Gt_array;
        a_Gr_array_et=a_Gr_array;
        for ofdmi=1:ofdm_num%获得总的Q
            %% 获得发射端相位噪声等效矩阵用于估计信道
            t_et=a_Gt_array_et(:,ofdmi);
            te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
            for i=1:N
                te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
            end
            G_eT_et=zeros(Nt*N,N);
            G_eT_et(:,1)=te_et;
            tee_et=te_et;
            for i=2:N;
                tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                G_eT_et(:,i)=tee_et;
            end
            %% 获得接收端相位噪声等效矩阵用于估计信道
            R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
            Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
            for i=1:N
                Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
            end
            G_eR_et=zeros(Nrf*N,Nr*N);
            G_eR_et(:,1:Nr)=Re_et;
            Ree_et=Re_et;
            for i=2:N
                Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
            end
            %% 获得sub_Q和Q
            matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
            T_h_k=zeros(Nr,Nr*Nt);
            T_h=zeros(Nr*N,Nr*Nt*N);
            for i=1:N
                T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
            end
            for i=1:N
                T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
            end
            sub_Q=G_eR_et*T_h*F_h*G*sqrt(rhos(iii));
            Q((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
        end
        vector_j_data=zeros(K,1);%储存最大值数据
        r=y;
        Q_j=zeros(Nrf*N*ofdm_num,1);
        for nn=1:K
            nn;
            z=Q'*r;%(M*M*L,1)
            for i=1:M*M*L
                if(abs(z(i,1))>abs(vector_j_data(nn,1)))
                    vector_j_data(nn,1)=z(i,1);
                    vector_j(nn,1)=i;%储存最大值的位置
                end
            end
            Q_j(:,nn)=Q(:,vector_j(nn));
            h_et_sparse_value=inv(Q_j'*Q_j)*Q_j'*y;
            r=y-Q_j*h_et_sparse_value;
            array2(nn,1)=r'*r;
        end
        h_vitural_et=zeros(M*M*L,1);
        for i=1:nn
            h_vitural_et(vector_j(i,1),1)=h_et_sparse_value(i,1);
        end
        ht_et=G*h_vitural_et;
        error_h=ht_et-ht;
        error_h_norm=(error_h'*error_h)/(ht'*ht);
        array_error_method6(kkk,iii)=error_h_norm;
        %方法六估计结束
        %% LS上界
        a_Gt_array_et=a_Gt_array_et0;
        a_Gr_array_et=a_Gr_array_et0;
        for ofdmi=1:ofdm_num%获得总的Q
            %% 获得发射端相位噪声等效矩阵用于估计信道
            t_et=a_Gt_array_et(:,ofdmi);
            te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
            for i=1:N
                te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
            end
            G_eT_et=zeros(Nt*N,N);
            G_eT_et(:,1)=te_et;
            tee_et=te_et;
            for i=2:N;
                tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                G_eT_et(:,i)=tee_et;
            end
            %% 获得接收端相位噪声等效矩阵用于估计信道
            R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
            Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
            for i=1:N
                Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
            end
            G_eR_et=zeros(Nrf*N,Nr*N);
            G_eR_et(:,1:Nr)=Re_et;
            Ree_et=Re_et;
            for i=2:N
                Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
            end
            %% 获得sub_Q和Q
            matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
            T_h_k=zeros(Nr,Nr*Nt);
            T_h=zeros(Nr*N,Nr*Nt*N);
            for i=1:N
                T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
            end
            for i=1:N
                T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
            end
            sub_Q=G_eR_et*T_h*F_h*sqrt(rhos(iii));
            Q_LS((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
        end
        ht_et=inv(Q_LS'*Q_LS)*Q_LS'*y;%估计信道
        error_h=ht_et-ht;
        error_h_norm=(error_h'*error_h)/(ht'*ht);
        LS_error_h_upper(kkk,iii)=error_h_norm;
        %% LS下界
        a_Gt_array_et=a_Gt_array;
        a_Gr_array_et=a_Gr_array;
        for ofdmi=1:ofdm_num%获得总的Q
            %% 获得发射端相位噪声等效矩阵用于估计信道
            t_et=a_Gt_array_et(:,ofdmi);
            te_et=zeros(Nt*N,1);%排成一列,从上往下第一块是te[0],第二块是te[1],......
            for i=1:N
                te_et((i-1)*Nt+1:i*Nt,:)=F_RF(:,1)*t_et(i);
            end
            G_eT_et=zeros(Nt*N,N);
            G_eT_et(:,1)=te_et;
            tee_et=te_et;
            for i=2:N;
                tee_et= circshift(tee_et,[Nt,0]);%由循环移位得到
                G_eT_et(:,i)=tee_et;
            end
            %% 获得接收端相位噪声等效矩阵用于估计信道
            R_et=kron(a_Gr_array_et(:,ofdmi),I_Nrf);
            Re_et=zeros(Nrf*N,Nr);%排成一列,从上往下第一块是Re[0],第二块是Re[1],......
            for i=1:N
                Re_et((i-1)*Nrf+1:i*Nrf,:)=R_et((i-1)*Nrf+1:i*Nrf,:)*(W_RF(:,(ofdmi-1)*Nrf+1:ofdmi*Nrf))';
            end
            G_eR_et=zeros(Nrf*N,Nr*N);
            G_eR_et(:,1:Nr)=Re_et;
            Ree_et=Re_et;
            for i=2:N
                Ree_et= circshift(Ree_et,[Nrf,0]);%由a循环移位得到A
                G_eR_et(:,(i-1)*Nr+1:i*Nr)=Ree_et;
            end
            %% 获得sub_Q和Q
            matrix_Gtx=G_eT_et*s_array(:,ofdmi);%x是导频
            T_h_k=zeros(Nr,Nr*Nt);
            T_h=zeros(Nr*N,Nr*Nt*N);
            for i=1:N
                T_h_k(:,:,i)=kron(eye(Nr),matrix_Gtx((i-1)*Nt+1:i*Nt).');
            end
            for i=1:N
                T_h((i-1)*Nr+1:i*Nr,(i-1)*Nr*Nt+1:i*Nr*Nt)=T_h_k(:,:,i);
            end
            sub_Q=G_eR_et*T_h*F_h*sqrt(rhos(iii));
            Q_LS((ofdmi-1)*Nrf*N+1:ofdmi*Nrf*N,:)=sub_Q;
        end
        ht_et=inv(Q_LS'*Q_LS)*Q_LS'*y;%估计信道
        error_h=ht_et-ht;
        error_h_norm=(error_h'*error_h)/(ht'*ht);
        LS_error_h_lower(kkk,iii)=error_h_norm;
    end
end
t=toc
for i=1:L_SNR
    mean_array_error_prop4_7_4(1,i)=mean(array_error_prop4_7_4(:,i));
    mean_array_error_method5(1,i)=mean(array_error_method5(:,i));%OMP上界
    mean_array_error_method6(1,i)=mean(array_error_method6(:,i));%OMP下界
    mean_LS_error_h_upper(1,i)=mean(LS_error_h_upper(:,i));
    mean_LS_error_h_lower(1,i)=mean(LS_error_h_lower(:,i));
end
mean_array_error_prop4_7_4_db=10*log10(mean_array_error_prop4_7_4);
mean_array_error_method5_db=10*log10(mean_array_error_method5);
mean_array_error_method6_db=10*log10(mean_array_error_method6);
mean_LS_error_h_upper_db=10*log10(mean_LS_error_h_upper);
mean_LS_error_h_lower_db=10*log10(mean_LS_error_h_lower);
figure
hold on
plot(rhos_db,mean_array_error_prop4_7_4_db,'r-o','Markersize',8,'Linewidth',1.5)
plot(rhos_db,mean_array_error_method5_db,'b-+','Markersize',8,'Linewidth',1.5)
plot(rhos_db,mean_array_error_method6_db,'k-V','Markersize',8,'Linewidth',1.5)
plot(rhos_db,mean_LS_error_h_upper_db,'y-V','Markersize',8,'Linewidth',1.5)
plot(rhos_db,mean_LS_error_h_lower_db,'g-V','Markersize',8,'Linewidth',1.5)
xlabel('SNR(dB)');
ylabel('NMSE(dB)');
title('Channel estimation error');
legend('Proposed algorithm','Upper bound (Proposed algorithm)','Lower bound (Proposed algorithm)','Upper bound(LS)','Lower bound(LS)');
hold off
%% 数据储存dB
data_store1=zeros(L_SNR,5);%共5条线
data_store1(:,1)=rhos_db.';
data_store1(:,2)=mean_array_error_prop4_7_4_db.';
data_store1(:,3)=mean_array_error_method5_db.';
data_store1(:,4)=mean_array_error_method6_db.';
data_store1(:,5)=mean_LS_error_h_upper_db.';
data_store1(:,6)=mean_LS_error_h_lower_db.';

figure
hold on
plot(array_error_prop4_7_4,'r-o','Markersize',8,'Linewidth',1.5);
plot(array_error_method5,'b-V','Markersize',8,'Linewidth',1.5);
legend('prop','不调整');
hold off


