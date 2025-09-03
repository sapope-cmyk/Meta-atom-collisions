%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculates effective collsion parameters and resulting velocities for single meta-atom collision using two methods %%%
%%
%% tc - collision length
%% CMa - CM_a in the paper
%% CRa - CR_a in the paper
%% CMe - CM_e in the paper
%% CRe - CR_e in the paper
%% CMr - CM_r in the paper
%% CRr - CR_r in the paper
%% Dx1 - Velocity of sphere 1 at separation
%% Dx2 - Velocity of sphere 2 at separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Run "Model_v3.m" first

space_min = -3; 
space_max = 3; 
space_steps = 1000;
A = logspace(space_min,space_max,space_steps);
wn = sqrt(2*k/m);

t = [1e-10:1e-10:2*2*pi/wn];

tc_nr = 0.5*2*pi/wn;

tc = [];
Ddx = [];

Dx1 = [];
Dx2 = [];

CMa = [];
CRa = [];

AE_list = [0;0.5;1;1.5;2];
phir_list = [0:1/100:2]*pi;

vlr_l = [];
dlr_l = [];

Dx1_list = [];
Dx2_list = [];
tc_list = [];

CMa_list = [];
CRa_list = [];
CMe_list = [];
CRe_list = [];
CMr_list = [];
CRr_list = [];

CMa_a_list = [];
CRa_a_list = [];
CMe_a_list = [];
CRe_a_list = [];
CMr_a_list = [];
CRr_a_list = [];

vlr_list = [];
dlr_list = [];

Model_data = [];

Model_data_AE = [];

s = tf('s');

Mss = [m 0 0 0;0 mr 0 0;0 0 m 0;0 0 0 mr];

%% Sets up ss model x=[du u];
Ass11 = -inv(Mss)*zeros(4,4);
Ass21 = eye(4);
Ass22 = zeros(4,4);
Bss = [-inv(Mss)*[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];zeros(4,4)];
%% output vector y = [[u1 u2 u2-u1 u1+u2] [du1 du2 du2-du1 du1+du2] [ur1 ur2 ur2-ur1 ur1+ur2] [dur1 dur2 dur2-dur1 dur1+dur2]]'
Css = [[zeros(4,4) [1 0 0 0;0 0 1 0;-1 0 1 0;1 0 1 0];[1 0 0 0;0 0 1 0;-1 0 1 0;1 0 1 0] zeros(4,4)];[zeros(4,4) [0 1 0 0;0 0 0 1;0 -1 0 1;0 1 0 1];[0 1 0 0;0 0 0 1;0 -1 0 1;0 1 0 1] zeros(4,4)]];
Dss = zeros(16,4);

for AE_i = 1:length(AE_list)
    
    AE_i

    %% For loop to iterate through set of resonator initial phases
    for philr_i = 1:length(phir_list)
        
        philr_i
        
        %% For loop to iterate through set of resonator initial total energies
        for Ai = 1:length(A)

            %% Sets resonator displacement and velocity initial conditions based on total energy (KE + PE) constant
            vlr = sqrt(m*AE_list(AE_i)/(mr*A(Ai)*wn^2))*sqrt(A(Ai)*wn^2)*cos(phir_list(philr_i));
            vlr_l(Ai) = vlr;
            
            dlr = sqrt(m*AE_list(AE_i)/(mr*A(Ai)*wn^2))*sin(phir_list(philr_i));
            dlr_l(Ai) = dlr;

            %% Sets up coefficients 
            wr = sqrt(A(Ai)*wn^2);
            kr = mr*wr^2;
            Mes = tf([mr*kr],[mr 0 kr]);
            Me = m+Mes;
            gr = vlr+s*dlr;
            gall = Mes*gr+m*1;
            Dxs = -gall/(Me*s^2+2*k);
            Sxs = gall/(Me*s^2);

            Kss = -[-k-kr kr k 0;kr -kr 0 0;k 0 -k-kr kr;0 0 kr -kr];
            Ass = [Ass11 -inv(Mss)*Kss;Ass21 Ass22];
            Gss = ss(Ass,Bss,Css,Dss);

            x_init = [1;vlr;0;0;0;dlr;0;0];
            [dx,t_dx] = initial(Gss,x_init,t);
            t_ii = find((dx(3:end,3)./dx(2:end-1,3))<=0,1,'first');
            tc(Ai) = t_dx(t_ii); %% Separation time

            %% Calculates effective parameter values from simulation outputs of collision
            CMe_ss(Ai) = dx(t_ii,8);
            CRe_ss(Ai) = dx(t_ii,7);
            CMr_ss(Ai) = dx(t_ii,16);
            CRr_ss(Ai) = dx(t_ii,15);
            Dx1_ss = dx(t_ii,5);
            Dx2_ss = dx(t_ii,6);
            Ddx_ss = dx(t_ii,7);

            CMa_ss(Ai) = (m*CMe_ss(Ai)+mr*CMr_ss(Ai))/(m+mr);
            CRa_ss(Ai) = (m*CRe_ss(Ai)+mr*CRr_ss(Ai))/(m+mr);

            %% Calculates effective parameter values directly from TF model
            [CMe_CMa,CMa_t] = impulse(s*Sxs,tc(Ai));
            CMe(Ai) = CMe_CMa(end);
            [CRe_CMa,CRa_t] = impulse(s*Dxs,tc(Ai));
            CRe(Ai) = CRe_CMa(end);

            Dxrs = (kr*Dxs-mr*gr)/(mr*s^2+kr);
            Sxrs = (kr*Sxs+mr*gr)/(mr*s^2+kr);
            [CMr_CMa,CMa_t] = impulse(s*Sxrs,tc(Ai));
            CMr(Ai) = CMr_CMa(end);
            [CRr_CMa,CRa_t] = impulse(s*Dxrs,tc(Ai));
            CRr(Ai) = CRr_CMa(end);

            CMa(Ai) = (m*CMe(Ai)+mr*CMr(Ai))/(m+mr);
            CRa(Ai) = (m*CRe(Ai)+mr*CRr(Ai))/(m+mr);

            Dx1 = 0.5*(CMe(Ai)-CRe(Ai));
            Dx2 = 0.5*(CMe(Ai)+CRe(Ai));
            Ddx = Dx2-Dx1;            
           
        end

        Model_data(philr_i).tc = tc;
        Model_data(philr_i).Ddx = Ddx;

        Model_data(philr_i).CMa_list = CMa;
        Model_data(philr_i).CRa_list = CRa;
        Model_data(philr_i).CMe_list = CMe;
        Model_data(philr_i).CRe_list = CRe;
        Model_data(philr_i).CMr_list = CMr;
        Model_data(philr_i).CRr_list = CRr;
        Model_data(philr_i).Dx1 = Dx1;
        Model_data(philr_i).Dx2 = Dx2;

        Model_data(philr_i).CMa_ss_list = CMa_ss;
        Model_data(philr_i).CRa_ss_list = CRa_ss;
        Model_data(philr_i).CMe_ss_list = CMe_ss;
        Model_data(philr_i).CRe_ss_list = CRe_ss;
        Model_data(philr_i).CMr_ss_list = CMr_ss;
        Model_data(philr_i).CRr_ss_list = CRr_ss;
        Model_data(philr_i).Dx1_ss = Dx1_ss;
        Model_data(philr_i).Dx2_ss = Dx2_ss;
        
        Dx1_list(philr_i,:) = Dx1;
        Dx2_list(philr_i,:) = Dx2;
        tc_list(philr_i,:) = tc;
        CMa_list(philr_i,:) = CMa;
        CRa_list(philr_i,:) = CRa;
        CMe_list(philr_i,:) = CMe;
        CRe_list(philr_i,:) = CRe;
        CMr_list(philr_i,:) = CMr;
        CRr_list(philr_i,:) = CRr;

        CMa_ss_list(philr_i,:) = CMa_ss;
        CRa_ss_list(philr_i,:) = CRa_ss;
        CMe_ss_list(philr_i,:) = CMe_ss;
        CRe_ss_list(philr_i,:) = CRe_ss;
        CMr_ss_list(philr_i,:) = CMr_ss;
        CRr_ss_list(philr_i,:) = CRr_ss;

        vlr_list(philr_i,:) = vlr_l;
        dlr_list(philr_i,:) = dlr_l;

        
    end
    
    Model_data_AE(AE_i).Dx1_list = Dx1_list;
    Model_data_AE(AE_i).Dx2_list = Dx2_list;
    Model_data_AE(AE_i).tc_list = tc_list;
    Model_data_AE(AE_i).CMa_list = CMa_list;
    Model_data_AE(AE_i).CRa_list = CRa_list;
    Model_data_AE(AE_i).CMe_list = CMe_list;
    Model_data_AE(AE_i).CRe_list = CRe_list;
    Model_data_AE(AE_i).CMr_list = CMr_list;
    Model_data_AE(AE_i).CRr_list = CRr_list;

    Model_data_AE(AE_i).CMa_ss_list = CMa_ss_list;
    Model_data_AE(AE_i).CRa_ss_list = CRa_ss_list;
    Model_data_AE(AE_i).CMe_ss_list = CMe_ss_list;
    Model_data_AE(AE_i).CRe_ss_list = CRe_ss_list;
    Model_data_AE(AE_i).CMr_ss_list = CMr_ss_list;
    Model_data_AE(AE_i).CRr_ss_list = CRr_ss_list;
    Model_data_AE(AE_i).Dx1_ss_list = Dx1_ss;
    Model_data_AE(AE_i).Dx2_ss_list = Dx2_ss;
    
    Model_data_AE(AE_i).vlr_list = vlr_list;
    Model_data_AE(AE_i).dlr_list = dlr_list;

    
end
