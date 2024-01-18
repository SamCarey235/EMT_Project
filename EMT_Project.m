%% 

clc
clear

% Part A 

num_nodes = 3; % number of nodes
num_R = 2;
num_L = 0;
num_C = 0;
num_T = 1;
num_Vs = 1;
num_Is = 0; 
freq = 60;
pi = 3.141592653;

% resistor parameters
G = [1 1e-6];
G_node1 = [1 2];
G_node2 = [3 4];

% transformer parameters
T = 0.001;
T_node1 = 1;
T_node2 = 2;

% voltage source parameters 
V_mag = 35e3;
V_phase = 0;
V_node1 = 3;
V_node2 = 4;

% current source paramters
I_mag = [];
I_phase = [];
I_node1 =[];
I_node2 =[];

% create the admittance matrix 
Yorg = zeros(num_nodes + 1, num_nodes +1);
Ytemp = 0;

% add resitances to the matrix

for i = 1:num_R
    Ytemp = G(i);
    Yorg(G_node1(i), G_node1(i)) = Yorg(G_node1(i), G_node1(i)) + Ytemp;
    Yorg(G_node1(i), G_node2(i)) = Yorg(G_node1(i), G_node2(i)) - Ytemp;
    Yorg(G_node2(i), G_node1(i)) = Yorg(G_node2(i), G_node1(i)) - Ytemp;
    Yorg(G_node2(i), G_node2(i)) = Yorg(G_node2(i), G_node2(i)) + Ytemp;
end

% add transformers to the matrix

for i = 1:num_T
    Ytemp = 1 / (1j * 2 * pi * freq * T(i));
    Yorg(T_node1(i), T_node1(i)) = Yorg(T_node1(i), T_node1(i)) + Ytemp;
    Yorg(T_node1(i), T_node2(i)) = Yorg(T_node1(i), T_node2(i)) - Ytemp;
    Yorg(T_node2(i), T_node1(i)) = Yorg(T_node2(i), T_node1(i)) - Ytemp;
    Yorg(T_node2(i), T_node2(i)) = Yorg(T_node2(i), T_node2(i)) + Ytemp;
end

% add the voltage sources
Ve = zeros(num_Vs + 1, 1);

for i = 1:num_Vs
    if (V_node1(i) == num_nodes + 1) % if its ground node 
        V_node1(i) = V_node2(i); % use the other node
        V_node2(i) = num_nodes + 1; % set the other node as ground
    end 

    Ve(i) = V_mag(i) * exp(1j * V_phase(i));
%    Ve_nodes(i) = V_node1(i);

end

% add the current sources

I = zeros(num_nodes + 1, 1);
for i = 1:num_Is
    I(I_node1(i)) = I_mag(i) * exp(1j * I_phase);
    I(I_node2(i)) = -I(I_node1(i));
end   

Vd_nodes = [1, 2];
Ve_nodes = [3, 4];

Ydd = Yorg(Vd_nodes, Vd_nodes);
Id = I(Vd_nodes);
Yde = Yorg(Vd_nodes, Ve_nodes);

I = Id - Yde * Ve;
Vd = Ydd \ I;

% calculate branch currents

V1 = Vd(1);
V2 = Vd(2);
Itr = 0;

% results

V1 = real(V1);
V2 = real(V2);

%% Part B 

% transformer leakage current and turns ratio
Llk = 0.01;
a = 50;

% initialisation
Ts = 1e-5; % timestep
T = 0:Ts:0.1; % duration
tickmax = length(T); 

% store state variables in the circuit VHV, IHV, VLV, ILV
num_state_vars = 4; 

% store voltages on each dependent node
state_vars = zeros(tickmax, num_state_vars);
voltages = zeros(tickmax, num_nodes + 1);

% fill in the state varibales and voltages
state_vars(1,:) = [V1 0 V2 0]; % [VLV IHV VLV ILV]
voltages(1,:) = [V1 V2 35e3 0]; % [V1 V2 230e3 0]

% create the admittance matrix 
Yorg = zeros(num_nodes + 1, num_nodes +1);
Ytemp = 0;

% add resitances to the matrix

for i = 1:num_R
    Ytemp = G(i);
    Yorg(G_node1(i), G_node1(i)) = Yorg(G_node1(i), G_node1(i)) + Ytemp;
    Yorg(G_node1(i), G_node2(i)) = Yorg(G_node1(i), G_node2(i)) - Ytemp;
    Yorg(G_node2(i), G_node1(i)) = Yorg(G_node2(i), G_node1(i)) - Ytemp;
    Yorg(G_node2(i), G_node2(i)) = Yorg(G_node2(i), G_node2(i)) + Ytemp;
end

% add transformers to the matrix

for i = 1:num_T
    Ytemp = Ts*0.5/T(i);
    Yorg(T_node1(i), T_node1(i)) = Yorg(T_node1(i), T_node1(i)) + Ytemp;
    Yorg(T_node1(i), T_node2(i)) = Yorg(T_node1(i), T_node2(i)) - Ytemp;
    Yorg(T_node2(i), T_node1(i)) = Yorg(T_node2(i), T_node1(i)) - Ytemp;
    Yorg(T_node2(i), T_node2(i)) = Yorg(T_node2(i), T_node2(i)) + Ytemp;
end

% start the loop

for k = 2:tickmax 
%k = 2;
    % update voltage source
    Ve = zeros(num_Vs + 1, 1);
    Ve(1) = V_mag * cos(2 * pi * freq * k * Ts + V_phase);

    % update current source 
    I = zeros(num_nodes + 1, 1);
    
    % update equivalent current source related to the transformer
    Ytemp = Ts * (0.5 / T(i));
    Itemp = Ytemp *(state_vars(k-1, 1) - a*(state_vars(k-1, 3)) + state_vars(k-1, 2));
    I(1) = I(1) - Itemp;
    I(2) = I(2) + a*Itemp;
    %I(4) = I(4) + Itemp - a*Itemp;

    % eliminate the excitation nodes
    
    Ydd = Yorg(Vd_nodes, Vd_nodes);
    Id = I(Vd_nodes);
    Yde = Yorg(Vd_nodes, Ve_nodes);

    I = Id - Yde * Ve;

    % dependent node voltage
    Vd = Ydd \ I;

    % second node voltages
    for i = 1:length(Vd_nodes)
        voltages(k, Vd_nodes(i)) = Vd(i);
    end
    
    % add in the excitation and ground voltages again
    voltages(k, 3) = Ve(1); % excitation voltage from source 
    voltages(k, 4) = 0; % ground voltage 

    % calculate state varibales
    
    % transformers
    for i = 1:num_T
        ItempHV = (0.5 * Ts/T(i)) * state_vars(k-1, 1) + state_vars(k-1, 2);
        ItempLV = (0.5 * Ts/T(i)) * state_vars(k-1, 3) + state_vars(k-1, 4);

        % VHV
        state_vars(k, 1) = voltages(k, 1);
        % VLV
        state_vars(k, 3) = voltages(k, 2);

        % IHV
        state_vars(k, 2) = (0.5 * Ts/T(i)) * state_vars(k, 1) + ItempHV;
        % ILV
        state_vars(k, 4) = (0.5 * Ts/T(i)) * state_vars(k, 3) + ItempLV;

    end

end

figure;
plot(T,state_vars(:,1))
xlabel('Time in s') 
ylabel('Vl in V') 
figure;
plot(T,state_vars(:,2))
xlabel('Time in s') 
ylabel('Il in A') 
figure;
plot(T,state_vars(:,3))
xlabel('Time in s') 
ylabel('Vc in V') 
figure;
plot(T,state_vars(:,4))
xlabel('Time in s') 
ylabel('Ic in A') 
