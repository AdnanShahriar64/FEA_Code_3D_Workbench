clear all
clc
 
addpath('Input_Files','Element_Library','Mesh_Generation_Library','System_Preprocessor_Files','Decomposition_Files','Runtime_Calculation_Files','Solver_Files','Post_Processing_Files','External_Subsystem_Calculation')
addpath('System_Preprocessor_Files\Solid')
addpath('System_Preprocessor_Files\Fluid')
addpath('Runtime_Calculation_Files\Solid')
addpath('Runtime_Calculation_Files\Fluid')
addpath('Runtime_Calculation_Files\Fluid_Thermal');
addpath('System_Preprocessor_Files\Thermal');
addpath('Runtime_Calculation_Files\Thermal');
addpath('Runtime_Calculation_Files\Fluid_Pressure');


fprintf('\nInitial System building starts---------------\n')

%% Input FE elements
    [Element_Information]=Element_Loader;

%% Model Input File
% Read Abaqus Input file
      
      Name='Tunnel_Test_Vel.inp';
      S_Names={'Pr_Z','Vel_Z'};
      Reg_Name={'Sx','S1'};
      BC_Name={'Zero_V'};%,'Zero_V_Y','Zero_V_Z'};

      [Mesh_Information_V]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information);
      % [Mesh_Information_V]=Element_Converter_20_to_75(Mesh_Information_V);
      % Mesh_Information_V.Nodes(:,4)=Mesh_Information_V.Nodes(:,4)*0.12/0.025;

      [Mesh_Information_V]=Geometry_Model_Parameters_Preprocessor_Fluid(Mesh_Information_V,Element_Information);
      Mesh_Information_V.Boundary_Cond=[1;1;1];
      
      
%% Read Abaqus Input file
      
      Name='Tunnel_Test_Pr.inp';
      S_Names={'Pr_Z','Vel_Z'};
      Reg_Name={'Sx','S1'};
      BC_Name={'Zero_V_X','Zero_V_Y','Zero_V_Z'};

      [Mesh_Information_PT]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information);
      % [Mesh_Information_PT]=Element_Converter_20_to_75(Mesh_Information_PT);
      % Mesh_Information_PT.Nodes(:,4)=Mesh_Information_PT.Nodes(:,4)*0.12/0.025;

      [Mesh_Information_PT]=Geometry_Model_Parameters_Preprocessor_Fluid(Mesh_Information_PT,Element_Information);
      Mesh_Information_PT.Boundary_Cond=[1;1;1];      

%% Material input      
      % Material_Information.MoE=0+68*10^9;
      % Material_Information.PoRat=0.3;
      Material_Information.rho=2;
      Material_Information.beta=10^-6;
      Material_Information.miu=1.48*10^-5;
      Material_Information.Cp=1000;

     %% Rotation
%      Alpha=-40*pi/180;
%      Beta=30*pi/180;
%      Gamma=80*pi/180;
%      
%      Rz=[cos(Gamma) -sin(Gamma) 0;sin(Gamma) cos(Gamma) 0; 0 0 1];
%      Ry=[cos(Beta) 0 sin(Beta);0 1 0;-sin(Beta) 0 cos(Beta)];
%      Rx=[1 0 0; 0 cos(Alpha) -sin(Alpha);0 sin(Alpha) cos(Alpha)];
%      R1=Rx*Ry*Rz;
%      
%      for i=1:length(Mesh_Information.Nodes)
%          Mesh_Information.Nodes(i,2:4)=inv(R1)*[Mesh_Information.Nodes(i,2);Mesh_Information.Nodes(i,3);Mesh_Information.Nodes(i,4)];
%      end
     
%% System preprocessor
    
    % [FE_Information_Region]=Region_Integration_3D_Solid(Element_Information,Mesh_Information);

    [FE_Information_Region_Velocity]=Region_Integration_3D_Fluid(Element_Information,Mesh_Information_V);
    [FE_Information_Region_Thermal]=Region_Integration_3D_Heat_Transfer(Element_Information,Mesh_Information_PT);
    [FE_Information_Surface_Thermal]=Surface_Integration_3D_Flux_Thermal(Element_Information,Mesh_Information_PT);    

    % [FE_Information_Surface]=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information);
    
    fprintf('Initial System building completed--------------\n')    
    
    
%% Only Fluid

% Diffisive matrix
[K11e,K21e,K31e,K12e,K22e,K32e,K13e,K23e,K33e]=Elemental_Diffisuive_Matrix_Generator(Element_Information,Mesh_Information_V,Material_Information,FE_Information_Region_Velocity);
[K11] = Assembly(K11e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K21] = Assembly(K21e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K31] = Assembly(K31e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K12] = Assembly(K12e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K22] = Assembly(K22e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K32] = Assembly(K32e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K13] = Assembly(K13e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K23] = Assembly(K23e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[K33] = Assembly(K33e,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);


%% Convective matrix




%% Fluid-Pressure
%[Q1e,Q2e,Q3e]=Elemental_Pressure_Matrix_Generator(Element_Information,Mesh_Information,FE_Information_Region);
[Q1e,Q2e,Q3e]=Elemental_Pressure_Matrix_Generator(Element_Information,Mesh_Information_PT,Mesh_Information_V,FE_Information_Region_Velocity);

[Q1] = Assembly_Multi_Model(Q1e,Mesh_Information_V.Type,Mesh_Information_PT.Type,1,Mesh_Information_V.IEN,Mesh_Information_PT.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Mesh_Information_PT.nnp_Model,Element_Information.nnp_Elemental);
[Q2] = Assembly_Multi_Model(Q2e,Mesh_Information_V.Type,Mesh_Information_PT.Type,1,Mesh_Information_V.IEN,Mesh_Information_PT.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Mesh_Information_PT.nnp_Model,Element_Information.nnp_Elemental);
[Q3] = Assembly_Multi_Model(Q3e,Mesh_Information_V.Type,Mesh_Information_PT.Type,1,Mesh_Information_V.IEN,Mesh_Information_PT.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Mesh_Information_PT.nnp_Model,Element_Information.nnp_Elemental);


K=[ 2*K11+K22+K33   K21              K31                 -Q1;
    K12             K11+2*K22+K33    K32                 -Q2;    
    K13             K23              K11+K22+2*K33       -Q3;
   -Q1'            -Q2'             -Q3'                 -Q3'*Q3*0];    
    


    %SpcH=900; % 
    %Density=2700;
    
    Type=Mesh_Information_V.Type;
    nel=Mesh_Information_V.nel_Model;    
    
    [MeF]=Elemental_Capacitance_Matrix_Generator(1,Material_Information.rho,Element_Information.Shape_Elemental,FE_Information_Region_Velocity.dJB,Element_Information.W_Elemental,Type,1:nel,nel);
    [MF1] = Assembly(MeF,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
%     for i=1:length(Mesh_Information_V.S_Nodes)
%         Node_No=Mesh_Information_V.S_Nodes(i,1);
%         MF1(:,Node_No)=MF1(:,Node_No)*0;
%         MF1(Node_No,:)=MF1(Node_No,:)*0;
%         MF1(Node_No,Node_No)=1;
%     end      
    
    
% Body force vector
[F_1e,F_2e,F_3e]=Fluid_Elemental_Body_Force_Vector_Calculator(Element_Information,Mesh_Information_V,FE_Information_Region_Velocity,Material_Information);
[F1]=Assembly_Force(F_1e,Mesh_Information_V.IEN,Mesh_Information_V.Type,1,Mesh_Information_V.nnp_Model,Mesh_Information_V.nel_Model,Element_Information.nnp_Elemental);
[F2]=Assembly_Force(F_2e,Mesh_Information_V.IEN,Mesh_Information_V.Type,1,Mesh_Information_V.nnp_Model,Mesh_Information_V.nel_Model,Element_Information.nnp_Elemental);
[F3]=Assembly_Force(F_3e,Mesh_Information_V.IEN,Mesh_Information_V.Type,1,Mesh_Information_V.nnp_Model,Mesh_Information_V.nel_Model,Element_Information.nnp_Elemental);

FF_B=K(:,1)*0;
FF_B(1:length(F1)*3)=[F1*0;F2;F3*0]; 

    
%% Thermal-Fluid

[F_M1e,F_M2e,F_M3e]=Fluid_Expansion_Matrix_Calculator(Element_Information,Mesh_Information_V,Mesh_Information_PT,FE_Information_Region_Velocity,Material_Information);

[F1_M] = Assembly_Multi_Model(F_M2e,Mesh_Information_V.Type,Mesh_Information_PT.Type,1,Mesh_Information_V.IEN,Mesh_Information_PT.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Mesh_Information_PT.nnp_Model,Element_Information.nnp_Elemental);
[F2_M] = Assembly_Multi_Model(F_M2e,Mesh_Information_V.Type,Mesh_Information_PT.Type,1,Mesh_Information_V.IEN,Mesh_Information_PT.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Mesh_Information_PT.nnp_Model,Element_Information.nnp_Elemental);
[F3_M] = Assembly_Multi_Model(F_M3e,Mesh_Information_V.Type,Mesh_Information_PT.Type,1,Mesh_Information_V.IEN,Mesh_Information_PT.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Mesh_Information_PT.nnp_Model,Element_Information.nnp_Elemental);

% % BC with zero velocity
%     FV_G=K(:,1)*0+1;
%     for i=1:length(Mesh_Information_V.S_Nodes)
%         Node_No=Mesh_Information_V.S_Nodes(i,1);
%         F1_M(Node_No)=0;
%         F2_M(Node_No)=0;
%         F3_M(Node_No)=0;
%         
%         FV_G(Node_No)=0;
%         FV_G(Mesh_Information_V.nnp_Model+Node_No)=0;
%         FV_G(Mesh_Information_V.nnp_Model*2+Node_No)=0;
%     end


FF_Ex=zeros(length(K),size(F1_M,2)*4);
FF_Ex(1:size(F1_M,1)*3,1:size(F1_M,2)*3)=[F1_M*0   F2_M*0  F3_M*0
                                          F1_M*0   F2_M    F3_M*0
                                          F1_M*0   F2_M*0  F3_M*0  ];
 
 
     % F1_M'*0   F2_M'*0  F3_M'*0   zeros(size(F3_M,2),size(F3_M,2))];

% [De]=Elemental_Convective_Matrix_Generator_Fluid_Thermal(Element_Information,Mesh_Information_PT,Mesh_Information_V,Material_Information,FE_Information_Region_Thermal);     
% [D] = Assembly(De,Mesh_Information_PT.Type,1,Mesh_Information_PT.IEN,Mesh_Information_PT.nel_Model,Mesh_Information_PT.nnp_Model,Element_Information.nnp_Elemental);


%% Only Thermal 

    Kx=30;%239;
    
    [Dm]=Material_Array_Generator_Heat_Transfer(Kx,Kx,Kx);  

% Region calculation
    Type=Mesh_Information_PT.Type;
    nel=Mesh_Information_PT.nel_Model;
    [Ke_Conductive]=Elemental_Conductivity_Matrix_Generator(FE_Information_Region_Thermal.B,Dm,FE_Information_Region_Thermal.dJB,Element_Information.W_Elemental,Type,1:nel,nel);   
    
    % SpcH=900; % 
    % Density=2700;
    [Me]=Elemental_Capacitance_Matrix_Generator(Material_Information.Cp,Material_Information.rho,Element_Information.Shape_Elemental,FE_Information_Region_Thermal.dJB,Element_Information.W_Elemental,Type,1:nel,nel);
    
% Convection calculation    
    h=0;
    Ke_surf_per_h=FE_Information_Surface_Thermal.Ke_surf;
    Ke_Convective=h*Ke_surf_per_h(:,:,:,1);
    
% Radiation calculation
    epsilon=0.09;
    sigma=5.67*10^-8;
    Ke_Radiation=epsilon*sigma*(Ke_surf_per_h(:,:,:,1)*0+Ke_surf_per_h(:,:,:,1)*1); % Internal surface omitted
    
    sum(sum(sum(sum(Ke_Radiation))))
    
% Total matrix    
    Le=Ke_Conductive+Ke_Convective;
    
% Assembly    
    [L] = Assembly_Thermal(Le,Mesh_Information_PT.Type,1,Mesh_Information_PT.IEN,nel,Mesh_Information_PT.nnp_Model);    
    [M] = Assembly_Thermal(Me,Mesh_Information_PT.Type,1,Mesh_Information_PT.IEN,nel,Mesh_Information_PT.nnp_Model);
    
    M_artificial=1/Material_Information.rho*M;
    
%% Assign boundary cond

%     FTK=zeros(Mesh_Information_PT.nnp_Model,length(Mesh_Information_PT.S_Nodes));
%     FTM=zeros(Mesh_Information_PT.nnp_Model,length(Mesh_Information_PT.S_Nodes));
%     for i=1:length(Mesh_Information_PT.S_Nodes)
%         
%         Node_No=Mesh_Information_PT.S_Nodes(i,1);
%         
%         FTK(:,i)=-L(:,Node_No);
%         FTM(:,i)=-M(:,Node_No);
%         
%         L(:,Node_No)=L(:,Node_No)*0;
%         L(Node_No,Node_No)=1;
%         
%         M(:,Node_No)=M(:,Node_No)*0;
%         M(Node_No,Node_No)=1;
%     
%     end

    
%% bACK TO FLUID

MF=[MF1      MF1*0      MF1*0       MF1(:,1)*M_artificial(1,:)*0
    MF1*0    MF1        MF1*0       MF1(:,1)*M_artificial(1,:)*0
    MF1*0    MF1*0      MF1         MF1(:,1)*M_artificial(1,:)*0
    M_artificial(:,1)*MF1(1,:)*0    M_artificial(:,1)*MF1(1,:)*0      M_artificial(:,1)*MF1(1,:)*0       M_artificial*-1];


%%

UP=K(:,1)*0+0*101325;
UP(1:Mesh_Information_V.nnp_Model*3)=0;

dt=0.01;
Theta=2/3;

F=FF_B*-1;

% S_Nodes=Mesh_Information_V.S_Nodes;
% Nodes=Mesh_Information_V.Nodes;
% 
% cnt=1;
% for n=1:length(S_Nodes)
%     if abs(Nodes(S_Nodes(n),3))<10^-10
%         Chosen(cnt)=n;
%         cnt=cnt+1;
%     end
% end
% 
% UnChosen=setdiff(1:length(Mesh_Information_V.S_Nodes),Chosen);
% 
% 
% S_Nodes=Mesh_Information_PT.S_Nodes;
% Nodes=Mesh_Information_PT.Nodes;
% 
% cnt=1;
% for n=1:length(S_Nodes)
%     if abs(Nodes(S_Nodes(n),3))<10^-10
%         Chosen2(cnt)=n;
%         cnt=cnt+1;
%     end
% end
% 
% UnChosen2=setdiff(1:length(Mesh_Information_PT.S_Nodes),Chosen2);
% dsa

% BC_Sets=setdiff(Mesh_Information_V.BC_Nodes,Mesh_Information_V.S_Nodes(:,1));
% BC_Sets=setdiff(BC_Sets,Mesh_Information_V.S_Nodes(:,2));

BC_Sets=Mesh_Information_V.BC_Nodes;

Vel_Nodes=Mesh_Information_V.S_Nodes(:,2);
Pr_Nodes=Mesh_Information_PT.S_Nodes(:,2);

% Vel_Nodes=setdiff(Mesh_Information_V.S_Nodes(:,2),Mesh_Information_V.BC_Nodes);
% Pr_Nodes=Mesh_Information_PT.S_Nodes(:,1);%setdiff(Mesh_Information_PT.S_Nodes(:,1),Mesh_Information_PT.BC_Nodes);

MFO=MF;

% BC_Sets_OZ


%%



for It=1:100

UP_Del=UP;
% if It>1
%    UP_Del(Zero_Ass)=0;
%    UP_Del(Ten_Ass)=10;
% end
    
    
Velocity_x=UP_Del(1:Mesh_Information_V.nnp_Model);
Velocity_y=UP_Del(Mesh_Information_V.nnp_Model+1:Mesh_Information_V.nnp_Model*2);
Velocity_z=UP_Del(Mesh_Information_V.nnp_Model*2+1:Mesh_Information_V.nnp_Model*3);

[Ce]=Elemental_Convective_Matrix_Generator(Element_Information,Mesh_Information_V,Material_Information,FE_Information_Region_Velocity,Velocity_x,Velocity_y,Velocity_z);
[Ci] = Assembly(Ce,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
 
C0=K*0;
C0(1:length(Ci)*3,1:length(Ci)*3)=[Ci       Ci*0        Ci*0 
                                   Ci*0     Ci          Ci*0        
                                   Ci*0     Ci*0        Ci  ];
                                         
[dR1due,dR1dve,dR1dwe,dR2due,dR2dve,dR2dwe,dR3due,dR3dve,dR3dwe]=Elemental_Convective_Matrix_Generator_Diff(Element_Information,Mesh_Information_V,Material_Information,FE_Information_Region_Velocity,Velocity_x,Velocity_y,Velocity_z);
[dR1du] = Assembly(dR1due,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR1dv] = Assembly(dR1dve,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR1dw] = Assembly(dR1dwe,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR2du] = Assembly(dR2due,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR2dv] = Assembly(dR2dve,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR2dw] = Assembly(dR2dwe,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR3du] = Assembly(dR3due,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR3dv] = Assembly(dR3dve,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
[dR3dw] = Assembly(dR3dwe,Mesh_Information_V.Type,1,Mesh_Information_V.IEN,Mesh_Information_V.nel_Model,Mesh_Information_V.nnp_Model,Element_Information.nnp_Elemental);
                              
C1=K*0;
C1(1:length(Ci)*3,1:length(Ci)*3)=[dR1du dR1dv dR1dw;
                                   dR2du dR2dv dR2dw 
                                   dR3du dR3dv dR3dw];
                                                          

KF=C0+K;                               
Jac=C0+C1+K;                               


F=K(:,1)*0;
FV_G=F+1;

count=1;
for i=1:length(BC_Sets)
        % i=Chosen(sl);
        Node_No=BC_Sets(i,1);
        
        if Node_No>0
            X_Node_No=Node_No;
            Y_Node_No=Mesh_Information_V.nnp_Model+Node_No;
            Z_Node_No=Mesh_Information_V.nnp_Model*2+Node_No;
            
            % Jac(:,X_Node_No)=Jac(:,X_Node_No)*0;
            Jac(X_Node_No,:)=Jac(X_Node_No,:)*0;
            Jac(X_Node_No,X_Node_No)=1;  
            
            % KF(:,X_Node_No)=KF(:,X_Node_No)*0;
            KF(X_Node_No,:)=KF(X_Node_No,:)*0;
            KF(X_Node_No,X_Node_No)=1;   
            
            % Jac(:,X_Node_No)=Jac(:,X_Node_No)*0;
            Jac(Y_Node_No,:)=Jac(Y_Node_No,:)*0;
            Jac(Y_Node_No,Y_Node_No)=1;  
            
            % KF(:,X_Node_No)=KF(:,X_Node_No)*0;
            KF(Y_Node_No,:)=KF(Y_Node_No,:)*0;
            KF(Y_Node_No,Y_Node_No)=1;  
            
            
            % Jac(:,X_Node_No)=Jac(:,X_Node_No)*0;
            Jac(Z_Node_No,:)=Jac(Z_Node_No,:)*0;
            Jac(Z_Node_No,Z_Node_No)=1;  
            
            % KF(:,X_Node_No)=KF(:,X_Node_No)*0;
            KF(Z_Node_No,:)=KF(Z_Node_No,:)*0;
            KF(Z_Node_No,Z_Node_No)=1;              
            
            
            FV_G(X_Node_No)=0;
            
            Zero_Ass(count,1)=X_Node_No;
            count=count+1;
            
        end
end 

% for i=1:length(BC_Sets)
%         % i=Chosen(sl);
%         Node_No=BC_Sets(i,2);
%        
%         if Node_No>0        
%             X_Node_No=Node_No;
%             Y_Node_No=Mesh_Information_V.nnp_Model+Node_No;
%             Z_Node_No=Mesh_Information_V.nnp_Model*2+Node_No;
% 
%             % Jac(:,Y_Node_No)=Jac(:,Y_Node_No)*0;
%             Jac(Y_Node_No,:)=Jac(Y_Node_No,:)*0;
%             Jac(Y_Node_No,Y_Node_No)=1;
% 
%             % KF(:,Y_Node_No)=KF(:,Y_Node_No)*0;
%             KF(Y_Node_No,:)=KF(Y_Node_No,:)*0;
%             KF(Y_Node_No,Y_Node_No)=1;
%             
%             FV_G(Y_Node_No)=0;
%             
%             Zero_Ass(count,1)=Y_Node_No;
%             count=count+1;
%         end
% end 
% 
% for i=1:length(BC_Sets)
%         % i=Chosen(sl);
%         Node_No=BC_Sets(i,3);
%         
%         if Node_No>0        
%             X_Node_No=Node_No;
%             Y_Node_No=Mesh_Information_V.nnp_Model+Node_No;
%             Z_Node_No=Mesh_Information_V.nnp_Model*2+Node_No;
% 
%             % Jac(:,Z_Node_No)=Jac(:,Z_Node_No)*0;
%             Jac(Z_Node_No,:)=Jac(Z_Node_No,:)*0;
%             Jac(Z_Node_No,Z_Node_No)=1;   
%             
%             % KF(:,Z_Node_No)=KF(:,Z_Node_No)*0;
%             KF(Z_Node_No,:)=KF(Z_Node_No,:)*0;
%             KF(Z_Node_No,Z_Node_No)=1;            
% 
%             FV_G(Z_Node_No)=0;
%             
%             Zero_Ass(count,1)=Z_Node_No;
%             count=count+1;
%         end
% end 


%
countVal=1;
for i=1:length(Vel_Nodes)
        % i=Chosen(sl);
        Node_No=Vel_Nodes(i);%Mesh_Information_V.S_Nodes(i,2);
        
        X_Node_No=Node_No;
        Y_Node_No=Mesh_Information_V.nnp_Model+Node_No;
        Z_Node_No=Mesh_Information_V.nnp_Model*2+Node_No;
        
        Vel_Fz=-KF(:,Z_Node_No);
        % Jac(:,Z_Node_No)=Jac(:,Z_Node_No)*0;
        Jac(Z_Node_No,:)=Jac(Z_Node_No,:)*0;
        Jac(Z_Node_No,Z_Node_No)=1;  

        % KF(:,Z_Node_No)=KF(:,Z_Node_No)*0;
        KF(Z_Node_No,:)=KF(Z_Node_No,:)*0;
        KF(Z_Node_No,Z_Node_No)=1;        
        
        % F(Z_Node_No)=-K(Z_Node_No,Z_Node_No)*10;
        % F=F+Vel_Fz*10;
        F(Z_Node_No)=10;
        Ten_Ass(countVal)=Z_Node_No;
        countVal=countVal+1;
end   


UPOLD=UP;
Residual=(KF*UP)-F;
UP=UP-Jac\Residual*10^-2;

Error=norm(UPOLD-UP);
RelError=Error/norm(UP)*100;
fprintf('Iteration=%d   Error=%d    Relative error=%d    Velocity norm=%d\n',It,Error,RelError,norm(UP))
UPS(:,It)=UP;
DSA


end

