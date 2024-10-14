clear all
clc
 
addpath('Input_Files','Element_Library','Mesh_Generation_Library','System_Preprocessor_Files','Decomposition_Files','Runtime_Calculation_Files','Solver_Files','Post_Processing_Files','External_Subsystem_Calculation')
addpath('System_Preprocessor_Files\Solid')
addpath('Runtime_Calculation_Files\Solid')
addpath('Abaqus_Input_Files')

fprintf('\nInitial System building starts---------------\n')

%% Input FE elements
    [Element_Information]=Element_Loader;

%% Model Input File
% Read Abaqus Input file

      Name='SPL_Section.inp';
      S_Names={'S4','S2','S3'};
      Reg_Name={'Sx','S1'};
      BC_Name={'S3'};

      [Mesh_Information]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information);
     

%%

  
% Assign proper dimension
    Outer_Radius=2.7;
    Inner_Radius=2.6;
    Middle_Radius=(Outer_Radius+Inner_Radius)/2;

    In_O_R=2.9;
    In_I_R=2.5;
    In_M_R=2.7;
    
    Nodes=Mesh_Information.Nodes;
    
    for i=1:length(Nodes)
        
       if Nodes(i,3)>=0
          Radius=norm(Nodes(i,2:4)); 
           if abs(Radius-In_O_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Outer_Radius;
           elseif abs(Radius-In_I_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Inner_Radius;
           elseif abs(Radius-In_M_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Middle_Radius;          
           else
               problem
           end
           
       else
           
          Radius=norm(Nodes(i,[2 4]));
          
          Z=abs(Nodes(i,3));
          Zmin=abs(min(Nodes(:,3)));
          Mult=1+Z/Zmin*0;
          
          
           if abs(Radius-In_O_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Outer_Radius*Mult;
           elseif abs(Radius-In_I_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Inner_Radius;
           elseif abs(Radius-In_M_R)<10^-4
               Nodes(i,2:4)=Nodes(i,2:4)/Radius*Middle_Radius;          
           else
               problem
           end          
       end
       
       

       
       
    end
    Mesh_Information.Nodes=Nodes;

    %% Rotation
     Alpha=-90*pi/180;
     Beta=0;
     Gamma=0;
     
     Rz=[cos(Gamma) -sin(Gamma) 0;sin(Gamma) cos(Gamma) 0; 0 0 1];
     Ry=[cos(Beta) 0 sin(Beta);0 1 0;-sin(Beta) 0 cos(Beta)];
     Rx=[1 0 0; 0 cos(Alpha) -sin(Alpha);0 sin(Alpha) cos(Alpha)];
     R1=Rx*Ry*Rz;
     
     for i=1:length(Nodes)
         Nodes(i,2:4)=inv(R1)*[Nodes(i,2);Nodes(i,3);Nodes(i,4)];
     end    
    
    
    Mesh_Information.Nodes=Nodes;
   
    
   
    
%%    
      [Mesh_Information]=Element_Converter_20_to_27(Mesh_Information,Element_Information);
      % [Mesh_Information]=Element_Converter_20_to_Transition_75_from_27(Mesh_Information,Element_Information);
      
      
      % [Mesh_Information]=Element_Converter_20_to_75(Mesh_Information);
      
      %Mesh_Information.Nodes(:,2)=Mesh_Information.Nodes(:,2)*01/0.3;
      %Mesh_Information.Nodes(:,3)=Mesh_Information.Nodes(:,3)*01/0.3;
      %Mesh_Information.Nodes(:,4)=Mesh_Information.Nodes(:,4)*01/0.025;

      [Mesh_Information]=Geometry_Model_Parameters_Preprocessor_Solid(Mesh_Information,Element_Information);
      Mesh_Information.Boundary_Cond=[1;1;1];
    
    
 
      [Mesh_Information]=Node_Rearranger(Element_Information,Mesh_Information,Mesh_Information.S_Nodes(:,1)); % Outer layer goes 1st
      [Mesh_Information]=Node_Rearranger(Element_Information,Mesh_Information,Mesh_Information.S_Nodes(:,2));
      [Mesh_Information]=Node_Rearranger(Element_Information,Mesh_Information,Mesh_Information.S_Nodes(:,3));
      
      Mesh_Information.BC_Nodes=Mesh_Information.S_Nodes(:,3);
          
    
%% Material input      
      Material_Information.MoE=68*10^9;
      Material_Information.PoRat=0.3;
      Material_Information.rho=2703;
      Material_Information.alpe=10^-6;

      [Dm]=Material_Array_Generator(Material_Information);
      
%% System preprocessor
    
    [FE_Information_Region]=Region_Integration_3D_Solid(Element_Information,Mesh_Information);

    [FE_Information_Surface]=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information);
    
    [IS_F]=Internal_Stress_to_Nodal_FV(Element_Information,Mesh_Information,FE_Information_Region,Dm);
    
    fprintf('Initial System building completed--------------\n')    
    
    
%% Mat property and Stiffness    
 
      LM=Mesh_Information.LM_Model;
      Target_BC_Node=min(setdiff(Mesh_Information.BC_Nodes,0));
      BC_x=Mesh_Information.BC_Nodes*3-2;
      BC_y=Mesh_Information.BC_Nodes*3-1;
      BC_z=Mesh_Information.BC_Nodes*3-0;
      
      for e=1:size(LM,2) 
          for n=1:size(LM,1)              
              if sum(intersect(LM(n,e),BC_x))>0
                  LM(n,e)=Target_BC_Node*3-2;
              elseif sum(intersect(LM(n,e),BC_y))>0
                  LM(n,e)=Target_BC_Node*3-1;
              elseif sum(intersect(LM(n,e),BC_z))>0
                  LM(n,e)=Target_BC_Node*3-0;
              end              
          end
      end
      neq=Target_BC_Node*3;
      Mesh_Information.LM_Model=LM;
      Mesh_Information.neq_Model_original=Mesh_Information.neq_Model;
      Mesh_Information.neq_Model=neq;
    
    [M,K]=Stiffness_Mass_Matrix_Assembled_Basic(Element_Information,Mesh_Information,Material_Information,FE_Information_Region,Dm);

    
    MO=M;
    KO=K;
    
    zxdcxz
%% Boundary input  
   
    [Bc]=Boundary_Condition_Set(Mesh_Information);
    f_BC=Bc*-1+1;
    
    % [K]=Boundary_Set(Bc,length(K),K);
    % [M]=Boundary_Set(Bc,length(M),M);
    
    [K,F_K]=Boundary_Set_Conserved(Bc,length(K),K);
    [M,F_M]=Boundary_Set_Conserved(Bc,length(M),M);
    
    K=K(1:1353,1:1353);
    M=M(1:1353,1:1353);
    F_K=F_K(1:1353,:);
    
    
    Mesh_Information.neq_Model=1353;
    
    % M(1354:1356,1354:1356)=0;
    

%% Pregen

    

    Tm=Mesh_Information.gp_X*0;
    Tm(37,12)=1;
    Impacted_Volume=FE_Information_Region.dJB(37,12);


   [f_Volumetric]=Internal_Stress_to_Equivalent_Nodal_Force(Element_Information,Mesh_Information,IS_F,Tm);

   Ps=Mesh_Information.S_Type(1,:)*100;
   
   [fS]=Pressure_to_Equivalent_Nodal_Force(Element_Information,Mesh_Information,FE_Information_Surface,Ps);

xcz
%% Dynamic

            dt=1*10^-4;%1E-05;
            delta=0.5; alpha=0.25;
            % gama    % beta
            a0=1/alpha/dt^2; a1=delta/alpha/dt;  
            a2=1/alpha/dt; 
            a3=1/2/alpha-1; 
            a4=delta/alpha-1;
            a5=dt/2*(delta/alpha-2); 
            a6=dt*(1-delta); 
            a7=delta*dt;

            C=K*0.00001;
            Kbar=K+a0*M+a1*C;
            ext_F=fS*0;
            IKbar=inv(Kbar);

%
           sigma_gp=zeros(6,size(Mesh_Information.IEN,1),size(Mesh_Information.IEN,2));
           strain_gp=sigma_gp;
           P_Eq=Mesh_Information.IEN*0;
            
            Nodes=Mesh_Information.Nodes;
            IEN=Mesh_Information.IEN;
            neq=Mesh_Information.neq_Model;

            U0=zeros(neq,1);Vel0=U0;Acc0=U0;Accal=U0;
            

            
            U_Store=U0;
            Acc_Store=U0;
            Vel_Store=U0;
            nn=IEN(16,12);%IEN(71,13);
            Time=0;
            
            Impactor_Velocity=2000;
            Impactor_Radius=50.E-03/2;
            Impactor_Density=2700;
            
            Impactor_Volume_Total=4/3*pi*Impactor_Radius^3;
            Impactor_Energy=0.5*Impactor_Volume_Total*Impactor_Density*Impactor_Velocity^2;
            F_Pressure=f_Volumetric/1.7000e+11;
            Pressure=0.5*2700*2000^2;
            Impact_Time=2*(50.E-03/2)/2000;     
            
            Pressure=0.5*Impactor_Density*Impactor_Velocity^2;
            
            
            

            count=1;
            Impactor_Volume_Old=0;
            for t=0:dt:50000*dt
                
                Impactor_Volume0=1/3*pi*Impactor_Velocity^2*t^2*(3*Impactor_Radius-Impactor_Velocity*t);
                Impactor_Volume1=1/3*pi*Impactor_Velocity^2*(t+dt)^2*(3*Impactor_Radius-Impactor_Velocity*(t+dt));
                
                dImpactor_Volume=(Impactor_Volume1-Impactor_Volume0);
                
                Pressure_Scaled=0.5*Impactor_Density*10^-2/Impacted_Volume*Impactor_Velocity^2;
                
                
%                 if dImpactor_Volume>0
%                     ext_F=F_Pressure*Pressure_Scaled;
%                     ext_F(nn*3)=ext_F(nn*3)*-0;
%                     tfixed=t;
%                 else
%                     ext_F=f_Volumetric*0;
%                     % ext_F(nn*3)=ext_F(nn*3)*(1-exp(-(t+10*dt)*100000));
%                     % (1-exp(-(t+10*dt)*100000))
%                     % ext_F=fT*10^5;
%                     % ext_F=fT*0;
%                     dImpactor_Volume=0;
%                     Impactor_Volume1=0;
%                     Impactor_Volume0=0;
%                     Pressure_Scaled=0;
%                 end


                %if t<100*dt
                %   ext_F=F_Pressure*10^10; 
                %else
                %   ext_F=f_Volumetric*0;
                %end
                
                ext_F=f_BC*0;
                
                ext_F(1:Mesh_Information.neq_Model)=F_K*ones(1356,1);%*sin(2*pi*10*t);
                
                dsa
                
                
                %ext_F(nn*3-2)=1000;
                %ext_F(nn*3-1)=1000;
                %ext_F(nn*3-0)=1000;

                 % ext_F(3:3:end)=0;
                 % ext_F(1173*3)=-0.0012;
                 
                 % [U1,Vel1,Acc1,strain_gp,sigma_gp,P_Eq]=Newmark_Beta_Plastic_Constant_K(Element_Information,Mesh_Information,Material_Information,FE_Information_Region,M,C,IKbar,a0,a1,a2,a3,a4,a5,ext_F,f_BC,sigma_gp,strain_gp,Dm,P_Eq,U0,Vel0,Acc0);
                 
                 
                 External_Force=ext_F.*f_BC;
                 Rbar=M*(a0*U0+a2*Vel0+a3*Acc0)+C*(a1*U0+a4*Vel0+a5*Acc0)+External_Force(1:Mesh_Information.neq_Model);  
                 U1 = IKbar*Rbar;

                % Solve for the corresponding velocities and accelerations
                 Acc1=a0*(U1-U0)-a2*Vel0-a3*Acc0;
                 Vel1=Vel0+a7*Acc1+a6*Acc0; 

                 Kinetic_Energy=0.5*Vel1'*M*Vel1;
                 Strain_Energy=0.5*U1'*K*U1;
                 Total_Energy=Kinetic_Energy+Strain_Energy;
                 
                 U0=U1;
                 Vel0=Vel1;
                 Acc0=Acc1;
                 
                 U_Store(:,count)=U0;
                 Vel_Store(:,count)=Vel0;
                 Acc_Store(:,count)=Acc0;
                 Time(count)=t;
                 Energy(1,count)=Strain_Energy;
                 Energy(2,count)=Kinetic_Energy;
                 Energy(3,count)=Total_Energy;
                 Force(1,count)=norm(ext_F);
                 Impactor(1,count)=Impactor_Volume1;
                 Impactor(2,count)=Pressure_Scaled;
                 % Force(2,count)=Impactor_Volume1;
                 
                 fprintf('Time=%1.6f Max plastic strain=%2.3f\n',t,max(max(P_Eq)))
                 % max(max(P_Eq))
                 count=count+1;

            end

            
            %% Explicit
            
            cD=sqrt(Material_Information.MoE/Material_Information.rho);
            IEN=Mesh_Information.IEN;
            Nodes=Mesh_Information.Nodes;
            DistanceO=10^9;
            for e=1:Mesh_Information.nel_Model
               for n1=1:Element_Information.nnp_Elemental(Mesh_Information.Type(1,e))-1
                   Node_n1=IEN(n1,e);
                   for n2=n1+1:Element_Information.nnp_Elemental(Mesh_Information.Type(1,e))
                       Node_n2=IEN(n2,e);
                       Distance=norm(Nodes(Node_n1,2:4)-Nodes(Node_n2,2:4));
                       if Distance<DistanceO
                           DistanceO=Distance;
                       end
                   end
               end
            end
            Le=DistanceO;
            
            dt=Le/cD/2;
            rt
            dt=10^-6;
            
            delta=0.5; alpha=0.25;
            % gama    % beta
            a0=1/alpha/dt^2; a1=delta/alpha/dt;  
            a2=1/alpha/dt; 
            a3=1/2/alpha-1; 
            a4=delta/alpha-1;
            a5=dt/2*(delta/alpha-2); 
            a6=dt*(1-delta); 
            a7=delta*dt;

            C=K*0;
            IM=inv(M);

            
            Nodes=Mesh_Information.Nodes;
            IEN=Mesh_Information.IEN;
            neq=Mesh_Information.neq_Model;

            U0=zeros(neq,1);Vel0h=U0;Vel0=U0;Acc0=U0;Accal=U0;
            U_Store=U0;
            Acc_Store=U0;
            Vel_Store=U0;
            nn=IEN(47,13);%IEN(71,13);
            Time=0;

            Impactor_Energy=0.5*4/3*pi*(50.E-03/2)^3*2703*10000^2;
            
            
            F_Pressure=f_Volumetric/1.7000e+11;
            Pressure=0.5*2700*2000^2;
            Impact_Time=2*(50.E-03/2)/2000;

            
            count=1;
            Total_Energy=0;
            Switch=0;
            for t=0:dt:1000*dt

                % if t<1*dt
                % if Total_Energy<Impactor_Energy && Switch==0                  
                if count<=Impact_Time/dt
                    ext_F=f_Volumetric*0.1;
                    ext_F=F_Pressure*Pressure;
                    ext_F(nn*3)=ext_F(nn*3)*0;
                    tfixed=t;
                    count
                else
                    % ext_F=f_Volumetric*0.15*0;
                    % ext_F=f_Volumetric*(exp(-(t-tfixed)*1000))*0.015;
                    % ext_F=F_Pressure*Pressure*(exp(-(t-tfixed)*1000));
                    % ext_F(nn*3)=ext_F(nn*3)*0;
                    % ext_F(nn*3)=ext_F(nn*3)*(1-exp(-(t-tfixed)*10));
                    % (1-exp(-(t+10*dt)*100000))
                    % ext_F=fT*10^5;
                    ext_F=ext_F*0;
                    Switch=1;
                end
                 % ext_F(3:3:end)=0;
                 % ext_F(1173*3)=-0.0012;

                 Acc1=IM*(ext_F.*f_BC-C*Vel0-K*U0);
                 Vel1h=Vel0h+dt*Acc1;
                 U1=U0+dt*Vel1h;

                 Vel0=(Vel1h+Vel0h)/2;
                 U0=U1;
                 Vel0h=Vel1h;
                 Acc0=Acc1;
                 
                 Kinetic_Energy=0.5*Vel0h'*M*Vel0h;
                 Strain_Energy=0.5*U1'*K*U1;
                 Total_Energy=Kinetic_Energy+Strain_Energy;
                 
                 U_Store(: ,count)=U0;
                 Vel_Store(:,count)=Vel0;
                 Acc_Store(:,count)=Acc0;
                 Time(count)=t;
                 Energy(1,count)=Kinetic_Energy;
                 Energy(2,count)=Strain_Energy;
                 Energy(3,count)=Total_Energy;
                 F_Store(count)=norm(ext_F);
                 
                 
                 

                 count=count+1;

            end            
            
    