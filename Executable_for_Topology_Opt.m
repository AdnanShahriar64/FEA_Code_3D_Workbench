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

      Name='Block.inp';
      S_Names={'S1'};
      Reg_Name={'Sx','S1'};
      BC_Name={'S2','P1'};

      [Mesh_Information]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information);
      
      Mesh_Information.Boundary_Cond=[1 0 0;
                                      1 0 0
                                      1 0 0]; 


     %% Rotation
     
     Nodes=Mesh_Information.Nodes;

     Alpha=90*pi/18*0;
     Beta=-90*pi/18*0;
     Gamma=45*pi/18*0;
     
     Rz=[cos(Gamma) -sin(Gamma) 0;sin(Gamma) cos(Gamma) 0; 0 0 1];
     Ry=[cos(Beta) 0 sin(Beta);0 1 0;-sin(Beta) 0 cos(Beta)];
     Rx=[1 0 0; 0 cos(Alpha) -sin(Alpha);0 sin(Alpha) cos(Alpha)];
     R1=Rx*Ry*Rz;

     for i=1:length(Nodes)
         Nodes(i,2:4)=inv(R1)*[Nodes(i,2);Nodes(i,3);Nodes(i,4)];
     end

     Mesh_Information.Nodes=Nodes;
     

%% Material input      
      Material_Information.MoE=68*10^9;
      Material_Information.PoRat=0.3;
      Material_Information.rho=2703;
      Material_Information.alpe=10^-6;

      [Dm]=Material_Array_Generator(Material_Information);     


%% System preprocessor
    
    [Mesh_Information]=Geometry_Model_Parameters_Preprocessor_Solid(Mesh_Information,Element_Information);

    [FE_Information_Region]=Region_Integration_3D_Solid(Element_Information,Mesh_Information);

    [FE_Information_Region]=Region_Integration_3D_Solid_Extended(Element_Information,Mesh_Information,FE_Information_Region,Material_Information);

    [FE_Information_Surface]=Surface_Integration_3D_Pressure_Type_Solid(Element_Information,Mesh_Information);

      
%% Stiffness matrix



    [K]=Elemental_Stiffness_Matrix_Generator_GP_Wise(Element_Information,Mesh_Information,Material_Information,FE_Information_Region,Dm);
    
    % [K]=Elemental_Stiffness_Matrix_Generator_GP_Wise(B,Dm,dJB,W,Type,1:nel,nel,BEx,LM,neq);
    

    [Bc]=Boundary_Condition_Set(Mesh_Information);
    f_BC=Bc*-1+1;

    [K]=Boundary_Set(Bc,Mesh_Information.neq_Model,K);

    Fs=zeros(Mesh_Information.neq_Model,1);
    Fs(Mesh_Information.BC_Nodes(1,2)*3)=-100;
      
    d=inv(K)*Fs;

