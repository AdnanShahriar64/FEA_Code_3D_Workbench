
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

      Name='FEM_Matlab_25.inp';
      S_Names={'Sx','S1'};
      Reg_Name={'Sx','S1'};
      BC_Name={'BCx','BCy','BCz'};

      [Mesh_Information]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information);
      
     % [Nodes,IEN,Type,S_Nodes,S_Type,Reg_Nodes,~,BC_Nodes]=Abaqus_File_Processor('Beam_Hole_Matlab_Small.inp',{'Sy','BCx'},{'Reg'},{'BCx','BCy','BCz'},{'Material-1'});
     
     % [Nodes,IEN,Type,S_Nodes,S_Type,Reg_Nodes,~,BC_Nodes]=Abaqus_File_Processor('FEM_Matlab_25.inp',{'Sx','S1'},{'Reg'},{'BCx','BCy','BCz'},{'Material-1'});

     % [Nodes,IEN,Type,S_Nodes,S_Type,Reg_Nodes,~,BC_Nodes]=Abaqus_File_Processor('FEM_Matlab_25_Linear.inp',{'Sx','Sy','S1'},{'Reg'},{'BC'},{'Material-1'});
     
     % [Nodes,IEN,Type,S_Nodes,S_Type,Reg_Nodes,~,BC_Nodes]=Abaqus_File_Processor('Pipe_Hole.inp',{'S1','S2'},{'Reg'},{'BC2','BC1','BC3'},{'Material-1'});
     
     % [Nodes,IEN,Type,S_Nodes,S_Type,Reg_Nodes,~,BC_Nodes]=Abaqus_File_Processor('Void_Matlab.inp',{'S1','S2'},{'Reg'},{'BCx','BCy','BCz'},{'Material-1'});

   

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


%% Update model element    
      % [Mesh_Information]=Element_Converter_20_to_27(Mesh_Information,Element_Information);
      Mesh_Information.Boundary_Cond=[1 0 0;
                                      0 1 0
                                      0 0 1];      

%% Update extendened element      

% 1st element
        EN=13;

        Type=Mesh_Information.Type;  
        IEN=Mesh_Information.IEN;
        nnp_el=Element_Information.nnp_Elemental(Type(1,EN));
        Nodes=Mesh_Information.Nodes;

        Type(1,EN)=11;
        Type(5,EN)=nnp_el;
        Type(6,EN)=2;
        
        LN=max(max(IEN));
        
        for i=1:Type(5,EN)
          % Nodes(LN+i,1:4)=[LN+i Nodes(IEN(i),2) Nodes(IEN(i),3) Nodes(IEN(i),4)]; 
          IEN(nnp_el+i,EN)=LN+i;
        end   

        Mesh_Information.Type=Type;
        Mesh_Information.IEN=IEN;

% 2nd element
        EN=14;

        nel=size(Mesh_Information.Type,2);
        Type=Mesh_Information.Type;  
        IEN=Mesh_Information.IEN;
        nnp_el=Element_Information.nnp_Elemental(Type(1,EN));
        Nodes=Mesh_Information.Nodes;

        Type(1,EN)=11;
        Type(5,EN)=nnp_el;
        Type(6,EN)=2;
        
        LN=max(max(IEN));
        
        count=1;
        for i=1:Type(5,EN)
          % Nodes(LN+i,1:4)=[LN+i Nodes(IEN(i),2) Nodes(IEN(i),3) Nodes(IEN(i),4)]; 
          
          nn=IEN(i,EN);
          
          Got_one=0;
          for e=setdiff(1:nel,EN)
              for n=1:nnp_el
                 nn_roam=IEN(n,e);
                 if nn==nn_roam && IEN(end,e)>0
                    IEN(nnp_el+i,EN)=IEN(nnp_el+n,e);
                    Got_one=1;
                    break
                 end 
              end
          end
             
          if Got_one==0
              IEN(nnp_el+i,EN)=LN+count;
              count=count+1;
          end
          
        end   

        Mesh_Information.Type=Type;
        Mesh_Information.IEN=IEN;



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




dcs

    
%% Solution    

    d=inv(K)*(fS.*f_BC);
    
    