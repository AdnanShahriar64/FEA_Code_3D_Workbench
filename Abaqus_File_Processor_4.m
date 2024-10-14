function[Mesh_Information,Nodes,IEN,Type,S_Nodes,S_Type,Reg_Nodes,Reg_Elems,BC_Nodes]=Abaqus_File_Processor_4(Name,S_Names,Reg_Name,BC_Name,Element_Information)
%%)

% Name='FEM_Matlab_25.inp';
% S_Names={'Sx','S1'};
% Reg_Name={'Sx','S1'};
% BC_Name={'Sx','S1'};

fprintf('\n\n------------------------------------------------------\n')
fprintf('Converting Abaqus to Matlab\n')
fprintf('------------------------------------------------------\n')
fprintf('Developer: Adnan Shahriar\n           Contact: adnan.shahriar@utsa.edu\n           Website: adnansh.com\n           Department of Mechanical Engineering, UTSA, Texas\n')
fprintf('------------------------------------------------------\n')

% fprintf(2,'Reading and Importing Abaqus Input File...\n')
fprintf('Reading and Importing Abaqus Input File...\n')

% Element assignment
E_Type={'C3D4','C3D4R','C3D10','C3D10R','C3D6','C3D6R','C3D15','C3D15R','C3D8','C3D8R','C3D20','C3D20R'
           2      2       3        3       4      4       5        5      6       6       7        7};



% E_Type=[1  2 3  4  5 5  6  8;   % Element Type
%         3 10 6 15 18 8 20 27;   % Node Number
%         1  1 6 21 21 8 27 20];  % Gauss Points
% 
% % 8 and 20 Node Brick Element Surface Data    
% Surf6=[2     3     6     7  % g=1
%        1     4     5     8  % g=-1
%        3     4     7     8  % h=1
%        1     2     5     6  % h=-1
%        5     6     7     8  % r=1
%        1     2     3     4];% r=-1    
%    
% % 6 and 15 Node Wedge Element Data   
% Surf3=[9 9 9;9 9 9; 9 9 9; % x, y surface are not incorporated
%        4 5 6;       % r=1
%        1 2 3];      % r=-1
    

   
% Import    
    F_Name = Importfile(Name);

% Node Import
    fprintf('Importing Nodes...')
    
    for i=1:length(F_Name)
        if F_Name(i,1)=='*Node'
            N_s=i+1;
        elseif F_Name(i,1)=='*Element'
            N_e=i-1;
            break
        end
    end

    Nodes=zeros(N_e-N_s,4);
    count=1;
    for i=N_s:N_e
        Nodes(count,1)=str2num(F_Name(i,1));
        Nodes(count,2)=str2num(F_Name(i,2));
        Nodes(count,3)=str2num(F_Name(i,3));
        Nodes(count,4)=str2num(F_Name(i,4));
        count=count+1;
    end
    fprintf('Done\n')

%% Element Import
fprintf('Importing Elements...')

    E_d=0;          % E_d store all Element array data
                    % [ Element1 start      Element2 start ...
                    %   Element1 end        Element2 end   ...
                    %   Element1 Type       Element2 Type  ...]
    count=1;
    for i=N_e:length(F_Name)
        
        % Get E_d
        if F_Name(i,1)=='*Element'
            E_d(1,count)=i+1;
            if count>1
                E_d(2,count-1)=i-1; 
            end

            for eT=1:size(E_Type,2)
                ET=sprintf('type=%s',E_Type{1,eT});
                if F_Name(i,2)==ET
                    E_d(3,count)=E_Type{2,eT};
                end

            end

%             if F_Name(i,2)=='type=C3D8' || F_Name(i,2)=='type=C3D8R'
%                 E_d(3,count)=6;
%             elseif F_Name(i,2)=='type=C3D6' || F_Name(i,2)=='type=C3D6R'
%                 E_d(3,count)=3;
%             elseif F_Name(i,2)=='type=C3D20' || F_Name(i,2)=='type=C3D20R'
%                 E_d(3,count)=7;
%             elseif F_Name(i,2)=='type=C3D15' || F_Name(i,2)=='type=C3D15R'
%                 E_d(3,count)=4;
%             end

            count=count+1;
        elseif F_Name(i,1)=='*Nset' || F_Name(i,1)=='*End Part'
            E_d(2,count-1)=i-1;
            break
        end
    end
    
%%

    % Generate 
    IEN=0;
    S=size(E_d);
    
    for eType=1:S(2)    % Go through each element types
       
        i=E_d(1,eType);
        nne=Element_Information.nnp_Elemental(E_d(3,eType));
        
        while i<=E_d(2,eType)
            
            enum=str2num(F_Name(i,1));
            
            n=1;
            j=2;
            while n<=nne
                IEN(n,enum)=str2num(F_Name(i,j));
                j=j+1;
                n=n+1;
                if j>16 && nne>16 %E_d(3,eType)>6 
                    j=1;
                    i=i+1;
                end
            end
            
            Type(1,enum)=E_d(3,eType);
            Type(2,enum)=Element_Information.nnp_Elemental(E_d(3,eType));
            Type(3,enum)=Element_Information.ngp_Elemental(E_d(3,eType));
            Type(4,enum)=1;
            Type(5,enum)=0;
            Type(6,enum)=0;
            
            i=i+1;
        end
        
        
    end
    fprintf('Done\n')


    
%% Surface Import   
    % Get Surface Nodes
    %fprintf('\n--------------------------------------\n')
    fprintf('Importing surface sets...\n')
    %fprintf('--------------------------------------\n')

    S_Nodes=0;
    for SN=1:length(S_Names)
       
        % Check Start and End
        Sd=sprintf('nset=%s',S_Names{SN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd
                St=i+1;
            elseif F_Name(i,1)=='*Elset' && St>0
                En=i-1;
                break
            end
        end
        
        if St==0 && En==0
            fprintf(2,'   Warning! Nodes of "%s" surface set can not be detected\n',S_Names{SN});
        else

                % Assign Nodes
                if F_Name(St-1,4)=='generate'
                    st=str2num(F_Name(St,1));
                    en=str2num(F_Name(St,2));
                    ns=[st:str2num(F_Name(St,3)):en];
                    S_Nodes(1:length(ns),SN)=ns;     
                else
                    count=1;
                    for i=St:En
                        for j=1:16
                            if F_Name(i,j)~=""
                                S_Nodes(count,SN)=str2num(F_Name(i,j));
                                count=count+1;
                            end
                        end
                    end
                end  
                
                fprintf('   Success! Nodes of "%s" surface set has been detected\n',S_Names{SN});
        end
    end
        

        % Get Surface Elements
        count=1;
        S_Type_B=0;
        for SN=1:length(S_Names)

            Sd=sprintf('elset=%s',S_Names{SN});

            % Check Start and End
            St=0;
            En=0;
            for i=1:length(F_Name)
                if F_Name(i,2)==Sd
                    St=i+1;
                elseif (F_Name(i,1)=='*Elset' || F_Name(i,1)=='*Nset' || F_Name(i,1)=='*End Assembly' || F_Name(i,1)=='*Surface') && St>0
                    En=i-1;
                    break
                end
            end

            % Breaking loop if no surface is detected
            if St==0 && En==0
                fprintf(2,'   Warning! Elements of "%s" surface set can not be detected\n',S_Names{SN});
            else

                % Assign Elements
                if F_Name(St-1,4)=='generate'
                    st=str2num(F_Name(St,1));
                    en=str2num(F_Name(St,2));
                    ns=[st:str2num(F_Name(St,3)):en];
                    S_Type_B(count:count+length(ns)-1,2)=ns; 
                    S_Type_B(count:count+length(ns)-1,1)=ns*0+SN;
                    count=count+length(ns);
                else

                    for i=St:En
                        for j=1:16
                            if F_Name(i,j)~=""
                                S_Type_B(count,2)=str2num(F_Name(i,j));
                                S_Type_B(count,1)=SN;
                                count=count+1;
                            end
                        end
                    end

                end 
                fprintf('   Success! Elements of "%s" surface set has been detected\n',S_Names{SN});
            end
        end
    
        
    %% Hardest Part! Surface type extraction!
    S=size(S_Type_B);
    
    c=1;
    for i=1:S(1)               % For each Element

        SN=S_Type_B(i,1);         % Extract Surface Number
        
            % Breaking loop if no surface is detected
            if S_Type_B==0
                break
            end

            e=S_Type_B(i,2);      % Extract the Element

            % Later this Stupid algorithm should be replaced!
            count=1;
            L=0;
            nn=Element_Information.nnp_Elemental(Type(1,e));
            for k=1:nn   % Extract the array no by checking each Nodes from List

                for l=1:length(S_Nodes)         % Run through all nodes and check!
                    if S_Nodes(l,SN)==IEN(k,e)
                        L(count,1)=k;
                        count=count+1;
                        break
                    end
                end

            end

            ss=size(Element_Information.Surf_ngp_Elemental,1);
                cs=1;
                for l=1:ss

                    SN_Limit=Element_Information.Surf_nnp_Elemental(l,Type(1,e));

                    if SN_Limit>0
                        SNS=Element_Information.Surf_Nodes_Elemental(1:SN_Limit,l,Type(1,e));
                        if length(intersect(SNS,L))==SN_Limit
                           S_Type_BU(c,3)=l;
                           S_Type_BU(c,1:2)=S_Type_B(i,1:2);
      
                           cs=cs+1;
                           c=c+1;
                        end
                    end
                    
                end

                if cs>3
                    fprintf(2,'Something wrong with Surface type extraction\n')
                end


    end

    if S_Type_B~=0
        S_Type_B=S_Type_BU;
    end

    S_Type=S_Type_B';


%% Boundary Extract
    BC_Nodes=0;
    %fprintf('\n--------------------------------------\n')
    fprintf('Importing Boundary sets...\n')
    %fprintf('--------------------------------------\n')
    for BN=1:length(BC_Name)
        
        Sd=sprintf('nset=%s',BC_Name{BN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd
                St=i+1;
            elseif (F_Name(i,1)=='*Elset' || F_Name(i,1)=='*Nset') && St>0
                En=i-1;
                break
            end
        end
        
        if St==0 && En==0
            fprintf(2,'   Warning! "%s" boundary set can not be detected\n',BC_Name{BN});
        else
            
        
        
        if F_Name(St-1,4)=='generate' 
            st=str2num(F_Name(St,1));
            en=str2num(F_Name(St,2));
            ns=[st:str2num(F_Name(St,3)):en];
            BC_Nodes(1:length(ns),BN)=ns;     
        else        
            count=1;
            for i=St:En
                for j=1:16
                    if F_Name(i,j)~=""
                        BC_Nodes(count,BN)=str2num(F_Name(i,j));
                        count=count+1;
                    end
                end
            end
        end  
        
        
        fprintf('   Success! "%s" boundary set has been detected\n',BC_Name{BN});
        
        end
        
    end
    
 Type(end-2,:)=1;   
 Type(end,:)=0;   
    
%% Material Extract

% for MN=1:length(Mat_Names)
%     
%     Sd=sprintf('name=%s',Mat_Names{MN});
%     
%     % Get element data
%     for i=1:length(F_Name)
% 
%         if F_Name(i,2)==Sd
%             rho(MN,1)=str2num(F_Name(i+2,1));
%             MoE(MN,1)=str2num(F_Name(i+4,1));
%             PoRat(MN,1)=str2num(F_Name(i+4,2));
%             break
%         end
%         
%     end
%     
%     Sd2=sprintf('material=%s',Mat_Names{MN});
%     % Update the Type 
%     for i=1:length(F_Name)
% 
%         % Look for the Section part
%         if F_Name(i,3)==Sd2
%             
%             Elem_search=F_Name(i,2);
%             
%             % Look with the element name from Section part
%             for j=1:length(F_Name)
%                 
%                 % Got the part! Now Update Type
%                 if F_Name(j,1)=='*Elset' && F_Name(j,2)==Elem_search
%                     
%                     if F_Name(j,4)=='generate' || F_Name(j,3)=='generate'
%                        
%                         EList=[str2num(F_Name(j+1,1)):str2num(F_Name(j+1,3)):str2num(F_Name(j+1,2))];
%                         
%                         Type(4,EList)=MN;
%                         
%                     else
%                         
%                         R_Up=1;
%                         while R_Up~=0
%                             for k=1:16
%                                 if F_Name(j+1,k)==""
%                                     break
%                                 end
%                                 EList=str2num(F_Name(j+1,k));
%                                 Type(4,EList)=MN;
%                                 
%                             end
% 
%                             j=j+1;
%                             
%                             if F_Name(j+1,1)=="*Nset"
%                                 R_Up=0;
%                             end
%                             if F_Name(j+1,1)=="*Elset"
%                                 R_Up=0;
%                             end
%                             
%                         end
%                     end
%                     
%                     break
%                 end
%   
%             end
%                break
%         end
%         
%     end
%           
% end


%% Region Extract

%% Region Element Extract
    %fprintf('\n--------------------------------------\n')
    fprintf('Importing Region sets...\n')
    %fprintf('--------------------------------------\n')
    
    Reg_Elems=0;
    for RN=1:length(Reg_Name)
        
        Sd=sprintf('elset=%s',Reg_Name{RN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd
                St=i+1;
            elseif sum(str2num(F_Name(i,1)))==0 && St>0
                En=i-1;
                break
            end
        end
        
        
       if St==0 && En==0  
       else
        
            if F_Name(St-1,4)=='generate' 
                st=str2num(F_Name(St,1));
                en=str2num(F_Name(St,2));
                ns=[st:str2num(F_Name(St,3)):en];
                Reg_Elems(1:length(ns),RN)=ns;     
            else        
                count=1;
                for i=St:En
                    for j=1:16
                        if F_Name(i,j)~=""
                            Reg_Elems(count,RN)=str2num(F_Name(i,j));
                            count=count+1;
                        end
                    end
                end
            end  
       end
        
    end
    

        

%% Region Node Extract
    c=0;
    Siz=size(Reg_Elems);
    for RN=1:length(Reg_Name)
        
        Sd=sprintf('nset=%s',Reg_Name{RN});
        Sd2=sprintf('elset=%s',Reg_Name{RN});
        St=0;
        En=0;
        for i=1:length(F_Name)
            if F_Name(i,2)==Sd %|| F_Name(i,2)==Sd2
                St=i+1;
            elseif (F_Name(i,1)=='*Elset' || F_Name(i,1)=='*Nset') && St>0
                En=i-1;
                break
            end
        end
        
        
        if St==0 && En==0
            fprintf(2,'   Warning! "%s" Region set can not be detected\n',Reg_Name{RN});
            c=c+1;
        else
 
            if St>0

                if F_Name(St-1,4)=='generate' 
                    st=str2num(F_Name(St,1));
                    en=str2num(F_Name(St,2));
                    ns=[st:str2num(F_Name(St,3)):en];
                    Reg_Nodes(1:length(ns),RN)=ns;     
                else        
                    count=1;
                    for i=St:En
                        for j=1:16
                            if F_Name(i,j)~=""
                                Reg_Nodes(count,RN)=str2num(F_Name(i,j));
                                count=count+1;
                            end
                        end
                    end
                end   

            elseif St==0

                Reg_Nodes=0;
                fprintf(2,'\n     WARNING! No Region Nodes were found! Probably, region has been created using Elements. Generate Externally if you want so much')
                fprintf('\n')            
            end
            
            
            fprintf('   Success! "%s" Region set has been detected\n',Reg_Name{RN});
        end
        
    end
    
    if c==length(Reg_Name)
       Reg_Nodes=0;
    end
    
    
    

    Mesh_Information.Nodes=Nodes;
    Mesh_Information.IEN=IEN;
    Mesh_Information.Type=Type;
    Mesh_Information.S_Nodes=S_Nodes;
    Mesh_Information.S_Type=S_Type;
    Mesh_Information.Reg_Nodes=Reg_Nodes;
    Mesh_Information.Reg_Elements=Reg_Elems;
    Mesh_Information.BC_Nodes=BC_Nodes;



    % fprintf('\n..Done\n')
    fprintf('------------------------------------------------------\n')

end    
        


%% Importing Function
function Ab2Mat = Importfile(filename)

startRow=1;
fid=fopen(filename);
g = textscan(fid,'%s','delimiter','\n');
fclose(fid);
endRow=length(g{1});


%IMPORTFILE Import numeric data from a text file as a matrix.
%   AB2MAT = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   AB2MAT = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   Ab2Mat = importfile('Ab2Mat.txt', 1, 4493);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/05/30 17:35:22

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: text (%s)
%	column8: text (%s)
%   column9: text (%s)
%	column10: text (%s)
%   column11: text (%s)
%	column12: text (%s)
%   column13: text (%s)
%	column14: text (%s)
%   column15: text (%s)
%	column16: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
Ab2Mat = [dataArray{1:end-1}];

end
