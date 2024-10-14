function [SysM, SysN]=ApplyBC55(nnp,Npoip,Neq,Ibcu,SysM,SysN,Uvel,Vvel,Pres)

%Apply boundary condition for nodal U-velocity
IEQ1=1;
IEQ2=nnp;
for IEQ=IEQ1:IEQ2
    IEQU=IEQ;
    if(Ibcu(IEQU) ==0)
        for IR=1:Neq
            if(IR==IEQ)
                SysN(IR)=SysN(IR)-SysM(IR,IEQ)*Uvel(IEQU);
                SysM(IR,IEQ)=0;
            end
        end
        for IC=1:Neq
            SysM(IEQ,IC)=0;
        end
        
        SysM(IEQ,IEQ)=1;
        SysN(IEQ)=Uvel(IEQU);
    end
    


% Apply boundary condition for nodal V-velocity
IEQ1=nnp+1;
IEQ2=2*nnp;
for IEQ=IEQ1:IEQ2
    IEQV=IEQ-nnp;
    if(Ibcv(IEQV) ==0)
        for IR=1:Neq
            if(IR==IEQ)
                SysN(IR)=SysN(IR)-SysM(IR,IEQ)*Vvel(IEQV);
                SysM(IR,IEQ)=0;
            end
        end
        for IC=1:Neq
            SysM(IEQ,IC)=0;
        end
        SysM(IEQ,IEQ)=1;
        SysN(IEQ)=Vvel(IEQV);
    end
end

% Apply boundary condition for nodal W-velocity
IEQ1=2*nnp+1;
IEQ2=3*nnp;
for IEQ=IEQ1:IEQ2
    IEQW=IEQ-nnp*2;
    if(Ibcv(IEQW) ==0)
        for IR=1:Neq
            if(IR==IEQ)
                SysN(IR)=SysN(IR)-SysM(IR,IEQ)*Wvel(IEQW);
                SysM(IR,IEQ)=0;
            end
        end
        for IC=1:Neq
            SysM(IEQ,IC)=0;
        end
        SysM(IEQ,IEQ)=1;
        SysN(IEQ)=Wvel(IEQW);
    end
end



% Apply boundary condition for nodal P-pressure
IEQ1=2*nnp+1;
IEQ2=Neq;
for IEQ=IEQ1:IEQ2
    IEQP=IEQ-2*nnp;
    if(Ibcp(IEQP) ==0)
        for IR=1:Neq
            if(IR==IEQ)
                SysN(IR)=SysN(IR)-SysM(IR,IEQ)*Pres(IEQP);
                SysM(IR,IEQ)=0;
            end
        end
    for IC=1:Neq
        SysM(IEQ,IC)=0;
    end
    SysM(IEQ,IEQ)=1;
    SysN(IEQ)=Pres(IEQP);
    end
end

end
    
    