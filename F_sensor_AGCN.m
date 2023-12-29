function [isensors, Rinv, time]=F_sensor_AGCN(Xorg,r1,r2,p)
    %pre
    [U, Ur, Sr, dSdiag, ~, ~, ~]=F_pre_SVD_corrnoise(Xorg,r1,r2);
    
    tic;
    [n,r1]=size(U);
    Cpp=zeros(1,r1);
    Sr_sq=Sr*Sr;
    
    initial=true;
    tr_test=zeros(p,1);
    
    isensors=zeros(0,1);
    Cpp=zeros(0,r1);
    CCinv=inv(Cpp*Cpp');
    RC=Cpp;
    Rinv=zeros(0,0);
    R=zeros(0,0);
    A1   = zeros(1,r1);
    A2   = eye(r1,r1);
    for pp=1:p
        if pp<=r1
            CCinv=inv(Cpp*Cpp');
            tr_vec=zeros(1,n);
        %% searching
            for nn=1:n
                u_i=U(nn,:);
                s=zeros(1,pp-1);
                t=0;
                for l=1:(pp-1)
                    s(1,l)=Ur(nn,:)*Sr_sq*Ur(isensors(l,1),:)';
                end
                t=Ur(nn,:)*Sr_sq*Ur(nn,:)'+dSdiag(nn);   % same as rank of Ur
                
                nume=u_i*u_i'-u_i*Cpp'*CCinv*Cpp*u_i';
                %%AGwR
                aR=(CCinv*Cpp*u_i'*u_i*Cpp'*CCinv)/nume*R;
                bp=-(u_i*Cpp'*CCinv)/(nume)*s';
                btpt=-(CCinv*Cpp*u_i')/(nume)*s;
                dq=t/nume;
                tr_vec(1,nn)=trace(aR)+trace(bp)+trace(btpt)+trace(dq);
                %%AGmod
                
                %%%%det
                %nume=u_i*u_i'-u_i*Cpp'*CCinv*Cpp*u_i';
                %dnm=t-s*Rinv*s';
                %det_vec(1,nn)=dnm\nume;         
                %%%%let
            end
        
            for l=1:(pp-1)
                tr_vec(1,isensors(l,1))=10000000000;
            end
            [tr_test(pp,1),isensors(pp,1)]=min(tr_vec);   % argmaxdet
        
%%   Update Rinv&C after we get pp-th sensor  
            s=zeros(1,pp-1);
            u_i=U(isensors(pp,1),:);
        
            for l=1:(pp-1)
                s(1,l)=Ur(isensors(pp,1),:)*Sr_sq*Ur(isensors(l,1),:)';
            end
            t=Ur(isensors(pp,1),:)*Sr_sq*Ur(isensors(pp,1),:)'+dSdiag(isensors(pp,1));   % same as rank of Ur
            R=[R s';s t];
            dnm=t-s*Rinv*s';

            Cpp=[Cpp;u_i]; %#ok<AGROW>

            sR=s*Rinv;
            Rinv_new=zeros(pp,pp);
            Rinv_new(1:pp-1,1:pp-1)=Rinv;
            Rinv=Rinv_new+[sR';-1]*(dnm\[sR -1]);

        elseif pp>r1 %=========
               RC=Rinv*Cpp;
                tr_vec=zeros(1,n);
            %% searching
            for nn=1:n
                u_i=U(nn,:);
                s=zeros(1,pp-1);
                t=0;
                for l=1:(pp-1)
                    s(1,l)=Ur(nn,:)*Sr_sq*Ur(isensors(l,1),:)';
                end
                t=Ur(nn,:)*Sr_sq*Ur(nn,:)'+dSdiag(nn);   % same as rank of Ur
                CTRCinv = inv(Cpp'*RC);
                B=Cpp'*Rinv*s'-u_i';
                C=s*RC-u_i;
                a=1/(t-s*Rinv*s');
                tr_vec(1,nn)=-trace(a*CTRCinv*B*inv(1+a*C*CTRCinv*B)*C*CTRCinv);
            end
            
            for l=1:(pp-1)
                tr_vec(1,isensors(l,1))=1000000;
            end
            [tr_test(pp,1),isensors(pp,1)]=min(tr_vec);   % argmaxdet
            
    %%   Update Rinv&C after we get pp-th sensor  
            s=zeros(1,pp-1);
            t=0;
            u_i=U(isensors(pp,1),:);
            
            for l=1:(pp-1)
                s(1,l)=Ur(isensors(pp,1),:)*Sr_sq*Ur(isensors(l,1),:)';
            end
            t=Ur(isensors(pp,1),:)*Sr_sq*Ur(isensors(pp,1),:)'+dSdiag(isensors(pp,1));   % same as rank of Ur
            dnm=t-s*Rinv*s';
   
            Cpp=[Cpp;u_i];
            
            sR=s*Rinv;
            Rinv_new=zeros(pp,pp);
            Rinv_new(1:pp-1,1:pp-1)=Rinv;
            Rinv=Rinv_new+[sR';-1]*(dnm\[sR -1]);
        end
    end
    time=toc;
end