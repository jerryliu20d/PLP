pars=csvread('dataset');
set_num=200;
for i = 1:(size(pars,1)/4)

    par_t=pars(i,:);
    for noz = length(par_t):-1:1
	if par_t(noz)==0
             par_t=par_t(1:(noz-1));
        else
             break
        end
    end
    dataset=par_t(1);
    rng('shuffle');
    dirname=strcat('data',num2str(dataset));
    if ~isempty(dir(dirname))
        rmdir(dirname,'s')
    end
    mkdir(dirname);
    cd(dirname);
    nc=floor((length(par_t)-3)/3);
    m=par_t(3);
    tau=par_t(4:(3+nc));
    par=[par_t((4+nc):(4+2*nc));par_t((4+2*nc+1):(4+3*nc+1))]';
    for j = 1:set_num
        [z,Nj,C]=latent_simu_f(tau,m,par);
	if dataset == '9'
	   print(z)
	end
        filename1=['censortm' num2str(j) '.csv'];
        filename2=['eventtm' num2str(j) '.csv'];
        filename3=['Nt' num2str(j) '.csv'];
        csvwrite(filename1,C);
        dlmwrite(filename2,z,'precision',10)
        csvwrite(filename3,Nj)
    end
    cd '..'
end

function l=Lambda(x,tau,par)
l=0;
for i=1:length(tau)
    if x<=tau(i)
        if i>1
            l=l+((x/par(i,2))^par(i,1)-(tau(i-1)/par(i,2))^par(i,1));
        else
            l=l+(x/par(i,2))^par(i,1);
        end
        break
    else
        if i>1
            l=l+(tau(i)/par(i,2))^par(i,1)-(tau(i-1)/par(i,2))^par(i,1);
        else
            l=l+(tau(i)/par(i,2))^par(i,1);
        end
        if i==length(tau)
            l=l+(x/par(i+1,2))^par(i+1,1)-(tau(i)/par(i+1,2))^par(i+1,1);
        end
    end
end
end

function [z,Nj,C]=latent_simu_f(tau,m,par)% simulate date with K unique values of change-points
C=unifrnd(450,500,m,1);% censoring time
Nj=zeros(1,m);% # of events for each driver
Ft=@(x,c,tau2,par2) (Lambda(x,tau2,par2)/Lambda(c,tau2,par2));
avgnum2=40;
z=-1000*ones(m,avgnum2);
for j=1:m
    n=poissrnd(Lambda(C(j),tau,par));
    X0=zeros(1,n);
    i=1;
    while X0(n)==0
        u=unifrnd(0,1);
        a=0;
        b=C(j);
        for s=1:20
            temp=Ft((a+b)/2,C(j),tau,par);
            if temp<=u 
                binf=(a+b)/2; bsup=b;
            elseif temp>u 
                bsup=(a+b)/2;binf=a;
            end   
            a=binf;
            b=bsup;
        end
        if i==1
            X0(i)=(a+b)/2;
            i=i+1;
        else
            tol=1e-4;flag=0;
            for test=1:(i-1)
                if abs(X0(test)-(a+b)/2)<tol
                    flag=1;
                end
            end
            if flag==0
                X0(i)=(a+b)/2;
                i=i+1;
            end
        end
    end
    z(j,1:n)=sort(X0);% data for one driver
    Nj(j)=n;
end
z(z<0)=0;
end
