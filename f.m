function mop = f( testname, dimension )
%Get test multi-objective problems from a given name. 
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The test problem will be encapsulated in a structure,
%   which can be obtained by function get_structure('testmop'). 
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.

    mop=get_structure('testmop');   
    switch lower(testname)
        case 'f1'
            mop=f1(mop, dimension);
        case 'f2'
            mop=f2(mop, dimension);
        case 'f3'
            mop=f3(mop, dimension);
        case 'f4'
            mop=f4(mop, dimension);
        case 'f5'
            mop=f5(mop, dimension);     
        case 'f6'
            mop=f6(mop, dimension); 
        case 'f7'
            mop=f7(mop, dimension); 
        case 'f8'
            mop=f8(mop, dimension); 
        case 'f9'
            mop=f9(mop, dimension);                              
        otherwise 
            error('Undefined test problem name');                
    end 
end

%%
function p = f1(p,dim)
 p.name = 'F1';
 p.pd = dim;% default dim=30
 p.od = 2;
 p.domain = [0,1;-1*ones(dim-1,1),ones(dim-1,1)];
 p.func = @evaluate;

 function y = evaluate(x)
    n = size(x,1);
    io = 3:2:n;
    ie = 2:2:n;
    t = x-0.9*sin(pi*(1:1:n)'/n); 
    g1 = 2.0/length(io)*(sum(abs(t(io)).^0.7));
    g2 = 2.0/length(ie)*(sum(abs(t(ie)).^0.7));
    y = zeros(2,1);
    y(1) = g1+(1-cos(pi*x(1)/2));
    y(2) = g2+(10-10*sin(pi*x(1)/2));
 end
end
%%
function p = f2(p,dim)
 p.name = 'F2';
 p.pd = dim;% default dim=30
 p.od = 2;
 p.domain = [0,1;-1*ones(dim-1,1),ones(dim-1,1)];
 p.func = @evaluate;

 function y = evaluate(x)
    n = size(x,1);
    io = 3:2:n;
    ie = 2:2:n;
    t = x-0.9*sin(pi*(1:1:n)'/n);
    g1 = 2.0/length(io)*(sum(abs(t(io)+sin(1*pi*(t(io)))/(1*pi)).^1));
    g2 = 2.0/length(ie)*(sum(abs(t(ie)+sin(1*pi*(t(ie)))/(1*pi)).^1)); 
    y = zeros(2,1);
    y(1) = g1+x(1);
    if x(1)<=0.05
       y(2) = g2+(1-19*x(1));
    else
       y(2) = g2+(1/19-x(1)/19);
    end
 end
end
%%
function p = f3(p,dim)
 p.name = 'F3';
 p.pd = dim;% default dim=30
 p.od = 2;
 p.domain = [0,1;-1*ones(dim-1,1),ones(dim-1,1)];
 p.func = @evaluate;

 function y = evaluate(x)
    n = size(x,1);
    io = 3:2:n;
    ie = 2:2:n;
    t = x-0.9*sin(pi*(1:1:n)'/n);
    g1 = 2.0/length(io)*(sum(1-exp(-1*abs(t(io)).^1)));
    g2 = 2.0/length(ie)*(sum(1-exp(-1*abs(t(ie)).^1)));
    y = zeros(2,1);
    y(1) = g1+x(1);
    if x(1)<=0.5
       y(2) = g2+1-0.5*(2*x(1)).^4;
    else
       y(2) = g2+0.5*(2*(1-x(1))).^4;
    end
 end
end
%%
function p = f4(p,dim)
 p.name = 'F4';
 p.pd = dim;% default dim=30
 p.od = 2;
 p.domain = [0,1;-1*ones(dim-1,1),ones(dim-1,1)];
 p.func = @evaluate;
 
 function y = evaluate(x)    
    n = size(x,1);
    io = 3:2:n;
    ie = 2:2:n;
    t1 = x(io)-1.8*abs(x(1)-0.55)*cos(5.0*pi*(x(1))+(io)'*pi/n);
    t2 = x(ie)-1.8*abs(x(1)-0.55)*sin(5.0*pi*(x(1))+(ie)'*pi/n);
    g1 = 2.0/length(io)*sum(t1.^2.0);
    g2 = 2.0/length(ie)*sum(t2.^2.0);
    y = zeros(2,1);
    y(1) = x(1)+g1;
    y(2) = (1.0-sqrt(x(1)))+g2;
 end
end

%%
function p = f5(p,dim)
 p.name = 'F5';
 p.pd = dim;% default dim=30
 p.od = 2;
 p.domain = [0,1;-1*ones(dim-1,1),ones(dim-1,1)];
 p.func = @evaluate;
 
 function y = evaluate(x)
    n = size(x,1);
    io = 3:2:n;
    ie = 2:2:n;
    t1 = x(io)-(0.75*(abs(x(1)-0.55).^2)*cos(16.0*pi*x(1)+4.0*(io)'*pi/n)+1.4*abs(x(1,:)-0.55)).*cos(4.0*pi*x(1)+(io)'*pi/n);
    t2 = x(ie)-(0.75*(abs(x(1)-0.55).^2)*cos(16.0*pi*x(1)+4.0*(ie)'*pi/n)+1.4*abs(x(1,:)-0.55)).*sin(4.0*pi*x(1)+(ie)'*pi/n);
    g1 = 2.0/length(io)*sum(t1.^2.0);
    g2 = 2.0/length(ie)*sum(t2.^2.0);   
    y = zeros(2,1);
    y(1) = x(1)+g1;
    y(2) = (1.0-sqrt(x(1)))+g2;
 end
end
%%
function p = f6(p,dim)
 p.name = 'F6';
 p.od = 2;
 p.pd = dim;% default dim=30
 p.domain = [0,1;-1*ones(dim-1,1),ones(dim-1,1)];
 p.func = @evaluate;
   
 function y = evaluate(x)
    n = size(x,1);
    io = 3:2:n;
    ie = 2:2:n;
    t = x-0.9*sin(2.0*pi*x(1)+pi*(1:1:n)'/n);
    t1 = t(io);
    t2 = t(ie);
    g1 = 0.02/length(io)*(length(io)+100*sum(t1.*t1)- sum(cos(10*pi*t1)));%0.25р╡©ирт
    g2 = 0.02/length(ie)*(length(ie)+100*sum(t2.*t2)- sum(cos(10*pi*t2)));
    y = zeros(2,1);
    y(1) = g1+x(1);
    y(2) = g2+(1.0-sqrt(x(1)));   
 end
end

%%
function p = f7(p,dim)
 p.name = 'F7';
 p.od = 2;
 p.pd = dim;% default dim=30
 p.domain = [0,1;-1*ones(dim-1,1),ones(dim-1,1)];
 p.func = @evaluate;
   
 function y = evaluate(x)
    y = zeros(2,1);
    n = size(x,1);
    io = 3:2:n;
    ie = 2:2:n;
    t = x-0.9*sin(2.0*pi*x(1)+pi*(1:1:n)'/n);
    t1 = t(io);
    t2 = t(ie);
    g1 = 0.05/length(io)*(100*sum(t1.*t1)-1*prod(cos(20*pi*t1./sqrt(io)'))+1);
    g2 = 0.05/length(ie)*(100*sum(t2.*t2)-1*prod(cos(20*pi*t2./sqrt(ie)'))+1);
    y(1) = g1+x(1);
    y(2) = g2+(1.0-sqrt(x(1)));   
 end
end

function p = f8(p,pdim)
 p.name ='F8';
 p.od = 3;
 p.pd = pdim;% default dim=10
 p.domain = [-1.0*ones(pdim,1) 1.0*ones(pdim,1)]; p.domain(1:2,1) = 0.0; p.domain(1:2,2) = 1.0;
 p.func = @evaluate;
 function y = evaluate(x)
    odim = 3;
    y = zeros(odim,1);
    xm = x(odim:pdim,1);
    d_vars = 0.9*ones(pdim-odim+1,1);
    for d_var = 1:odim-1
     d_vars = d_vars.*sin(2*pi*x(d_var)+pi*(odim:1:pdim)'/pdim);
    end
    t= xm-d_vars;
    g= sum(abs(t).^2);
    y(1) = (1-prod(x(1:odim-1,1)))+g;
    y(odim) = x(1)+g;
    for i = 2:odim-1
        y(i) = (1-prod(x(1:odim-i,1))*(1-x(odim-i+1)))+g;
    end               
 end
end
%%
function p = f9(p,pdim)
 p.name = 'F9';
 p.od = 3;
 p.pd = pdim;% default dim=10
 p.domain = [-1.0*ones(pdim,1) 1.0*ones(pdim,1)]; p.domain(1:2,1) = 0.0; p.domain(1:2,2) = 1.0;
 p.func = @evaluate;
 function y = evaluate(x)
     odim=3;
     y=zeros(odim,1);
     xm=x(odim:pdim,1);
     d_vars=0.9*ones(pdim-odim+1,1);
     for d_var=1:odim-1
      d_vars = d_vars.*sin(2*pi*x(d_var)+pi*(odim:1:pdim)'/pdim);
     end
     t=xm-d_vars;
     g=0.01*(sum(1+100*(t).^2-cos(4*pi*(t))));
     y(1)=(1-prod(cos(pi*x(1:odim-1,1)./2)))+g;
     y(odim)=(1-sin(pi*x(1)/2))+g;
     for i=2:odim-1
         y(i)=(1-prod(cos(pi*x(1:odim-i,1)./2))*sin(pi*x(odim-i+1)./2))+g;
     end  
 end
end

