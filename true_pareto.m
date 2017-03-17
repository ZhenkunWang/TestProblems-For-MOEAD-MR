% pareto.m
% 
% The Matlab source codes to generate the PF and the PS of the test
% instances in the paper "On the use of two reference points in 
% decomposition based multiobjective evolutionary algorithms".
% 
% Usage: [pf, ps] = pareto(problem_name, no_of_points, variable_dim)
% 
% Please refer to the report for more information.
% 
% History:
%   v1 Feb.07 2017

function [pf, ps] = true_pareto(name, no, dim)

    if nargin<3, dim = 3; end
    if nargin<2, no  = 500; end
    switch name
        case {'F1'}
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = 0.9*sin(pi*repmat((2:dim)',[1,no])/dim);           
            pf          = zeros(2,no);
            pf(1,:)     = 1-cos(pi*ps(1,:)/2);
            pf(2,:)     = 10-10*sin(pi*ps(1,:)/2);
        case {'F2'}
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = 0.9*sin(pi*repmat((2:dim)',[1,no])/dim);               
            pf(1,1:round(0.5*no)) = linspace(0,0.05,round(0.5*no));
            pf(1,(round(0.5*no)+1):no) = linspace(0.05,1,no-round(0.5*no));
            pf(2,1:round(0.5*no)) = 1-19*pf(1,1:round(0.5*no));
            pf(2,(round(0.5*no)+1):no) = 1/19-pf(1,(round(0.5*no)+1):no)/19;
        case {'F3'}
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = 0.9*sin(pi*repmat((2:dim)',[1,no])/dim);        
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);   
            pf(2,1:round(0.5*no)) = 1-0.5*(2*ps(1,1:round(0.5*no))).^4;
            pf(2,(round(0.5*no)+1):no) = 0.5*(2*(1-ps(1,(round(0.5*no)+1):no))).^4;        
        case {'F4'}
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            io = 3:2:dim;
            ie = 2:2:dim;
            ps(io,:) = 1.8*abs(repmat(ps(1,:),[length(io),1])-0.55).*cos(5*pi*repmat(ps(1,:),[length(io),1])+pi*repmat((io)',[1,no])/dim);
            ps(ie,:) = 1.8*abs(repmat(ps(1,:),[length(ie),1])-0.55).*sin(5*pi*repmat(ps(1,:),[length(ie),1])+pi*repmat((ie)',[1,no])/dim);
        case {'F5'}
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            io = 3:2:dim;
            ie = 2:2:dim;
            ps(io,:) = (0.75*abs(repmat(ps(1,:),[length(io),1])-0.55).^2.*cos(16*pi*repmat(ps(1,:),[length(io),1])+4*pi*repmat((io)',[1,no])/dim)+1.4*abs(repmat(ps(1,:),[length(io),1])-0.55)).*cos(4*pi*repmat(ps(1,:),[length(io),1])+pi*repmat((io)',[1,no])/dim);
            ps(ie,:) = (0.75*abs(repmat(ps(1,:),[length(ie),1])-0.55).^2.*cos(16*pi*repmat(ps(1,:),[length(ie),1])+4*pi*repmat((ie)',[1,no])/dim)+1.4*abs(repmat(ps(1,:),[length(ie),1])-0.55)).*sin(4*pi*repmat(ps(1,:),[length(ie),1])+pi*repmat((ie)',[1,no])/dim);            
         case {'F6','F7'}
            pf          = zeros(2,no);
            pf(1,:)     = linspace(0,1,no);
            pf(2,:)     = 1-sqrt(pf(1,:));
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = 0.9*sin(2*pi*repmat(ps(1,:),[dim-1,1])+pi*repmat((2:dim)',[1,no])/dim);            
         case {'F8'}
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = repmat(ps(2,:),[dim-2,1]).*repmat(ps(1,:),[dim-2,1]);             
            pf          = zeros(3,no);
            pf(1,:)     = 1-ps(1,:).*ps(2,:);
            pf(2,:)     = 1-ps(1,:).*(1-ps(2,:));
            pf(3,:)     = ps(1,:);   
            clear s t;
         case {'F9'}
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = repmat(ps(2,:),[dim-2,1]).*repmat(ps(1,:),[dim-2,1]);             
            pf          = zeros(3,no);
            pf(1,:)     = 1-cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
            pf(2,:)     = 1-cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
            pf(3,:)     = 1-sin(0.5*pi*ps(1,:));   
            clear s t;          
    end
end