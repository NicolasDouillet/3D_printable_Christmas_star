function [V, T] = printable_3D_Christmas_star(nb_samples, option_display)
%
% Author & support : nicolas.douillet (at) free.fr, 2016-2022
%
%
% Input arguments
%
% - nb_samples : positive integer scalar double, sampling > 2.
%
% - option_display : logical *true (1) / false (0).
%
%
% Output arguments
%
%       [|  |  | ]
% - V = [Vx Vy Vz], real matrix double, the point set. Size(V) = [nb_vertices,3].
%       [|  |  | ]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]


% Input parsing
if nargin < 2
    
    option_display = true;
    
    if nargin < 1
        
        nb_samples = 32;
        
    end
    
end


% Geometric parameters
phi_n = 0.5*(1+sqrt(5));
a = 2/sqrt(phi_n*sqrt(5)); % edge length
e = 0.2; % radial sinusoidal variation amplitude

% Z axis 2*pi/5 rotation matrix 
Mrz = [cos(0.4*pi) -sin(0.4*pi) 0;
       sin(0.4*pi)  cos(0.4*pi) 0;
       0            0           1];

centre_angle = 2*asin(1/sqrt(phi_n*sqrt(5)));
           
% 1st equilateral triangle
V0 = [0 0 1]';
V1 = [sin(centre_angle) 0 cos(centre_angle)]';
V2 = Mrz*V1;

% Bidirectional (u,v) sampling + compute corresponding squared distances vector
[P0,T0] = sample_and_shape_triangle(V0,V1,V2,a,nb_samples,e);

% Replicate / rotation -> upper crown
Pu = P0';

for k = 1:4    
    Pu = cat(1,Pu,(Mrz^k*P0)');    
end

% Lower base triangle with /O symetry
V3 = -V0;
V4 = -V1;
V5 = -V2;

P1 = sample_and_shape_triangle(V3,V4,V5,a,nb_samples,e);

% Lower crown
Pl = P1';

for k = 1:4    
    Pl = cat(1,Pl,(Mrz^k*P1)');    
end

P = cat(1,Pu,Pl);

% 1st belt triangle
V6 = V1;
V7 = V2;
V8 = Mrz^2*V5;

P2 = sample_and_shape_triangle(V6,V8,V7,a,nb_samples,e);
P2 = P2';

% 2nd belt triangle
V9  = -V6;
V10 = -V7;
V11 = -V8;

P3 = sample_and_shape_triangle(V9,V10,V11,a,nb_samples,e);
P3 = P3';

% Full belt = centre crown
P4 = cat(1,P2,P3);
Pc = P4;

for k = 1:4    
    Pc = cat(1,Pc,(Mrz^k*P4')');    
end

V = cat(1,P,Pc);

% Triangulation
T = T0;
S0 = size(P0,2);

for k = 1:19    
    T = cat(1,T,T0+k*S0);
end
              

% Preparation for 3D printing

% Remove duplicated triangles        
T = sort(T,2);
T = unique(T,'rows','stable');        
        
% Remove duplicated vertices
[V,T] = remove_duplicated_vertices(V,T);


TRI = triangulation(T,V(:,1),V(:,2),V(:,3));


% Display
if option_display
    
    figure;
    set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0]);
    trisurf(TRI), shading faceted;
    colormap([1 1 0]);
    axis square, axis equal, axis off;
    ax = gca;
    ax.Clipping = 'off';
    camlight left;
    view(3);
    zoom(1.8);
    
end


end % printable_3D_Christmas_star


% sample_and_shape_triangle subfunction
function [V, T] = sample_and_shape_triangle(V1, V2, V3, a, nbstep, e)
%
% Author & support : nicolas.douillet (at) free.fr, 2016-2022


% (V1V2, V1V0) basis
u = (V3 - V2);
v = (V1 - V2);
nu = u / norm(u);
nv = v / norm(v);

stepu = a / nbstep;
stepv = a / nbstep;

% Sampling & points creation
V = [];

for m = 0:nbstep
    
    for n = 0:nbstep
        
        if(m+n <= nbstep) % in (V1,V2,V3) triangle conditions ; indices # nb segments                        
       
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;                                          
       
            % Basis (V2V1, V2V3)
            M = V2 + tv;            
            min_dst = min([norm(tv) norm(tv+V2-V1) norm(tv+V2-V3)]);
            
            % Radial vector            
            rv = M + e*cos(2*pi*min_dst/norm(u))*M/norm(M);       
            V = cat(2,V,rv);
            
        end
           
    end

end


% Index triplets list construction
T = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while p <= nbstep^2 && row_length > 1
    
     i = p;
    
    if p < 2 % "right" triangle serie only
        
        while i < cum_row_length
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while i < cum_row_length % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i+1 i-row_length]; % + upside-down triangles serie
            row_idx = row_idx + 1;
            
            i = i +1;
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end


end % sample_and_shape_triangle


% Remove_duplicated_vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)


tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);


end % remove_duplicated_vertices