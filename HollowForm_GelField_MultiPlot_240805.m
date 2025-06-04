%% Simulation codes

rng('shuffle')

N_component = 3; % Number of components, 1: p53, 2: Random DNA, 3: p21 DNA

% Parameter
kT = 1.0;
kappa = 1;

% Size
L = 101; % Number of grids for one side
L2 = L^2; % Number of grids of the system
dL = 1; % Length between grids

% Mobility
D = ones(L,L);

% Interaction for N_component = 3
chi = zeros(N_component);
chi(1,2) = -5;  chi(2,1) = chi(1,2);
chi(1,3) = -5;  chi(3,1) = chi(1,3);

lambda_solvent = 1.8;

% Initial condition
phi_mat = zeros(L,L,N_component);
% for idx_component = 1:3
%     phi_mat(:,:,idx_component) = 0.5 - 0.01 + 0.02*rand(L,L);
% end
phi_mat(:,:,1) = 0.084*ones(L,L);%phi_out_s 
phi_mat(:,:,2) = 0.084*ones(L,L);%psi_out_s 
R0 = 25;
for i = 1:L
    for j = 1:L
        if (i-(L-1)/2-1)^2+(j-(L-1)/2-1)^2 <= R0^2
            phi_mat(i,j,1) = 1.62;
            phi_mat(i,j,2) = 1.62;
        end
    end
end
phi_mat(:,:,3) = 0.001 + 0.001*rand(L,L);

gel_field = zeros(L,L);

% Simulation time info
t_total = 10000;
dt = 0.005;
dt_rec = 10;

t_add = 500;%1000;
add_judge = false;  

N_rec = round(t_total/dt_rec);
t_rec = zeros(1,N_rec);
phi_mat_rec = zeros(L2,N_rec,N_component);

t = 0;
t_next = 0;
count = 1;
tic;
while t <= t_total
    if t >= t_next
        disp(t);

        t_rec(count) = t;
        for idx_component = 1:N_component
            phi_temp = phi_mat(:,:,idx_component);
            phi_mat_rec(:,count,idx_component) = phi_temp(:);
        end

        count = count + 1;
        t_next = t_next + dt_rec;
    end

    % Add p21 DNA as phi_3
    if ~add_judge && t >= t_add
        phi_3_add_mean = 0.4;
        phi_3_add = phi_3_add_mean + 0.02*rand(L,L);
        phi_3_add(phi_mat(:,:,1)>1) = 0.001;
        phi_mat(:,:,2) = phi_mat(:,:,2) + phi_3_add;
        add_judge = true;
    end

    % Time step change
    if t >= t_add
        dt = 0.001;
    end
    if t > t_add + 1000
        dt = 0.005;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    phi_tot = sum(phi_mat,3);

    dphi_mat = zeros(L,L,N_component);
    for idx_component_i = 1:N_component
        phi = phi_mat(:,:,idx_component_i);
        
        U = zeros(L,L);
        for idx_component_j = 1:N_component
            U = U + chi(idx_component_i,idx_component_j).*phi_mat(:,:,idx_component_j);
        end
        
        delta2phi = delta2(phi,dL);
        [dmu_phi_x,dmu_phi_y] = delta(log(phi) + lambda_solvent*phi.^2 + U - kappa.*delta2phi, dL);
        v_x = -D.*dmu_phi_x + sqrt(2*D*dt).*randn(L,L);
        v_y = -D.*dmu_phi_y + sqrt(2*D*dt).*randn(L,L);
        dphi = dt*(-div(phi.*v_x,phi.*v_y,dL));

        dphi_mat(:,:,idx_component_i) = dphi;
    end

    % Gel field dynamics
    M_gel = 2; % M_g in the text
    phi_thresh = 1.0;
    gel_judge = (phi_mat(:,:,1).*phi_mat(:,:,3)-phi_thresh)./(2.0-phi_thresh); % a in the text
    dgel_field = dt*(M_gel.*(gel_judge.*gel_field-gel_field.^2)) + sqrt(2*M_gel/1e5*dt).*randn(L,L);
    gel_field = gel_field + dgel_field;
    % gel_field(gel_field <= 0) = 0;
    gel_field = abs(gel_field);
    gel_field_c = 0.3;
    D = exp(-gel_field/gel_field_c);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = t + dt;

    phi_mat = phi_mat + dphi_mat;

    phi_check = sum(phi_mat,3);
    if max(abs(imag(phi_check(:)))) > 0 || ~(max(abs(phi_check(:))) > 0)
        disp('Result Not Real');
        break;
    end

end
toc;

% Plot results
for i = 5:5:N_rec
    imagesc(reshape(phi_mat_rec(:,i,2),L,L));colormap sky;colorbar;clim([0 1]);
    hold on;
    title(['t=',num2str(i*dt_rec)])
    hold off;
    axis('square')
    set(gca,'Fontsize',15)
    pause(0.1)
end

color_base = lines(7); color_base = color_base(2,:);
color_list_plot = [linspace(1,color_base(1),200)',linspace(1,color_base(2),200)',linspace(1,color_base(3),200)'];
for i = 5:5:N_rec
    imagesc(reshape(phi_mat_rec(:,i,3),L,L));colormap(color_list_plot);colorbar;clim([0 1]);
    hold on;
    title(['t=',num2str(i*dt_rec)])
    hold off;
    axis('square')
    set(gca,'Fontsize',15)
    pause(0.1)
end




%% Make 1D video
v = VideoWriter('output1d','MPEG-4');
open(v);
color_red = [230,0,18]/255;
color_green = [0,153,68]/255;
color_purple = [126,49,142]/255;
color_list = [color_red;color_green;color_purple];
for i = 1:300
    cla;
    ti = i*dt_rec;
    for idx_component = 1:N_component       
        phi_plot = reshape(phi_mat_rec(:,i,idx_component),L,L);
        plot(phi_plot((L+1)/2,:),'Linewidth',1.5,'Color',color_list(idx_component,:))
        hold on
    end
    set(gca,'FontSize',15)
    xlim([0 L])
    ylim([0 3])
    xlabel('Position','Fontsize',20)
    ylabel('Concentration','Fontsize',20)
    legend('Protein','Random DNA','p21')
    % title(['t = ',num2str(i*dt_rec)],'Fontsize',15)

    if ti == t_add
        title(['t = ',num2str(i*dt_rec), ', Add 3\timesp21 DNA'],'Fontsize',15)
        for count=1:60
        frame = getframe(gcf);
        tmp = frame2im(frame);
        writeVideo(v,tmp);
        end
    elseif mod(i,1)==0
        title(['t = ',num2str(i*dt_rec)],'Fontsize',15)
        frame = getframe(gcf);
        tmp = frame2im(frame); 
        writeVideo(v,tmp);
    end

end
close(v)




%% Make 2D video
v = VideoWriter('output2d','MPEG-4');
open(v);
for i = 1:500
    cla;
    ti = i*dt_rec;
    set(gcf,'position',[460 527 800 233])
    ax = gobjects(1,3);
    for idx_component = 1:N_component
        ax(idx_component) = subplot(1,3,idx_component);
        color_red = [230,0,18]/255;
        color_green = [0,153,68]/255;
        color_purple = [126,49,142]/255;
        color_list = [color_red;color_green;color_purple];
        % color_base = lines(7); color_base = color_base(idx_component,:);
        color_base = color_list(idx_component,:);
        color_list_plot = [linspace(0,color_base(1),100)',linspace(0,color_base(2),100)',linspace(0,color_base(3),100)'];
        imagesc((1:L)*dL,(1:L)*dL,reshape(phi_mat_rec(:,i,idx_component),L,L));
        colorbar; colormap(ax(idx_component),color_list_plot); clim([0 2]);
        
        set(gca,'XTick',[])
        set(gca,'YTick',[])

        axis('square')
        set(gca,'position',[0.0400+0.33*(idx_component-1), 0.050, 0.25, 0.8])
        set(gca,'Fontsize',12)
    end
    
    sgtitle(['t = ',num2str(i*dt_rec)],'Fontsize',15)

    if ti == t_add
        sgtitle(['t = ',num2str(i*dt_rec), ', Add 3\timesp21 DNA'],'Fontsize',15)
        for count=1:60
        frame = getframe(gcf);
        tmp = frame2im(frame);
        writeVideo(v,tmp);
        end
    elseif mod(i,1)==0
        sgtitle(['t = ',num2str(i*dt_rec)],'Fontsize',15)
        frame = getframe(gcf);
        tmp = frame2im(frame); 
        writeVideo(v,tmp);
    end

end
close(v)


%% Plot 2D Figure

color_red = [230,0,18]/255;
color_green = [0,153,68]/255;
color_purple = [126,49,142]/255;
color_list = [color_red;color_green;color_purple];

i = 500;

figure;
idx_component = 1;
color_base = color_red;
color_list_plot = [linspace(0,color_base(1),100)',linspace(0,color_base(2),100)',linspace(0,color_base(3),100)'];
imagesc((1:L)*dL,(1:L)*dL,reshape(phi_mat_rec(:,i,idx_component),L,L));
colorbar; colormap(color_list_plot); clim([0 2]);
axis('square')
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'FontSize',20)
figure;
idx_component = 2;
color_base = color_green;
color_list_plot = [linspace(0,color_base(1),100)',linspace(0,color_base(2),100)',linspace(0,color_base(3),100)'];
imagesc((1:L)*dL,(1:L)*dL,reshape(phi_mat_rec(:,i,idx_component),L,L));
colorbar; colormap(color_list_plot); clim([0 2]);
axis('square')
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'FontSize',20)
figure;
idx_component = 3;
color_base = color_purple;
color_list_plot = [linspace(0,color_base(1),100)',linspace(0,color_base(2),100)',linspace(0,color_base(3),100)'];
imagesc((1:L)*dL,(1:L)*dL,reshape(phi_mat_rec(:,i,idx_component),L,L));
colorbar; colormap(color_list_plot); clim([0 2]);
axis('square')
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'FontSize',20)


%% Plot 1D Figure

color_red = [230,0,18]/255;
color_green = [0,153,68]/255;
color_purple = [126,49,142]/255;
color_list = [color_red;color_green;color_purple];

i = 500;
for idx_component = 1:N_component
    phi_plot = reshape(phi_mat_rec(:,i,idx_component),L,L);
    plot(0:L-1,phi_plot((L+1)/2,:),'Linewidth',1.5,'Color',color_list(idx_component,:))
    hold on
end
set(gca,'FontSize',25)
ylim([0 3])
xlabel('Position','Fontsize',30)
ylabel('Concentration','Fontsize',30)
%legend('p53 protein','Random DNA','3\timesp21 DNA')
title(['t = ',num2str(i*dt_rec)])
xlim([0 L])
set(gca,'XTick',[0 50 100])


%% Plot Min Eigenvalue inside DPIC
GaussCurve_list = zeros(1,N_rec);
phi_center_list = zeros(1,N_rec);
psi_center_list = zeros(1,N_rec);
xi_center_list = zeros(1,N_rec);
eigmin_list = zeros(1,N_rec);
det_list = zeros(1,N_rec);
for i = 1:N_rec
    phi_plot = reshape(phi_mat_rec(:,i,1),L,L);
    psi_plot = reshape(phi_mat_rec(:,i,2),L,L);
    xi_plot = reshape(phi_mat_rec(:,i,3),L,L);

    center_idx1 = (L+1)/2;
    center_idx2 = (L+1)/2;
    phi_center = phi_plot(center_idx1,center_idx2);
    psi_center = psi_plot(center_idx1,center_idx2);
    xi_center = xi_plot(center_idx1,center_idx2);

    phi_center_list(i) = phi_center;
    psi_center_list(i) = psi_center;
    xi_center_list(i) = xi_center;

    GaussCurve_list(i) = GaussCurve(phi_center,psi_center,chi(1,2),lambda_solvent);

    dd_mat = [1/phi_center+2*lambda_solvent*phi_center,chi(1,2); ...
        chi(2,1),1/psi_center+2*lambda_solvent*psi_center];
    eigmin_list(i) = min(eig(dd_mat));
    det_list(i) = det(dd_mat);
end

figure;
yyaxis left
plot((1:N_rec)*dt_rec,eigmin_list,'LineWidth',1.5)
hold on
plot((1:N_rec)*dt_rec,zeros(1,N_rec),'b--','LineWidth',1.5)
set(gca,'FontSize',15)
xlabel('Time','Fontsize',20)
ylabel('Min Eigenvalue','Fontsize',20)

yyaxis right
color_red = [230,0,18]/255;
color_green = [0,153,68]/255;
color_purple = [126,49,142]/255;
color_list = [color_red;color_green;color_purple];
p1 = plot((1:N_rec)*dt_rec,phi_center_list,'Color',color_list(1,:),'LineWidth',1.5);
p2 = plot((1:N_rec)*dt_rec,psi_center_list,'Color',color_list(2,:),'LineWidth',1.5);
ylabel('Concentration','Fontsize',20)

legend([p1,p2],'p53 protein','Random DNA')

xlim([400 2200])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Demanded function
function divA = div(Ax,Ay,d)
    Axp1 = circshift(Ax,[0,-1]);
    Axm1 = circshift(Ax,[0,1]);
    Axp2 = circshift(Ax,[0,-2]);
    Axm2 = circshift(Ax,[0,2]);
    dAx = (-Axp2+Axm2+8.*Axp1-8.*Axm1)/(12*d);
    Ayp1 = circshift(Ay,[-1,0]);
    Aym1 = circshift(Ay,[1,0]);
    Ayp2 = circshift(Ay,[-2,0]);
    Aym2 = circshift(Ay,[2,0]);
    dAy = (-Ayp2+Aym2+8.*Ayp1-8.*Aym1)/(12*d);
    divA = dAx + dAy;
end
function [dAx,dAy] = delta(A,d)
    Axp1 = circshift(A,[0,-1]);
    Axm1 = circshift(A,[0,1]);
    Axp2 = circshift(A,[0,-2]);
    Axm2 = circshift(A,[0,2]);
    dAx = (-Axp2+Axm2+8.*Axp1-8.*Axm1)/(12*d);
    Ayp1 = circshift(A,[-1,0]);
    Aym1 = circshift(A,[1,0]);
    Ayp2 = circshift(A,[-2,0]);
    Aym2 = circshift(A,[2,0]);
    dAy = (-Ayp2+Aym2+8.*Ayp1-8.*Aym1)/(12*d);
end
function d2A = delta2(A,d)
    Axp1 = circshift(A,[0,-1]);
    Axm1 = circshift(A,[0,1]);
    Ayp1 = circshift(A,[-1,0]);
    Aym1 = circshift(A,[1,0]);
    Axp1yp1 = circshift(A,[-1,-1]);
    Axm1yp1 = circshift(A,[-1,1]);
    Axp1ym1 = circshift(A,[1,-1]);
    Axm1ym1 = circshift(A,[1,1]);
    d2A = (4*Axm1+4*Axp1+4*Aym1+4*Ayp1+Axp1yp1+Axm1yp1+Axp1ym1+Axm1ym1-20*A)/(6*d^2);
end


function f = f_double(phi,psi,lambda,chi)
    f = phi.*log(phi) + psi.*log(psi) + lambda/3*phi.^3 + lambda/3*psi.^3 + chi.*phi.*psi;
end
function df = df_phi(phi,psi,lambda,chi)
    df = log(phi) + 1 + lambda*phi.^2 + chi*psi;
end
function df = df_psi(phi,psi,lambda,chi)
    df = log(psi) + 1 + lambda*psi.^2 + chi*phi;
end
function Pi = Pi_double(phi,psi,lambda,chi)
    Pi = phi.*df_phi(phi,psi,lambda,chi) + psi.*df_psi(phi,psi,lambda,chi) - f_double(phi,psi,lambda,chi);
end
function K = GaussCurve(phi,psi,lambda,chi)
    du_phi = df_phi(phi,psi,lambda,chi);
    du_psi = df_psi(phi,psi,lambda,chi);
    ddu_phiphi = 1./phi + 2*lambda.*phi;
    ddu_psipsi = 1./psi + 2*lambda.*psi;
    ddu_phipsi = chi;
    K = (ddu_phiphi.*ddu_psipsi-ddu_phipsi.^2)./((1+du_phi.^2+du_psi.^2).^2);
end