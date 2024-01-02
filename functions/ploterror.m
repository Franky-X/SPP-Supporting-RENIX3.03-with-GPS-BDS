function ploterror(year,month,date,station,S_path)
% Plot the results
% Noted: This function has the command that starting in R2018b: sgtitle()
close all
% Setting#3
graph_size = [10 10 700 400]; % figure size
gray   = [.5 .5 .5];
red    = [1 0 0];
blue   = [0 0 1];
ltgray = [.8 .8 .8];
% Test
try
    filename = [S_path 'mode_1\pos_m1_' station '_' year '_' month '_' date];
    load(filename)
catch
    error('Calculate positioning first')
end

%% Figure#1 Positioning
main1 = figure('Renderer', 'painters', 'Position', graph_size);
% Starting in R2018b
try
    sgtitle(['Positioning at ' station ' station date:' year '/' month '/' date])
catch
	 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
     set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
     text( 0.5, 0, ['Positioning at ' station ' station date:' year '/' month '/' date], 'FontSize', 14', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end
for mode = 1:4
    % download from save path
    m = num2str(mode);
    filename = [S_path 'mode_' m ...
        '\pos_m' m '_' station '_' year '_' month '_' date];
    load(filename)
    % reference position in LLA
    refpos_lla  = xyz2lla(refpos(1),refpos(2),refpos(3));
%     refpos_lla = [31.025026951 121.439449785 37.0669];
    namepos1 = ['userpos_mode_' m];
    namepos2 = ['disterr_mode_' m];
    namepos3 = ['model_mode_' m];
    eval(['userpos = ' namepos1 ';'])
    eval(['disterr = ' namepos2 ';'])
    eval(['model   = ' namepos3 ';'])
    
    if mode == 1 
        subplot(2,2,1);
        mess = 'mode 1: no atmos correction';
        title(mess)
        hold on
    elseif mode == 2 
        subplot(2,2,2);
        mess = 'mode 2: Iono delay correction';
        title(mess)
        hold on
    elseif mode == 3 
        subplot(2,2,3);
        mess = 'mode 3: Tropo delay correction';
        title(mess)
        hold on
    elseif mode == 4 
        subplot(2,2,4);
        mess = 'mode 4: Tropo + Iono delay correction';
        title(mess)
        hold on
    end
    % position in xyz
    xyz1 = userpos.xyz;
    xyz2 = refpos;
    Med  = nanmedian(userpos.llh); % median llh position
    Med2 = nanmedian(userpos.xyz); % median xyz position    
    refXyz  = (xyz1+(ones(size(xyz1(:,1)))*xyz2))/2;     
    refXyz2 = (Med2+xyz2)/2;
    
    llaDeg  = xyz2lla(refXyz(:,1),refXyz(:,2),refXyz(:,3));
    llaDegm = xyz2lla(refXyz2(1),refXyz2(2),refXyz2(3));
    % prepare matrix
    north  = zeros(length(llaDeg),1);
    east   = zeros(length(llaDeg),1);
    % User
    for i = 1:length(llaDeg)
        Ce2n = RotEcef2Ned(llaDeg(i,1),llaDeg(i,2));
        v1 = Ce2n*(xyz1(i,:)-xyz2)';
        north(i)=v1(1);
        east(i)=v1(2);
        % remove outlier
%         if north(i)>=10 || north(i)<= -10
%             north(i)=0;
%         end
%         if east(i)>=10 || east(i)<= -10
%             east(i)=0;
%         end
    end
    % Ref
    Ce2n2 = RotEcef2Ned(llaDegm(1),llaDegm(2));
    v2 = Ce2n2*(Med2-xyz2)';
    north_m   = v2(1);
    east_m    = v2(2);
    
    % print user position
    h0 = plot(east,north,'.');
    set(h0,'LineStyle','-','LineWidth',0.1,'Color',ltgray)
    % print median position
    lls = sprintf(' [%.6f^o, %.6f^o]',Med(1:2));
    h1 = plot(east_m,north_m,'+k','MarkerSize',18);
    ts1 = ['  median  ',lls];
    ht1 = text(east_m,north_m,ts1,'color','k');
    
    % print reference position
    lls2 = sprintf(' [%.6f^o, %.6f^o]',refpos_lla(1:2));
    h2 = plot(0,0,'+r','MarkerSize',18);
    ts2 = [' true pos ',lls2];
    ht2 = text(0,0,ts2,'color','r');
    
    if north_m <= 0
        set(ht1,'VerticalAlignment','top');     %moves the 'median' label  down
        set(ht2,'VerticalAlignment','bottom');
    else
        set(ht1,'VerticalAlignment','bottom');
        set(ht2,'VerticalAlignment','top');
    end
    axis equal, grid on
    ylabel('North (m)'),xlabel('East (m)')
    clear nedM
    nedM(:,1) = north;
    nedM(:,2) = east;
    pz = true(length(nedM),1);
    % compute error distribution and plot circle
    distM = sqrt(sum(nedM(pz,1:2).^2,2)); % use only finite values here
    medM = nanmedian(distM);
    % plot a circle using 'rectangle' (yes really :)
    hr=rectangle('Position',[-1 -1 2 2]*medM,'Curvature',[1 1]);
    set(hr,'EdgeColor',gray)
    ts3 = sprintf('50%% distribution = %.4f m',medM);
    text(0,medM,ts3,'VerticalAlignment','bottom','Color',gray)
    hold off
    % displayed error
    % 95 percentile
    nedM95(1) = prctile(abs(disterr.NS),95);
    nedM95(2) = prctile(abs(disterr.EW),95);
    nedM95(3) = prctile(abs(disterr.horizontal),95);
    nedM95(4) = prctile(abs(disterr.height),95);
    disp(['Positioning error ' mess])
    disp('=============================================================')
    fprintf('north-south 95 %% error  = %.4f m\n',nedM95(1));
    fprintf('east-west 95 %% error    = %.4f m\n',nedM95(2));
    fprintf('horizontal 95 %% error   = %.4f m\n',nedM95(3));
    fprintf('vertical 95 %% error     = %.4f m\n',nedM95(4));
    % Max
    [m1,tm1] = max(abs(disterr.NS));
    [m2,tm2] = max(abs(disterr.EW));
    [m3,tm3] = max(abs(disterr.horizontal));
    [m4,tm4] = max(abs(disterr.height));
    fprintf('Max north-south error   = %.2f m timeUTC = %.1f hour\n',m1,round(tm1/3600));
    fprintf('Max east-west error     = %.2f m timeUTC = %.1f hour\n',m2,round(tm2/3600));
    fprintf('Max horizontal error    = %.2f m timeUTC = %.1f hour\n',m3,round(tm3/3600));
    fprintf('Max vertical error      = %.2f m timeUTC = %.1f hour\n',m4,round(tm4/3600));
    disp('-------------------------------------------------------------')
    clear disterr userpos north east
end
movegui(main1,'center');

%% Figure#2 Horizontal error
main2 = figure('Renderer', 'painters', 'Position', graph_size);
stt     = 0;    % start UTC time(hour)
stp     = 25;   % stop  UTC time(hour)
int     = 30;   % interval second
% stt     = 0;    % start UTC time(hour)
% stp     = 0.25;   % stop  UTC time(hour)
% int     = 1;   % interval second
len_matrix = (stp - stt) * 3600 - 1;
t_ref = (1:len_matrix)/3600;
% Starting in R2018b
try
    sgtitle(['Horizontal error at ' station ' station date:' year '/' month '/' date])
catch
	 axes( 'Position', [0, 0.95, 1, 0.05] ) ;
     set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
     text( 0.5, 0, ['Horizontal error at ' station ' station date:' year '/' month '/' date], 'FontSize', 14', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end
% 95 percentile
nedM95_hori(1) = prctile(abs(disterr_mode_1.horizontal),95);
nedM95_hori(2) = prctile(abs(disterr_mode_2.horizontal),95);
nedM95_hori(3) = prctile(abs(disterr_mode_3.horizontal),95);
nedM95_hori(4) = prctile(abs(disterr_mode_4.horizontal),95);
hs1 = sprintf('95%% error [Horizontal] = %.3f m',nedM95_hori(1));
hs2 = sprintf('95%% error [Horizontal] = %.3f m',nedM95_hori(2));
hs3 = sprintf('95%% error [Horizontal] = %.3f m',nedM95_hori(3));
hs4 = sprintf('95%% error [Horizontal] = %.3f m',nedM95_hori(4));
% Max
hm1 = sprintf('maximum error [Horizontal] = %.3f m',max(disterr_mode_1.horizontal));
hm2 = sprintf('maximum error [Horizontal] = %.3f m',max(disterr_mode_2.horizontal));
hm3 = sprintf('maximum error [Horizontal] = %.3f m',max(disterr_mode_3.horizontal));
hm4 = sprintf('maximum error [Horizontal] = %.3f m',max(disterr_mode_4.horizontal));

subplot(4,1,1);plot(t_ref,disterr_mode_1.horizontal,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),hs1,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),hm1,'Color',blue,'FontWeight','bold');...
    title('mode 1: no atmos correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on


subplot(4,1,2);plot(t_ref,disterr_mode_2.horizontal,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),hs2,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),hm2,'Color',blue,'FontWeight','bold');...
    title('mode 2: Iono delay correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on

subplot(4,1,3);plot(t_ref,disterr_mode_3.horizontal,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),hs3,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),hm3,'Color',blue,'FontWeight','bold');...
    title('mode 3: Tropo delay correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on

subplot(4,1,4);plot(t_ref,disterr_mode_4.horizontal,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),hs4,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),hm4,'Color',blue,'FontWeight','bold');...
    title('mode 4: Tropo + Iono delay correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on

movegui(main2,'northwest');
%% Figure#3 Vertical error
main3 = figure('Renderer', 'painters', 'Position', graph_size);
% Starting in R2018b
try
    sgtitle(['Vertical error at ' station ' station date:' year '/' month '/' date])
catch
    axes( 'Position', [0, 0.95, 1, 0.05] ) ;
    set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
    text( 0.5, 0, ['Vertical error at ' station ' station date:' year '/' month '/' date], 'FontSize', 14', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end
% 95 percentile
nedM95_verti(1) = prctile(abs(disterr_mode_1.height),95);
nedM95_verti(2) = prctile(abs(disterr_mode_2.height),95);
nedM95_verti(3) = prctile(abs(disterr_mode_3.height),95);
nedM95_verti(4) = prctile(abs(disterr_mode_4.height),95);
vs1 = sprintf('95%% error [Vertical] = %.3f m',nedM95_verti(1));
vs2 = sprintf('95%% error [Vertical] = %.3f m',nedM95_verti(2));
vs3 = sprintf('95%% error [Vertical] = %.3f m',nedM95_verti(3));
vs4 = sprintf('95%% error [Vertical] = %.3f m',nedM95_verti(4));

% Max
vm1 = sprintf('maximum error [Vertical] = %.3f m',max(disterr_mode_1.height));
vm2 = sprintf('maximum error [Vertical] = %.3f m',max(disterr_mode_2.height));
vm3 = sprintf('maximum error [Vertical] = %.3f m',max(disterr_mode_3.height));
vm4 = sprintf('maximum error [Vertical] = %.3f m',max(disterr_mode_4.height));

subplot(4,1,1);plot(t_ref,disterr_mode_1.height,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),vs1,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),vm1,'Color',blue,'FontWeight','bold');...
    title('mode 1: no atmos correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on

subplot(4,1,2);plot(t_ref,disterr_mode_2.height,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),vs2,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),vm2,'Color',blue,'FontWeight','bold');...
    title('mode 2: Iono delay correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on

subplot(4,1,3);plot(t_ref,disterr_mode_3.height,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),vs3,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),vm3,'Color',blue,'FontWeight','bold');...
    title('mode 3: Tropo delay correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on

subplot(4,1,4);plot(t_ref,disterr_mode_4.height,'.');yt = get(gca,'ylim');...
    text(16,yt(2)+(diff(yt)*0.3),vs4,'Color',red,'FontWeight','bold');...
    text(16,yt(2)+(diff(yt)*0.1),vm4,'Color',blue,'FontWeight','bold');...
    title('mode 4: Tropo + Iono delay correction');...
    xlabel('Time(UTC)');ylabel('error (m)');xlim([0 stp]);grid on
movegui(main3,'southwest');

%% Figure#4 Satellite and receiver bias
main4 = figure('Renderer', 'painters', 'Position', graph_size);
try
    sgtitle(['instrument biases at ' station ' station date:' year '/' month '/' date])
catch
	axes( 'Position', [0, 0.95, 1, 0.05] ) ;
    set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
    text( 0.5, 0, ['instrument biases at ' station ' station date:' year '/' month '/' date], 'FontSize', 14', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end
c  = 299792458;                         %   light speed = 299792458 m/s
subplot(2,1,1);
plot(t_ref,userpos_mode_1.bs./(10^-9),'.');
title('Satellite clock delay');
xlabel('Time(UTC)');
ylabel('delay (ns)');
yt = get(gca,'ylim');

grid on;xlim([0 stp]);ylim([0 yt(2)+(diff(yt)*0.15)]);
subplot(2,1,2);
hold on
plot(t_ref,userpos_mode_1.br./(c*10^-9),'.');
plot(t_ref,userpos_mode_2.br./(c*10^-9),'.');
plot(t_ref,userpos_mode_3.br./(c*10^-9),'.');
plot(t_ref,userpos_mode_4.br./(c*10^-9),'.');
hold off
legend('mode#1','mode#2','mode#3','mode#4')
title('Receiver clock delay');
xlabel('Time(UTC)');
ylabel('delay (ns)');
yt = get(gca,'ylim');

grid on;xlim([0 stp]);ylim([yt(1)-1 yt(2)+(diff(yt)*0.15)]);
movegui(main4,'northeast');

%% Figure#5 delay model
main5 = figure('Renderer', 'painters', 'Position', graph_size);
try
    sgtitle(['Atmospheric delay model at ' station ' station date:' year '/' month '/' date])
catch
	axes( 'Position', [0, 0.95, 1, 0.05] ) ;
    set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
    text( 0.5, 0, ['Atmospheric delay model at ' station ' station date:' year '/' month '/' date], 'FontSize', 14', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end
c  = 299792458;                         %   light speed = 299792458 m/s
subplot(2,1,1);
plot(t_ref,model.tropo./(c*10^-9),'.');
title('Tropospheric delay');
xlabel('Time(UTC)');
ylabel('delay (ns)');
yt = get(gca,'ylim');
text(18,yt(2)+(diff(yt)*0.1),'CSSRG Laboratory@KMITL, Thailand.','Color',gray,'FontSize',6);
grid on;xlim([0 stp]);ylim([0 yt(2)+(diff(yt)*0.15)]);
subplot(2,1,2);
plot(t_ref,model.iono./(c*10^-9),'.');
title('Ionospheric delay');
xlabel('Time(UTC)');
ylabel('delay (ns)');
yt = get(gca,'ylim');
text(18,yt(2)+(diff(yt)*0.1),'CSSRG Laboratory@KMITL, Thailand.','Color',gray,'FontSize',6);
grid on;xlim([0 stp]);ylim([0 yt(2)+(diff(yt)*0.15)]);
movegui(main5,'southeast');

end
