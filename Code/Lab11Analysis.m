% PIV Data Analysis
% AER E 344 Spring 2024 - Section 3 Group 3
clear; clc; close all;

%% Constants
rho = 1.225; % [kg/m^3]
figure_dir = "../Figures/";
data_dir = "./Data/";
zip_dir = "./Data/AER E 344 PIV LAB 11 Data.zip";
file_dir = "./PIV LAB/";
space = 5.209344;
chord = 101; % [mm]

%% Instantaneous PIV Measurements
instant_title = ["B00202.dat","B00245.dat"];
source = ["AOA8","AOA16"];
source_title_str = ["AoA 8°","AoA 16°"];
for t = 1:length(instant_title)
    data = readtable(data_dir+instant_title(t),"VariableNamingRule", ...
        "modify","NumHeaderLines",3);
    data.Properties.VariableNames = {'x','y','Vx','Vy'};
    % Initalize matrices for vorticity calculations.
    X = min(data.x):space:max(data.x)+1;
    Y = max(data.y):-space:min(data.y)-1;
    U = zeros(64,64);
    V = zeros(64,64);
    % Put .dat columns into matrix format with 1,1 in the top left
    count = 1;
    for i = 1:64
        for j = 1:64
            U(i,j) = data.Vx(count);
            V(i,j) = data.Vy(count);
            count = count +1;
        end
    end
    % Plot velocity field on top of vorticity
    figure;
    colormap("jet");
    [curlz,cav] = curl(X,Y,U,V);
    pcolor(X,Y,curlz);
    hold on;
    quiver(data.x,data.y,data.Vx,data.Vy,"Color","black","AutoScale", ...
        "off");
    title_str = "Instantaneous PIV Measurements";
    title(title_str + " (" + source_title_str(t) + ")");
    xlabel("X Position [mm]");
    ylabel("Y Position [mm]");
    fontname("Times New Roman");
    fontsize(12, "points");
    c = colorbar;
    c.Label.String = "Vorticity [s^{-1}]";
    saveas(gcf, figure_dir + title_str + " - " + source(t) + ".svg");
    %% Slide 30 Method (does not look right compared to TecPlot)
    % figure
    % [PartialVx,PartialVy] = gradient(V,space,-space);
    % [PartialUx,PartialUy] = gradient(U,space,-space);
    % w = PartialVy - PartialUx;
    % pcolor(X,Y,w);
end

%% Ensemble Average Measurements 
unzip(zip_dir);
dirs = ["AOA4","AOA8","AOA12","AOA16"];
aoa_str = ["AoA 4°","AoA 8°","AoA 12°","AoA 16°"];
for d = 1:length(dirs)
    sum = table(data.x,data.y,zeros(4096,1),zeros(4096,1),zeros(4096,1));
    sum.Properties.VariableNames = {'x','y','U','V','TKE'};
    % Iterate through all 250 frames to get ensemble average
    for f = 1:250
        filename = file_dir+dirs(d)+"/B00"+num2str(f,"%03d")+".dat";
        data = readtable(filename,"VariableNamingRule","modify", ...
            "NumHeaderLines",3);
        data.Properties.VariableNames = {'x','y','Vx','Vy'};
        sum.U = sum.U + data.Vx;
        sum.V = sum.V + data.Vy;
    end
    sum.U = sum.U/250;
    sum.V = sum.V/250;
    count = 1;
    for i = 1:64
        for j = 1:64
            U(i,j) = sum.U(count);
            V(i,j) = sum.V(count);
            count = count +1;
        end
    end
    % Plot ensemble average data
    figure;
    colormap("jet");
    [curlz,cav] = curl(X,Y,U,V);
    pcolor(X,Y,curlz);
    hold on;
    quiver(sum.x,sum.y,sum.U,sum.V,"Color","black","AutoScale","off");
    title_str = "Ensemble PIV Measurements";
    title(title_str + " (" + aoa_str(d) + ")");
    xlabel("X Position [mm]");
    ylabel("Y Position [mm]");
    fontname("Times New Roman");
    fontsize(12, "points");
    c = colorbar;
    c.Label.String = "Vorticity [s^{-1}]";
    saveas(gcf, figure_dir + title_str + " - " + dirs(d) + ".svg");

    %% Turbulence Distributions
    avg_U = mean(sum.U);
    avg_V = mean(sum.V);
    prime = table(data.x,data.y,zeros(4096,1),zeros(4096,1));
    prime.Properties.VariableNames = {'x','y','u','v'};
    for f = 1:250
        filename = file_dir+dirs(d)+"/B00"+num2str(f,"%03d")+".dat";
        data = readtable(filename,"VariableNamingRule","modify", ...
            "NumHeaderLines",3);
        data.Properties.VariableNames = {'x','y','Vx','Vy'};
        prime.u = prime.u + ((data.Vx - avg_U).^2)/250;
        prime.v = prime.v + ((data.Vy - avg_V).^2)/250;
    end
    prime.u = sqrt(prime.u);
    prime.v = sqrt(prime.v);
    sum.TKE = .5 * rho * (prime.u.^2 + prime.v.^2);
    TKE = zeros(64,64);
    % RE = zeros(64,64);
    % Put .dat columns into matrix format with 1,1 in the top left
    count = 1;
    for i = 1:64
        for j = 1:64
            TKE(i,j) = sum.TKE(count);
            %RE(i,j) = sum.tau(count);
            count = count +1;
        end
    end
    figure;
    colormap("jet");
    pcolor(X,Y,TKE);
    hold on;
    quiver(sum.x,sum.y,sum.U,sum.V,"Color","black","AutoScale","off");
    title_str = "Ensemble Turbulent Kinetic Energy";
    title(title_str + "  (" + aoa_str(d) + ")");
    xlabel("X Position [mm]");
    ylabel("Y Position [mm]");
    fontname("Times New Roman");
    fontsize(12, "points");
    c = colorbar;
    c.Label.String = "Total Kinetic Energy";
    saveas(gcf, figure_dir + title_str + " - " + dirs(d) + ".svg");

    %% Reynolds Stress (add tau column to sum table)
    % sum.tau = -rho* prime.u .* prime.v;
    % figure;
    % colormap("jet");
    % pcolor(X,Y,RE);
    % hold on;
    % quiver(sum.x,sum.y,sum.U,sum.V,"Color","black","AutoScale","off");
    % title_str = "Reynold"s Stress Distribution (" + dirs(d) + ")";
    % title(title_str);
    % colorbar;

    %% Wake Profile at 1/2 chord downstream
    figure;
    % Estimating half chord is at x position of 0 (about index 31)
    % Half Chord downstream is about -100 millimeters (about index 11)
    plot(Y,sqrt(U(:,11).^2 + V(:,11).^2));
    title_str = "Wake Profile at Half Chord";
    title(title_str + "  (" + aoa_str(d) + ")");
    xlabel("Y Position [mm]");
    ylabel("Velocity Magnitude [m/s]");
    fontname("Times New Roman");
    fontsize(12, "points");
    grid on;
    saveas(gcf, figure_dir + title_str + " - " + dirs(d) + ".svg");

    %% Write to File
    filename = "./Data/Ensemble/"+dirs(d)+".dat";
    fid = fopen(filename, "wt");
    fprintf(fid, 'TITLE = "%s"\n',dirs(d));
    fprintf(fid, 'VARIABLES = "x", "y", "Vx", "Vy", "TKE"\n');
    fprintf(fid, 'ZONE T="Frame 0", I=64, J=64\n');
    fclose(fid);
    writetable(sum, filename,"Delimiter","\t","WriteVariableNames", ...
        false,"WriteMode","append");
end

%% Clean Up
rmdir(file_dir,"s");
