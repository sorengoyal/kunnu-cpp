clear all; close all; clc
%*************************************************************************
files = dir('*.csv');
new_files = cell(1,length(files));
for i=1:length(files)
   new_files{i} = csvread(files(i).name,1,0);
   new_files{i} = unique(new_files{i},'rows');
   new_files{i} = sortrows(new_files{i},[5 4]);
end
firstout = new_files(1);  % First file
out = cell2mat(firstout); % Converting cell to matrix 
[m,n0] = size(out);  % Get the total no. of nodes on mesh as 'm'
[m1,n01] = size(files);    % Total CSV files in the folder
t = 3;        % No. of time steps saved for each run (enter manually)
n2 = m1/t;  % n2 is the no. of samples = total no. of csv files/no. of 'dt'
load xi_n.dat
xi_n = sort(xi_n, 'ascend');
xi_1 = reshape(xi_n(:,1),[1,length(xi_n)]); % Extract the vector of std. normal variable
Samples = length(xi_1);
% Script Check
if Samples == n2
   disp('All Ok')
else
   disp('xi_1 and n2 not matching')
end
%Bring X and Y coordinate Columns
X1 = out(:, 4);  Y1 = out(:, 5);  TRI = delaunay(X1,Y1);
%*************************************************************************
% For all given .csv files in the folder, we extract the column '1' and '2'
% separately, and store them in x_sol, and y_sol respectively. They both
% contain the solution vectors for all samples and all time steps (sorted
% in an ascending order as per the file name), for directions 'x' and 'y'
x_sol = cell(1,m1);
for i=1:m1
    tx = new_files(i);
    tempx = cell2mat(tx);
    x_sol{i}=tempx(:,1);  % create cell from a part of matrix in other cell
end
x_sol1 = cell2mat(x_sol); % Converting cell to matrix 
y_sol = cell(1,m1);
for j=1:m1
    ty = new_files(j);
    tempy = cell2mat(ty);
    y_sol{j}=tempy(:,2); 
end
y_sol1 = cell2mat(y_sol); % Converting cell to matrix 
%*************************************************************************
% The description here is same as that given right below for y_solf, expect
% that here we create the final solution matrix for x-direction
x_solf = [];
for kx = 1:t
    x_sol2 = [];
    for gx = kx:t:m1
        x_sol2 = [x_sol2 x_sol1(:,gx)];         
    end
    x_solf = [x_solf;x_sol2];
end
% Here we create an empty matrix y_solf, which represents the final
% solution matrix in 'y' direction. The variable 'ky' is a representative
% of time step, and 'm1' is the total no. of CSV files. Therefore, what
% 'gy' represents basically is all the solution vectors (generated for
% n2 random samples), for a single time step,i.e., value of ky. When ky=1
% gy will extract all solution vectors for n2 samples for time step 1. The
% final structure of y_solf will be (m*t,n2), where all solution vectors
% for 1st time step will be stored side by side in y_solf(1:m,n2). Below
% that will be the solutions for time step 2 stored in y_solf(m:2m,n2).
y_solf = [];
for ky = 1:t
    y_sol2 = [];
    for gy = ky:t:m1
        y_sol2 = [y_sol2 y_sol1(:,gy)];         
    end
    y_solf = [y_solf;y_sol2];
end
% Hermite Polynomial generated from file PC_examples_1D for 2RV + 3rd order
% It is essential that xi_1 and xi_2 are same as the std. normal variable
% generated at the time of running the model through the bash script. That
% vector has to be imported here from that bash script.
Psi_0 = ones(1,n2); % 1st Gauss Hermite Polynomial for single RV
Psi_1 = xi_1;       % 2nd Gauss Hermite Polynomial for single RV
Psi_2 = (xi_1.^2 - 1);
Psi_3 = (xi_1.^3 - 3*xi_1);
% Coefficients for solution in x-direction; Size of x_solf is (m*t,n2),
% therefore size of coefficient vector is (m*t,1), because each value of
% coefficient is calculated across each row of x_solf. Hence, a vector of
% coefficients is calulated for all nodes on the mesh and all time steps.
% The structure of this array will be such that it contains the coefficient
% vector for first time step from (1:m,1), then next time step from (m:2m,1)
ax0 = zeros(m*t,1); 
ax1 = zeros(m*t,1); 
ax2 = zeros(m*t,1); 
ax3 = zeros(m*t,1); 
% As explained above, here we finally calculate the coefficients in a loop
% as per the MC based PCE formula. Here 'n2' is the total no. of samples,
% and also the no. of columns of matrices - x_solf, and later y_solf
for ii = 1:m*t
 ax0(ii) = (sum(Psi_0(:).*x_solf(ii,:)')/n2)/(sum(Psi_0(:).*Psi_0(:))/n2); 
 ax1(ii) = (sum(Psi_1(:).*x_solf(ii,:)')/n2)/(sum(Psi_1(:).*Psi_1(:))/n2); 
 ax2(ii) = (sum(Psi_2(:).*x_solf(ii,:)')/n2)/(sum(Psi_2(:).*Psi_2(:))/n2); 
 ax3(ii) = (sum(Psi_3(:).*x_solf(ii,:)')/n2)/(sum(Psi_3(:).*Psi_3(:))/n2);  
end
% Same description as given before ax(i) terms, except that these are
% coefficient vectors for solution in 'y-direction', and now we will
% consider y_solf, instead of x_solf.
ay0 = zeros(m*t,1); 
ay1 = zeros(m*t,1); 
ay2 = zeros(m*t,1); 
ay3 = zeros(m*t,1); 
% As explained above, here we finally calculate the coefficients in a loop
% as per the MC based PCE formula. Here 'n2' is the total no. of samples,
% and also the no. of columns of matrices y_solf
for jj = 1:m*t
  ay0(jj) = (sum(Psi_0(:).*y_solf(jj,:)')/n2)/(sum(Psi_0(:).*Psi_0(:))/n2); 
  ay1(jj) = (sum(Psi_1(:).*y_solf(jj,:)')/n2)/(sum(Psi_1(:).*Psi_1(:))/n2); 
  ay2(jj) = (sum(Psi_2(:).*y_solf(jj,:)')/n2)/(sum(Psi_2(:).*Psi_2(:))/n2); 
  ay3(jj) = (sum(Psi_3(:).*y_solf(jj,:)')/n2)/(sum(Psi_3(:).*Psi_3(:))/n2); 
end
% Calculating the final coefficients based on 'x' and 'y' projection
a0 = zeros(m*t,1); 
a1 = zeros(m*t,1); 
a2 = zeros(m*t,1); 
a3 = zeros(m*t,1);
for mm = 1:m*t
    a0(mm) = sqrt(ax0(mm).^2 + ay0(mm).^2);
    a1(mm) = sqrt(ax1(mm).^2 + ay1(mm).^2);
    a2(mm) = sqrt(ax2(mm).^2 + ay2(mm).^2);
    a3(mm) = sqrt(ax3(mm).^2 + ay3(mm).^2);
end
% Simple Check to see if everything is all right - Construct the solution
% matrix again from the coefficients obtained, and check the norm of the
% original solution and constructed solution
y_check = zeros(m*t,n2);
x_check = zeros(m*t,n2);
for kk = 1:m*t
  x_check(kk,:) = Psi_0*ax0(kk) + Psi_1*ax1(kk) + Psi_2*ax2(kk) + Psi_3*ax3(kk);
  y_check(kk,:) = Psi_0*ay0(kk) + Psi_1*ay1(kk) + Psi_2*ay2(kk) + Psi_3*ay3(kk);
end 
ynorm1 = norm(y_check(:,5));
ynorm2 = norm(y_solf(:,5));
xnorm1 = norm(x_check(:,5));
xnorm2 = norm(x_solf(:,5));
% Plot of 3rd coefficient in x direction for all time steps 
figure(1)
hold on;
for g = 1:m:m*t
    v = g+m-1;
    plot(ax3(g:v))
end   
hold off
% Exporting the resulting coefficients back to PARAVIEW for viewing
% Coefficient a0 for all 3 time steps 
Z0_1 = a0(1:m);
Z0_2 = a0(m+1:2*m);
Z0_3 = a0(2*m+1:3*m);
vtktrisurf(TRI,X1,Y1,Z0_1,'Coefficient_a0_1','data_a0_1.vtu')
vtktrisurf(TRI,X1,Y1,Z0_2,'Coefficient_a0_2','data_a0_2.vtu')
vtktrisurf(TRI,X1,Y1,Z0_3,'Coefficient_a0_3','data_a0_3.vtu')
% Coefficient a1 for all 3 time steps 
Z1_1 = a1(1:m);
Z1_2 = a1(m+1:2*m);
Z1_3 = a1(2*m+1:3*m);
vtktrisurf(TRI,X1,Y1,Z1_1,'Coefficient_a1_1','data_a1_1.vtu')
vtktrisurf(TRI,X1,Y1,Z1_2,'Coefficient_a1_2','data_a1_2.vtu')
vtktrisurf(TRI,X1,Y1,Z1_3,'Coefficient_a1_3','data_a1_3.vtu')
% Coefficient a2 for all 3 time steps 
Z2_1 = a2(1:m);
Z2_2 = a2(m+1:2*m);
Z2_3 = a2(2*m+1:3*m);
vtktrisurf(TRI,X1,Y1,Z2_1,'Coefficient_a2_1','data_a2_1.vtu')
vtktrisurf(TRI,X1,Y1,Z2_2,'Coefficient_a2_2','data_a2_2.vtu')
vtktrisurf(TRI,X1,Y1,Z2_3,'Coefficient_a2_3','data_a2_3.vtu')
% Coefficient a3 for all 3 time steps 
Z3_1 = a3(1:m);
Z3_2 = a3(m+1:2*m);
Z3_3 = a3(2*m+1:3*m);
vtktrisurf(TRI,X1,Y1,Z3_1,'Coefficient_a3_1','data_a3_1.vtu')
vtktrisurf(TRI,X1,Y1,Z3_2,'Coefficient_a3_2','data_a3_2.vtu')
vtktrisurf(TRI,X1,Y1,Z3_3,'Coefficient_a3_3','data_a3_3.vtu')
