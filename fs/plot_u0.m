% Load the data (space-separated)
data = load('u0.dat');

% Extract columns
x = data(:,1);
y = data(:,2);
z = data(:,3);
w = data(:,4);

% Example 1: 2D plot of first two columns
figure;
plot(x, y, '-o');
xlabel('x');
ylabel('y');
title('Plot of y vs x');
grid on;

% Example 2: 3D trajectory using first three columns
figure;
plot3(x, y, z, '-o');
xlabel('x');
ylabel('y');
zlabel('z');
title('3D plot of (x, y, z)');
grid on;

% Example 3: If you want to visualize all four columns as curves vs index
figure;
plot(data);
legend('col1','col2','col3','col4');
xlabel('Index');
ylabel('Value');
title('All columns vs index');
grid on;
