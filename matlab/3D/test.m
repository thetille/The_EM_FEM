clc
clear
clf(1)
clf(2)
%close all
% physical triangle
%Ph_triangle = rand(3,2);
Ph_triangle = [0.5, 0.1; 0.5, 0.5; 0.1, 0.5];
Ref_triangle = [0, 0; 1, 0; 0, 1];



for i = 1:3
    x = Ph_triangle(1,1) + (Ph_triangle(2,1)-Ph_triangle(1,1))*Ref_triangle(i,1)+(Ph_triangle(3,1)-Ph_triangle(1,1))*Ref_triangle(i,2);
    y = Ph_triangle(1,2) + (Ph_triangle(2,2)-Ph_triangle(1,2))*Ref_triangle(i,1)+(Ph_triangle(3,2)-Ph_triangle(1,2))*Ref_triangle(i,2);
    sol_triangle(i,:) =  [x,y]
end

%grid
figure(1)
hold on
[x_grid, y_grid] =  meshgrid((0:0.05:1),(0:0.05:1));
u_inc = sin((pi*x_grid/1));

quiver(x_grid,y_grid,zeros(length(x_grid),length(y_grid)), u_inc);


plot([sol_triangle(1,1), sol_triangle(2,1)],[sol_triangle(1,2), sol_triangle(2,2)],'r')
plot([sol_triangle(2,1), sol_triangle(3,1)],[sol_triangle(2,2), sol_triangle(3,2)],'g')
plot([sol_triangle(3,1), sol_triangle(1,1)],[sol_triangle(3,2), sol_triangle(1,2)],'b')

figure(2)
hold on
[u_grid, v_grid] =  meshgrid((0:0.05:1),(0:0.05:1));
trix = Ph_triangle(:,1);
triy = Ph_triangle(:,2);
u_inc_u = (triy(1)+(triy(2)-triy(1))*v_grid)* sin((pi*(trix(1)+(trix(2)-trix(1))*u_grid+(trix(3)-trix(1))*v_grid )/1));
u_inc_v = (triy(1)+(triy(3)-triy(1))*u_grid)* sin((pi*(trix(1)+(trix(2)-trix(1))*u_grid+(trix(3)-trix(1))*v_grid )/1));

quiver(u_grid,v_grid,u_inc_v,u_inc_u);

plot([Ref_triangle(1,1), Ref_triangle(2,1)],[Ref_triangle(1,2), Ref_triangle(2,2)],'r')
plot([Ref_triangle(2,1), Ref_triangle(3,1)],[Ref_triangle(2,2), Ref_triangle(3,2)],'g')
plot([Ref_triangle(3,1), Ref_triangle(1,1)],[Ref_triangle(3,2), Ref_triangle(1,2)],'b')
