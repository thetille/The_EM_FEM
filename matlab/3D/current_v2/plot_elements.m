clear all

load cylinder_waveguide
figure
hold on

for i = 1:length(el2no)
    color = rand(1,3);
    nodes = no2xyz(:,el2no(:,i));
    line1 = [nodes(:,1),nodes(:,2)];
    line2 = [nodes(:,1),nodes(:,3)];
    line3 = [nodes(:,1),nodes(:,4)];
    line4 = [nodes(:,2),nodes(:,3)];
    line5 = [nodes(:,2),nodes(:,4)];
    line6 = [nodes(:,3),nodes(:,4)];
    
    
    plot3(line1(1,:),line1(2,:),line1(3,:),'color',color)
    plot3(line2(1,:),line2(2,:),line2(3,:),'color',color)
    plot3(line3(1,:),line3(2,:),line3(3,:),'color',color)
    plot3(line4(1,:),line4(2,:),line4(3,:),'color',color)
    plot3(line5(1,:),line5(2,:),line5(3,:),'color',color)
    plot3(line6(1,:),line6(2,:),line6(3,:),'color',color)
end