tetra = rand(3,4);
point = rand(3,1);


plot3(tetra(1,[1,2]),tetra(2,[1,2]),tetra(3,[1,2]),'k')
hold on
plot3(tetra(1,[1,3]),tetra(2,[1,3]),tetra(3,[1,3]),'k')
plot3(tetra(1,[1,4]),tetra(2,[1,4]),tetra(3,[1,4]),'k')
plot3(tetra(1,[2,3]),tetra(2,[2,3]),tetra(3,[2,3]),'k')
plot3(tetra(1,[2,4]),tetra(2,[2,4]),tetra(3,[2,4]),'k')
plot3(tetra(1,[3,4]),tetra(2,[3,4]),tetra(3,[3,4]),'k')


plot3(point(1),point(2),point(3),'x');


in_tetra = true;
i = 0;
while in_tetra || i == 4
    %normal := cross(v2 - v1, v3 - v1)
    normal := cross(v2 - v1, v3 - v1)
    
end