
function [total_nodes, nodes, total_elements, triangles, tetrahedron] = import_msh_v2(filename)
    file_id = fopen(filename,'r');

    while ~strcmp(fgetl(file_id),'$PhysicalNames')
    end
    total_groups = fscanf(file_id,'%d',1);
    for i = 1:total_groups
       group(i).dimention = fscanf(file_id, ['%d'],1);
       group(i).number = fscanf(file_id, ['%d'],1);
       group(i).name = fscanf(file_id, ['%s'],1);
    end
    % Nodes
    while ~strcmp(fgetl(file_id),'$Nodes')
    end
    total_nodes = fscanf(file_id,'%d',1);   

    for i = 1:total_nodes
       nodes(i,:)= fscanf(file_id, ['%*f', '%f', '%f', '%f'],3);
    end

    % Elements
    while ~strcmp(fgetl(file_id),'$Elements')
    end
    total_elements = fscanf(file_id,'%d',1);   


    temp = fgetl(file_id); %start on new line
    i = 0;
    line = fgetl(file_id);
    A = sscanf(line,'%d')';
    while length(A) == 8
        i = i+1; 
        triangles(i,:) = A;
        line = fgetl(file_id);
        A = sscanf(line,'%d')';
    end

    i = 0;
    while length(A) == 9
        i = i+1;
        tetrahedron(i,:) = A;
        line = fgetl(file_id);
        A = sscanf(line,'%d')';
    end

    fclose(file_id);
end