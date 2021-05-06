function plot_edges(edges,no2xyz,plotNr,color)
    figure(plotNr), hold on
    for edge = edges
        cords = [no2xyz(:,edge(1)),no2xyz(:,edge(2))];
        plot3(cords(1,:),cords(2,:),cords(3,:),color);
        hold on
    end


end