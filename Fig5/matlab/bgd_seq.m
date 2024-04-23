    lea  = [0 0; 15000 0;15000 1; 0 1];
    v  = [15000 0; 30000 0;30000 1; 15000 1];
    lea1  = [30000 0; 45000 0;45000 1; 30000 1];
    v2 = [45000 0; 60000 0;60000 1; 45000 1];
    lea2  = [60000 0; 75000 0;75000 1; 60000 1];
    v3 = [75000 0; 90000 0;90000 1; 75000 1];
    lea3  = [90000 0; 105000 0;105000 1; 90000 1];
    
    v4 = [105000 0; 120000 0;120000 1; 105000 1];
    lea4  = [120000 0; 135000 0;135000 1; 120000 1];
    v5 = [135000 0; 150000 0;150000 1; 135000 1];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices', lea, 'FaceColor', color_learning, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', lea1, 'FaceColor', color_learning, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', lea2, 'FaceColor', color_learning, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', lea3, 'FaceColor', color_learning, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', lea4, 'FaceColor', color_learning, 'EdgeColor', 'none')
    
    patch('Faces',f,'Vertices', v, 'FaceColor', color_pink, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', v2, 'FaceColor', color_pink, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', v3, 'FaceColor', color_pink, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', v4, 'FaceColor', color_pink, 'EdgeColor', 'none')
    patch('Faces',f,'Vertices', v5, 'FaceColor', color_pink, 'EdgeColor', 'none')
    