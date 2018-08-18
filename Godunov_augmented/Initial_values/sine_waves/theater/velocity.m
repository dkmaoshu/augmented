% 	program velocity

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,3);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\show\exact\exact_solution.dat
    xx=exact_solution(:,1);
    yy=exact_solution(:,3);
    plot(xx,yy,'-');
    hold off
            
    axis([-0.01 1.01 -0.5 1.5]);
   