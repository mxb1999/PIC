
func makeplots(np){

	window,0;fma;
        var=phi
	p6;	// this is the color palette
	limits,-0.015,0.135,-0.075,0.075;
        plf,var,y_var,x_var;
        plf,var,-y_var,x_var;
        plfc,var,y_var(1:0,1:0),x_var(1:0,1:0),levs=span(min(var),max(var),50);
        plfc,var,-y_var(1:0,1:0),x_var(1:0,1:0),levs=span(min(var),max(var),50);
        plm,y_var,x_var
        colorbar,min(var),max(var);
        xytitles,"X (m)","Y (m)";
        pltitle,"Potential !F (V)";
        plmk,part_x(1:np,2),part_x(1:np,1),marker=4,width=11,msize=0.05


        window,1;fma;
        var=sqrt(efx^2+efy^2);
	p7;	// this is the color palette
	limits,-0.015,0.135,-0.075,0.075;
        plf,var,y_var,x_var;
        plf,var,-y_var,x_var;
        plfc,var,y_var(1:0,1:0),x_var(1:0,1:0),levs=span(min(var),max(var),50);
        plfc,var,-y_var(1:0,1:0),x_var(1:0,1:0),levs=span(min(var),max(var),50);
        plm,y_var,x_var
        colorbar,min(var),max(var);
        xytitles,"X (m)","Y (m)";
        pltitle,"|E| (V/m)";
        plmk,part_x(1:np,2),part_x(1:np,1),marker=4,width=11,msize=0.05

        window,2;fma;
        var=den;
	p9;	// this is the color palette
	limits,-0.015,0.135,-0.075,0.075;
        plf,var,y_var,x_var;
        plf,var,-y_var,x_var;
        plfc,var,y_var(1:0,1:0),x_var(1:0,1:0),levs=span(1.0e4,1.0e12,50);
        plfc,var,-y_var(1:0,1:0),x_var(1:0,1:0),levs=span(1.0e4,1.0e12,50);
        plm,y_var,x_var
        colorbar,min(var),1.0e12;
        xytitles,"X (m)","Y (m)";
        pltitle,"!r (#/m^3^)";

}
