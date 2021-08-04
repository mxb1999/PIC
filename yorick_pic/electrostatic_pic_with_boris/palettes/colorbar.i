func colorbar(cmin, cmax){
	plsys, 0;
	pli, span(0,1,200)(-,), .625, .46, .67, .84, legend="";
	plg, [.46,.84,.84,.46],[.67,.67,.625,.625], closed=1, marks=0, color="fg", width=1, type=1, legend="";
	plsys, 1;
	if (!is_void(cmin)){;
		c1=(cmax+cmin)/2.

		c2=(cmax+c1)/2.
		c3=(cmax+c2)/2.
		c4=(c1+c2)/2.

		c5=(cmin+c1)/2.
		c6=(cmin+c5)/2.
		c7=(c1+c5)/2.

		loc1=(0.84+0.46)/2.0

		loc2=(0.84+loc1)/2.0
		loc3=(0.84+loc2)/2.0
		loc4=(loc1+loc2)/2.0

		loc5=(0.46+loc1)/2.0
		loc6=(0.46+loc5)/2.0
		loc7=(loc1+loc5)/2.0

		plt, pr1(cmin), .6475, .46, justify="LH";
		plt, pr1(cmax), .6475, .84, justify="LH";
		plt, pr1(c1), 0.6475, loc1, justify="LH"
		plt, pr1(c2), 0.6475, loc2, justify="LH"
		plt, pr1(c3), 0.6475, loc3, justify="LH"
		plt, pr1(c4), 0.6475, loc4, justify="LH"
		plt, pr1(c5), 0.6475, loc5, justify="LH"
		plt, pr1(c6), 0.6475, loc6, justify="LH"
		plt, pr1(c7), 0.6475, loc7, justify="LH"
	};
};

