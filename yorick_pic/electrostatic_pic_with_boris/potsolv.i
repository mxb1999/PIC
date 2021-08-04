/* Potential solver for 2D PIC program
        Based on the Gauss-Seidel method */

func eval_2dpot_GS(phi){
extern b0, b, x
tol = solvertol               // solver tolerance

// Get nx and ny from size of density
nx = dimsof(den)(2)
ny = dimsof(den)(3)
nn = numberof(den)      // these are not necessary but do it anyway

// Convert density and potential into column vectors
b0 = reform(den, numberof(den),1)
x = reform(phi, numberof(phi),1)

// Solve
for (it=1; it<=20000; ++it){
        // recalculate RHS
        b = b0 - n0*exp((x-phi0)/Te)    // add boltzmann term for electrons
        b = -b*QE/EPS0;

        // set boundaries
        b(indgen(1:nx)) = 0             // zero electric field on y = 0
        b(indgen(nn-nx+1:nn)) = 0       // zero electric field on y = L
        b(indgen(nx:nn:nx)) = 0         // zero electric field on x = L
        b(indgen(1:nn:nx)) = phi0       // fixed potential on x=0

        // set potential on fixed nodes
        for (j=box(2,1); j<=box(2,2); ++j){
                // wall potential
                b((box(1,1)+(j-1)*nx):(box(1,2)+(j-1)*nx)) = phi_p*array(1.0,box(1,2)-box(1,1)+1,1)
        }

        // update nodes
        for (i=1; i<=nn-1; ++i){
//                x(i) = (b(i) - A(i,1:i-1)*x(1:i-1) - A(i,i+1:nn)*x(i+1:nn))/A(i,i)
                newA1 = A(i,1:i-1)
                newA2 = A(i,i+1:nn)
                newx1 = x(1:i-1)
                newx2 = x(i+1:nn)
                x(i) = (b(i) - newA1(+)*newx1(+) - newA2(+)*newx2(+))/A(i,i)
        }

        for (i=nn; i<=nn; ++i){
//                x(i) = (b(i) - A(i,1:i-1)*x(1:i-1) - A(i,i+1:nn)*x(i+1:nn))/A(i,i)
                newA1 = A(i,1:i-1)
                newx1 = x(1:i-1)
                x(i) = (b(i) - newA1(+)*newx1(+))/A(i,i)
        }

        // compute residue to check for convergence, do only every 5 iterations
        if ( it%5 == 0){
                R = sqrt(sum((b(,1)-A(,+)*x(+))^2.0))-phi0		// residue
                RR = sqrt(sum((b(,1)-A(,+)*x(+))^2.0))/sqrt(sum(b(,1)^2.0))	// relative residue

                if ( RR <= tol ){
                        print, "  GS converged in",it,"iterations with norm",R,"and rel. R",RR
                        break;
                }
        }

	if ( it%200 == 0){
		if ( RR > tol ){
			print, "  ... iteration:",it,"R is",R,"and RR is",RR
		}
	}
}


// Check if the solver converged to the specified tolerance
if ( RR > tol ){
	print, "  GS failed to converge!! it is",it-1,"and RR is ",RR;
}

// return solution as nx*ny array
x=reform(x,nx,ny);
return x;

}
