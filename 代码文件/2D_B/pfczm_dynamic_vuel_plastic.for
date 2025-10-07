
!**********************************************************************************************************
!
      module NumKind
!
!**********************************************************************************************************
        implicit none
        integer (kind(1)), parameter :: ikind = kind(1), 
     &                                  rkind = kind(0.D0), 
     &                                  lkind = kind(.true.)
!       
      end module Numkind
      

!**********************************************************************************************************
!
      module ModelParam
!
!**********************************************************************************************************
        use NumKind
        implicit none

        ! Constants

        ! Flag of initilization
        logical (lkind) :: bInitialized = .false.

        ! Tolerance
        real    (rkind), parameter :: TOL = 1.0d-12 
        ! number of guass points
        integer (ikind), parameter :: ngp = 2
        
        ! geometric function parameter
        real(rkind), parameter :: c0 = 3.1415926535897932384626433832d0
		
		! user defined variables
		integer (ikind), parameter :: NodeNum = 95116, NumEle=94500
        real(rkind), save :: allP(NumEle)
        
        ! 
        real(rkind) :: thk, EA, nu, Gf, ft, sigma_y, lb, rho
		real(rkind) :: eta, time_step, h_min, c_0, c_1
		real(rkind) :: c, S_1, gama_0, alpha_0, dec_type
		real(rkind) :: lamda, mu
        real(rkind) :: De(3, 3)        
        real(rkind) :: p, a1, a2, a3
        real(rkind) :: gp(ngp), gw(ngp)
        real(rkind) :: QQ(12,12) 

        !
        contains

        !===================================
          subroutine Initialize(props, nprops, istype)
        !===================================

            integer (ikind), intent (in) :: nprops, istype
            real    (rkind), intent (in) :: props(nprops)

            !********************************************
            real(rkind) :: G0, K11, K12
            integer(ikind) :: indexq(12), i

            ! material properties
            EA         =  props(1)  ! props(1) -- Young's modulus
            nu         =  props(2)  ! props(2) -- Poisson's ratio
            ft         =  props(3)  ! props(3) -- failure strength
			sigma_y    =  props(4)  ! props(4) -- yield stress
            Gf         =  props(5)  ! props(5) -- fracture energy
            lb         =  props(6)  ! props(6) -- length scale            
            thk        =  props(7)  ! props(7) -- thickness			
			rho        =  props(8)  ! props(8) -- density 
			eta        =  props(9)  ! props(9) -- viscous coffient 
			time_step  =  props(10) ! props(10) -- time step
			h_min      =  props(11) ! props(11) -- minimun mesh size
			c_0        =  props(12) ! props(12) -- c0 for artificial bulk viscosity method
			c_1        =  props(13) ! props(13) -- c1 for artificial bulk viscosity method
			c          =  props(14) ! props(14) -- c for Grüneisen equation
			S_1        =  props(15) ! props(15) -- S1 for Grüneisen equation
			gama_0     =  props(16) ! props(16) -- gama0 for Grüneisen equation
			alpha_0    =  props(17) ! props(17) -- alpha for Grüneisen equation
			dec_type   =  props(18) ! props(18) -- type of the decomposition method
			
            if (thk < TOL) thk = 1.0
            
            ! elastic stiffness matrix
			lamda   = EA * nu / (1.d0 + nu) / (1.d0 - 2.d0*nu)
			mu      = EA / (2.d0 * (1.d0 + nu))
            G0      = EA / (2.d0 * (1.d0 + nu))
            K11     = EA * (1.d0 - nu) / (1.d0 + nu) / (1.d0 - 2.d0*nu)
            K12     = EA * nu / (1.d0 + nu) / (1.d0 - 2.d0*nu)
            De(:,1) = (/ K11,  K12, 0.D0/)
            De(:,2) = (/ K12,  K11, 0.D0/)
            De(:,3) = (/0.D0, 0.D0,   G0/)
            
            ! softening parameters
            a1   =  4.d0/(c0*lb)*EA*Gf/(ft*ft)
            if      (istype == 1) then  ! linear softening
              p  =  2.d0
              a2 = -0.5d0
              a3 =  0.0d0
            else if (istype == 2) then  ! exponential softening
              p  =  2.5d0
              a2 =  2.0d0**(5.d0/3.d0) - 3.0d0
              a3 =  0.0d0
            else if (istype == 3) then  ! blinear softening
              p  =  2.0d0
              a2 =  0.03687d0
              a3 =  20.8343d0
            else if (istype == 4) then  ! concrete softening
              p  =  2.0d0
              a2 =  1.3868d0
              a3 =  0.6567d0
            else if (istype == 5) then  ! hyperbolic softening
              p  =  4.0d0 
              a2 =  2.0d0**(7.d0/3.d0) - 4.5d0
              a3 =  0.0d0
            else
              write (*,*) '**error: Softening law No. ', istype, 
     &                    'does not exist!'
            end if
            
            ! integration points
            gp = (/ -1.d0, 1.d0 /) / dsqrt(3.d0)
            gw = (/  1.d0, 1.d0 /)
            
            ! dof interchange
            indexq = (/ 1,2,9, 3,4,10, 5,6,11, 7,8,12 /)
            ! interchange the locations of dofs
            QQ = 0.d0
            do i = 1, 12
              QQ(indexq(i),i) = 1.d0
            end do             
            
            bInitialized = .true.
            
            return
          end subroutine Initialize
      !========================================================================= 
      end module ModelParam

!**********************************************************************************************************
!
      module FEM
!
!**********************************************************************************************************
        use NumKind
        implicit none

        contains      
          !==================shape function and its derivative with xi and eta======================    
          subroutine shapefuc(n, dn_xieta, xi, eta)
          
            implicit none      
            real(rkind) :: n(4), dn_xieta(2, 4), xi, eta

            n(1) = 0.25d0*(1.d0 - xi)*(1.d0 - eta)
            n(2) = 0.25d0*(1.d0 + xi)*(1.d0 - eta)
            n(3) = 0.25d0*(1.d0 + xi)*(1.d0 + eta)
            n(4) = 0.25d0*(1.d0 - xi)*(1.d0 + eta)
            
            dn_xieta(1, 1) = -0.25d0*(1.d0 - eta)
            dn_xieta(1, 2) =  0.25d0*(1.d0 - eta)
            dn_xieta(1, 3) =  0.25d0*(1.d0 + eta)
            dn_xieta(1, 4) = -0.25d0*(1.d0 + eta)
            
            dn_xieta(2, 1) = -0.25d0*(1.d0 - xi)
            dn_xieta(2, 2) = -0.25d0*(1.d0 + xi)
            dn_xieta(2, 3) =  0.25d0*(1.d0 + xi)
            dn_xieta(2, 4) =  0.25d0*(1.d0 - xi)
            
            return 
          end subroutine shapefuc

          !===============traditional b matrix==============================================      
          subroutine b_matrix(nd,bd,nn,b,det_jacb, coords,xi,eta)
          
            implicit none
            real(rkind) :: nd(4), bd(2,4), nn(2,8), b(3,8)
            real(rkind) :: jacb(2,2), inv_jacb(2,2), coords(2, 4)
            real(rkind) :: det_jacb, xi, eta
            
            !local varibles
            real(rkind) :: n(4), dn_xieta(2,4), dn_x(4), dn_y(4)
            integer(ikind) :: i, j
             
            ! shape functions 
            call shapefuc(n,dn_xieta,xi,eta)
            nd = n
			nn = 0.d0
			do i = 1,4
			  nn(1,2*i-1) = n(i)
			  nn(2,2*i) = n(i)
			end do
            
            ! jacob matrix
            jacb = matmul(dn_xieta, transpose(coords))            
            det_jacb = jacb(1,1)*jacb(2,2) - jacb(1,2)*jacb(2,1)
            inv_jacb(1, 1) = jacb(2, 2)
            inv_jacb(1, 2) =-jacb(1, 2)
            inv_jacb(2, 1) =-jacb(2, 1)
            inv_jacb(2, 2) = jacb(1, 1)
            inv_jacb = 1.d0/det_jacb*inv_jacb            
            
            !initialize varibles
            do i = 1,4
              dn_x(i) = inv_jacb(1,1)*dn_xieta(1,i)
     &                + inv_jacb(1,2)*dn_xieta(2,i)
              dn_y(i) = inv_jacb(2,1)*dn_xieta(1,i)
     &                + inv_jacb(2,2)*dn_xieta(2,i)
            end do
			
            ! B matrix for displacement
            b = 0.d0
            do j = 1, 4
              b(1, 2*(j-1) + 1) = dn_x(j)
              b(2, 2*(j-1) + 2) = dn_y(j)
              b(3, 2*(j-1) + 1) = dn_y(j)
              b(3, 2*(j-1) + 2) = dn_x(j)
            end do
            
            ! B matrix for damage
            do j = 1,4
              bd(1,j) = dn_x(j)
              bd(2,j) = dn_y(j)
            end do
          
            return
          end subroutine b_matrix
      
        !********************************************************************
        ! define the dyadic function
          function dyadic(vector1,vector2, vlen)
        !********************************************************************
            integer (ikind) :: vlen, i, j
            real    (rkind) :: vector1(vlen),vector2(vlen)
            real    (rkind) :: dyadic(vlen,vlen)
          
            do i = 1, vlen
              do j = 1, vlen
                dyadic(i,j) = vector1(i) * vector2(j)
              end do
            end do

            return
          end function dyadic

      end module FEM
	  
!**********************************************************************************************************
!
      subroutine computemass(amass,coords)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        use FEM
        implicit none

        real(rkind):: mass_u(8,8), mass_d(4,4), mm(12,12), amass(12,12), coords(2,4)     
        ! local varibles
        real(rkind):: nn(2,8), b(3,8), nd(4), bd(2,4)
        real(rkind):: det_jacb, dvol       
        integer(ikind):: i, j  
		
		mass_u  = 0.d0
		mass_d  = 0.d0
        do i = 1, ngp
          do j = 1, ngp      
            call b_matrix(nd,bd,nn,b,det_jacb, coords,gp(i),gp(j))
			dvol=  gw(i)*gw(j)*det_jacb*thk
			mass_u =  mass_u + dvol*rho*matmul(transpose(nn), nn)	
            mass_d =  mass_d + dvol*eta*dyadic(nd, nd, 4)		
          end do
        end do
		
		!lumped mass matrix
		mm = 0.d0
		do i=1,8
          do j=1,8
            mm(i,i) = mm(i,i) + mass_u(i,j)
          enddo
        enddo
		do i=1,4
          do j=1,4
            mm(i+8,i+8) = mm(i+8,i+8) + mass_d(i,j)
          enddo
        enddo
		
		!remap matrix
		amass = matmul(matmul(transpose(QQ),mm),QQ)
        return 
      end subroutine computemass
	  
!**********************************************************************************************************
!
      subroutine computerhs(rhs,coords,u,svars,steptime,elemNum)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        use FEM
        implicit none

        real(rkind):: rhs(12), coords(2,4)
        real(rkind):: svars(40), u(12)
        integer(ikind):: elemNum	
        ! local varibles
        real(rkind):: nn(2,8), b(3,8), nd(4), bd(2,4), mass_d(4,4), dd_mass(4,4)
        real(rkind):: uu(8), dd(4), rd(4), ru(8), d_0(4), d_1(4)
        real(rkind):: rr(12), r1(4), r0(4)
        real(rkind):: stressEff(4), pre_stressEff(4), stress(3), s_stress(4)
		real(rkind):: strain(3), pre_strain(3), delta_strain(3)
        real(rkind):: det_jacb, pre_energy_crk, energy_crk, phi, omega, domega, ddomega
        real(rkind):: dalpha, ddalpha, phi_source, dvol
        integer(ikind):: i, j, k, l
        real(rkind):: delta_epsilonv, v_epsilonv, steptime   
        real(rkind):: interEnergy, pre_interEnergy, stress_p, pre_stress_p, deg_stress_p
		real(rkind):: energyEff, pre_energyEff, mod_q, Hp
        
		! extrat nodal displacement and damage dofs
        do i = 1, 4
          uu(2*i - 1) = u(3*i - 2)
          uu(2*i    ) = u(3*i - 1)
          dd(i)       = u(3*i)
        end do
		
        ! initialize varibles
		mass_d  = 0.d0
        rd  = 0.d0
		ru  = 0.d0
		allP(elemNum) = 0.d0
        do i = 1, ngp
          do j = 1, ngp
            k = (i - 1) * 2 + j		  
            call b_matrix(nd,bd,nn,b,det_jacb, coords,gp(i),gp(j))              
            strain = matmul(b, uu)  ! strain field
			
			! Return mapping algorithm (get the deviatoric stress)
			pre_stressEff = 0.d0 
			do l = 1, 3
			  pre_strain(l) = svars(3*(k-1)+4+l)    ! effective strain for the last step
			end do
			do l = 1, 4
			  pre_stressEff(l) = svars(4*(k-1)+16+l)  ! effective stress for the last step
			end do
			delta_strain = strain-pre_strain         ! increment of total strain
			call RMA(pre_stressEff,delta_strain,s_stress)
			
			! Grüneisen equation of state (get the hydrostatic stress)
			pre_stress_p = svars(k+32)
			pre_interEnergy = svars(k+36)
			delta_epsilonv = delta_strain(1)+delta_strain(2)
			call EOS(strain,delta_strain,s_stress,stress_p,interEnergy,
     &               pre_stressEff,pre_stress_p,pre_interEnergy,delta_epsilonv)
			! Or through linear elasticity (get the hydrostatic stress)
			! call linearElasticP(pre_stressEff,delta_strain,stress_p)
			
			! update stressEff
			stressEff(1) = s_stress(1)-stress_p
			stressEff(2) = s_stress(2)-stress_p
			stressEff(3) = s_stress(3)-stress_p
			stressEff(4) = s_stress(4)			
						
			! crack driving force
			pre_energy_crk = svars(k)
			pre_energyEff = svars(k+40)
			if(dec_type == 1.d0) then
		      call decMaxPrin(stressEff,energyEff)
		    else if(dec_type == 2.d0)then
		      call decVolEnergy(stress_p,delta_epsilonv,energyEff,pre_energyEff)
			else
              write (*,*) '**error: Decomposition method No. ', dec_type, 
     &                    'does not exist!'
			end if
            energy_crk = max(max(energyEff,pre_energy_crk),0.5d0*ft**2.d0/EA)
			
			! limit the phase field
			phi  = dot_product(nd,dd) ! crack phase-field
			if (phi>1.d0) then
              phi=1.d0
            else if (phi<0.d0) then
              phi=0.d0
            endif
			
			! artificial bulk viscosity method
			v_epsilonv = (delta_strain(1)+delta_strain(2))/time_step
			call modifystressEff(v_epsilonv,mod_q,steptime)
			
			! degradated stress
			call geometricFunc(dalpha,ddalpha,phi) ! geometric function
            call energeticFunc(omega,domega,ddomega,phi) ! energetic function 	
			call degradateStress(omega,stress_p,s_stress,mod_q,stress)
			
			! residual for damage
            phi_source  = domega *energy_crk + Gf/(c0*lb)*dalpha
            dvol =  gw(i)*gw(j)*det_jacb*thk
			mass_d =  mass_d + dvol*eta*dyadic(nd, nd, 4)
            rd  =  rd  + dvol*(phi_source*nd + 2.d0*lb*Gf/c0
     &          *  matmul(transpose(bd), matmul(bd, dd)))	        
            ru  =  ru + dvol*matmul(transpose(b),stress)
			
			! update svars
			svars(k) = energy_crk
			do l = 1, 3
			  svars(3*(k-1)+4+l) = strain(l)
			end do
			do l = 1, 4
			  svars(4*(k-1)+16+l) = stressEff(l)
			end do
			svars(k+32) = stress_p
			svars(k+36) = interEnergy
	        svars(k+40) = energyEff
			
			! updata userdefined variables
			if (stress_p > 0) then
		      Hp = 1.d0
		    else
		      Hp = 0.d0
		    end if
			deg_stress_p = omega*(1.d0-Hp)*stress_p + Hp*stress_p
			allP(elemNum) = allP(elemNum) + deg_stress_p
          end do
        end do
		
		allP(elemNum) = allP(elemNum) / 4.d0
		
		! constraint of phase field
		do i=1,4
          do j=1,4
              dd_mass(i,i) = dd_mass(i,i)+mass_d(i,j)
          enddo
		  d_0(i) = 0.d0-dd(i)
		  d_1(i) = 1.d0-dd(i)
        enddo
		r1 = matmul(dd_mass,d_1)/time_step
		r0 = matmul(dd_mass,d_0)/time_step
	    do i=1,4
          if (-rd(i) > r1(i)) then
              rd(i) = -r1(i)
          endif
		  if (-rd(i) < r0(i)) then
              rd(i) = -r0(i)
          endif
        enddo

		! remap
		rr = 0.d0
        rr(1:8 ) = ru
        rr(9:12) = rd
		rhs = matmul(transpose(QQ),rr)
        return 
      end subroutine computerhs
      
!**********************************************************************************************************
!
      subroutine energeticFunc(omega,domega,ddomega,phi)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        implicit none
      
        real(rkind) :: omega, domega, ddomega, phi
        real(rkind) :: fac1, dfac1, ddfac1, fac2, dfac2, ddfac2
      
        fac1    =  (1.d0 - phi)**p
        dfac1   = -p*(1.d0 - phi)**(p - 1.d0); 
        ddfac1  =  p*(p - 1.d0)*(1.d0 - phi)**(p - 2.d0)
        
        fac2    =  fac1   + a1*phi + a1*a2*phi**2.d0 + a1*a2*a3*phi**3.d0
        dfac2   =  dfac1  + a1 + 2.d0*a1*a2*phi + 3.d0*a1*a2*a3*phi**2.d0
        ddfac2  =  ddfac1 + 2.d0*a1*a2 + 6.d0*a1*a2*a3*phi
        
        omega   =  fac1/fac2 + 1.0d-6        
        domega  =  (dfac1*fac2  - fac1*dfac2)/(fac2**2.d0)
        ddomega = ((ddfac1*fac2 - fac1*ddfac2)*fac2 - 2.d0*
     &             (dfac1*fac2 - fac1*dfac2)*dfac2)/(fac2**3.d0)
     
        return
      end subroutine energeticFunc
      
!**********************************************************************************************************
!
      subroutine geometricFunc(dalpha,ddalpha,phi)
!
!**********************************************************************************************************
        use NumKind
        implicit none
        
        real(rkind) :: dalpha, phi, ddalpha
        
        dalpha  = 2.d0 - 2.d0*phi
        ddalpha =-2.d0
        
        return 
      end subroutine geometricFunc  

!**********************************************************************************************************
!
      subroutine decMaxPrin(stressEff,energyEff)
!
!**********************************************************************************************************	  
	    use NumKind
		use ModelParam
        implicit none
		
		real(rkind):: savg, sdif, sdev, smax, smin, energyEff
		real(rkind):: stressEff(3)
		
	    savg = 0.5*(stressEff(1) + stressEff(2))
        sdif = 0.5*(stressEff(1) - stressEff(2))
        sdev = sqrt(sdif*sdif + stressEff(3)*stressEff(3))
        smax = savg + sdev
        smin = savg - sdev
		energyEff = 0.5d0*smax**2.d0/EA	  
        return 
      end subroutine decMaxPrin
	  
!**********************************************************************************************************
!
      subroutine decVolEnergy(stress_p,delta_epsilonv,energyEff,pre_energyEff)
!
!**********************************************************************************************************	  
	    use NumKind
		use ModelParam
        implicit none
		
		real(rkind):: stress_p, delta_epsilonv, energyEff, pre_energyEff
		real(rkind):: hp, delta_energyEff
		
		if(stress_p > 0.d0)then
		  hp = 1.d0
		else
		  hp = 0.d0
		end if
		delta_energyEff = -(1.d0-hp)*stress_p*delta_epsilonv
		energyEff = pre_energyEff + delta_energyEff
        return 
      end subroutine decVolEnergy	  
	  
	  
!**********************************************************************************************************
!
      subroutine RMA(pre_stressEff,delta_strain,s_stress)
!
!**********************************************************************************************************	  
	    use NumKind
		use ModelParam
        implicit none
		
		real(rkind):: pre_stressEff(4), delta_strain(3), s_stress(4)
		real(rkind):: sigmaEq_trial, p_trial, traceInc
		real(rkind):: sigma_trial(4), s_trial(4)
		integer(ikind):: i
		
		! plane strain condition
		traceInc = delta_strain(1) + delta_strain(2)
		sigma_trial(1) = pre_stressEff(1) + 2.d0*mu*delta_strain(1) + lamda*traceInc
	    sigma_trial(2) = pre_stressEff(2) + 2.d0*mu*delta_strain(2) + lamda*traceInc
        sigma_trial(3) = pre_stressEff(3) + lamda*traceInc		
		sigma_trial(4) = pre_stressEff(4) + 2.d0*mu*delta_strain(3)		
		p_trial = -(sigma_trial(1)+sigma_trial(2)+sigma_trial(3))/3.d0
		s_trial(1) = sigma_trial(1) + p_trial
	    s_trial(2) = sigma_trial(2) + p_trial
		s_trial(3) = sigma_trial(3) + p_trial
		s_trial(4) = sigma_trial(4)
		sigmaEq_trial = sqrt((3.0/2.0)*(s_trial(1)*s_trial(1)+s_trial(2)*s_trial(2)+s_trial(3)*s_trial(3))
     &            +3.d0*(s_trial(4)*s_trial(4)))
	    
		! check yield
		if (sigmaEq_trial <= sigma_y) then
		  s_stress = s_trial
		else
		  do i = 1, 4		    
			s_stress(i) = s_trial(i) - s_trial(i)/sigmaEq_trial*(sigmaEq_trial-sigma_y)			
		  end do
		end if
	     
        return 
      end subroutine RMA	  

!**********************************************************************************************************
!
      subroutine degradateStress(omega,stress_p,s_stress,mod_q,stress)
!
!**********************************************************************************************************	  
	    use NumKind
		use ModelParam
        implicit none
		
		real(rkind):: s_stress(4), stress(3)
		real(rkind):: omega, stress_p, mod_q
		integer(ikind):: i
		real(rkind):: Hp, mod_stress_p
				
		if (stress_p > 0) then
		  Hp = 1.d0
		else
		  Hp = 0.d0
		end if
		
		mod_stress_p = stress_p + mod_q
		stress(1) = omega*(s_stress(1)-(1.d0-Hp)*mod_stress_p) - Hp*mod_stress_p
	    stress(2) = omega*(s_stress(2)-(1.d0-Hp)*mod_stress_p) - Hp*mod_stress_p
		stress(3) = omega*s_stress(4)
	     
        return 
      end subroutine degradateStress	

!**********************************************************************************************************
      ! artificial bulk viscosity method
      subroutine modifystressEff(v_epsilonv,mod_q,steptime)
!
!**********************************************************************************************************	  
	    use NumKind
		use ModelParam
        implicit none
		
		real(rkind):: v_epsilonv, mod_q, C_b, steptime
		integer(ikind):: i

		C_b = sqrt(EA/(3.d0*rho*(1.d0-2.d0*nu)))
		if (v_epsilonv < 0.d0) then
		  mod_q = rho*h_min*(c_0*h_min*v_epsilonv**2.d0-c_1*C_b*v_epsilonv)
	    else
		  mod_q = 0.d0
		end if
		
		if (steptime == 0.d0) then
		  mod_q = 0.d0
		end if
		
        return 
      end subroutine modifystressEff
!**********************************************************************************************************
!
      subroutine EOS(strain,delta_strain,s_stress,stress_p,interEnergy,
     &               pre_stressEff,pre_stress_p,pre_interEnergy,delta_epsilonv)
!
!**********************************************************************************************************	  
	    use NumKind
		use ModelParam
        implicit none
		
		real(rkind):: stress_p, pre_stress_p, pre_interEnergy, delta_epsilonv, interEnergy
		real(rkind):: strain(3), delta_strain(3), s_stress(4), pre_stressEff(4)
		real(rkind):: pre_s_stress(4), delta_e(3) 
		integer(ikind):: i
        real(rkind):: A11, A12, A21, A22, b1, b2, delta
		real(rkind):: epsilonv, kapa, beta, Hk, delta_s_energy
		
		epsilonv = strain(1)+strain(2)
		kapa = -epsilonv/(1.d0+epsilonv)
		if (kapa > 0) then
		  Hk = 1.d0
		else
		  Hk = 0.d0
		end if
		beta = (1.d0-Hk)+Hk*(1.d0+(1.d0-0.5d0*gama_0)*kapa-0.5d0*alpha_0*kapa**2.d0)
     &          /(1.d0-(S_1-1.d0)*kapa)**2.d0
		pre_s_stress(1) = pre_stressEff(1)+pre_stress_p
		pre_s_stress(2) = pre_stressEff(2)+pre_stress_p
		pre_s_stress(3) = pre_stressEff(3)+pre_stress_p
		pre_s_stress(4) = pre_stressEff(4)
		delta_e(1) = delta_strain(1)-delta_epsilonv/3.d0
		delta_e(2) = delta_strain(2)-delta_epsilonv/3.d0
		delta_e(3) = delta_strain(3)
		delta_s_energy = (s_stress(1)+pre_s_stress(1))*delta_e(1) + (s_stress(2)+pre_s_stress(2))*delta_e(2)
     &          + 2.d0*(s_stress(4)+pre_s_stress(4))*delta_e(3)		
		
		A11 = 1.d0
		A12 = -(gama_0+alpha_0*kapa)
		A21 = delta_epsilonv/2.d0
		A22 = 1.d0
		b1 = beta*rho*c**2.d0*kapa
		b2 = pre_interEnergy - pre_stress_p/2.d0*delta_epsilonv + delta_s_energy/2.d0
		
		delta = A11*A22-A12*A21
		stress_p = (b1*A22-b2*A12)/delta
		interEnergy = (A11*b2-A21*b1)/delta
		
        return 
      end subroutine EOS
	  
!**********************************************************************************************************

      subroutine linearElasticP(pre_stressEff,delta_strain,stress_p)
!
!**********************************************************************************************************	  
	    use NumKind
		use ModelParam
        implicit none
		
		real(rkind):: pre_stressEff(4), delta_strain(3)
		real(rkind):: p_trial, traceInc, stress_p
		real(rkind):: sigma_trial(4)
		integer(ikind):: i
		
		! plane strain condition
		traceInc = delta_strain(1) + delta_strain(2)
		sigma_trial(1) = pre_stressEff(1) + 2.d0*mu*delta_strain(1) + lamda*traceInc
	    sigma_trial(2) = pre_stressEff(2) + 2.d0*mu*delta_strain(2) + lamda*traceInc
        sigma_trial(3) = pre_stressEff(3) + lamda*traceInc		
		sigma_trial(4) = pre_stressEff(4) + 2.d0*mu*delta_strain(3)		
		p_trial = -(sigma_trial(1)+sigma_trial(2)+sigma_trial(3))/3.d0
	    stress_p = p_trial
	    
        return 
      end subroutine linearElasticP	  

!**********************************************************************************************************
!

      subroutine VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars,
     &                energy,nnode,ndofel,props,nprops,jprops,
     &                njprops,coords,mcrd,u,du,v,a,jtype,jElem,
     &                time,period,dtimeCur,dtimePrev,kstep,kinc,
     &                lflags,dMassScaleFactor,predef,npredef,
     &                ndload, adlmag)
!**********************************************************************************************************

        use NumKind
        use ModelParam
        use FEM
        implicit none

       ! operational code keys
        integer (ikind), parameter :: jMassCalc = 1, jIntForceAndDtStable = 2, jExternForce = 3
       ! flag indices
        integer (ikind), parameter :: iProcedure = 1, iNlgeom = 2, iOpCode = 3, nFlags = 3
       ! energy array indices
        integer (ikind), parameter :: iElPd = 1, iElCd = 2, iElIe = 3, iElTs = 4, iElDd = 5, iElBv = 6,
     &            iElDe = 7, iElHe = 8, iUnused = 9, iElTh = 10, iElDmd = 11, iElDc = 12, nElEnergy = 12
       ! predefined variables indices
        integer (ikind), parameter :: iPredValueNew = 1, iPredValueOld = 2, nPred = 2     
       ! time indices
        integer (ikind), parameter :: iStepTime  = 1, iTotalTime = 2, nTime      = 2

!**********************************************************************************************************
      ! interface of vuel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! variables passed in
        integer (ikind), intent (in    ) :: nblock, nsvars, nnode, ndofel,
     &    nprops, njprops, mcrd, jtype, kstep, kinc, npredef, ndload 	 
        real (rkind), intent (in    ) :: dtimeCur, period, dtimePrev	 
        integer (ikind), intent (in    ) :: jElem(nblock),  lflags(*),
     &    jprops(njprops)
        real    (rkind), intent (in    ) :: props(nprops), coords(nblock,nnode,mcrd),
     &    u(nblock,ndofel), du(nblock,ndofel), v(nblock,ndofel), a(nblock,ndofel),
     &    time(nTime), dMassScaleFactor(nblock), adlmag(nblock),
     &    predef(nblock, nnode, npredef, nPred)
  
        ! variables to be updated (the update of energy(8) is optional)
        real    (rkind), intent (in out) :: dtimeStable(nblock), rhs(nblock,ndofel), 
     &    amass(nblock,ndofel,ndofel), svars(nblock,nsvars), energy(nblock,nElEnergy)
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !********************************************************************************************************               
      ! user coding to define rhs, amatrx, svars, energy and pnewdt (optional for the last two)
        ! initialize parameters, etc. 
		! meaning of svars: 1-4 H of integration points; 5-16 previous strain of integration points; 
		! 17-32 previous effective stress of integration points; 33-36 previous effective hydrostatic stress of integration points;
		! 37-40 previous internal energy of integration points; 41-44 previous effective driven energy of integration points;
		integer (ikind) :: kblock
        if (.not. bInitialized) call Initialize(props, nprops, jtype)
		if (lflags(iOpCode).eq.jMassCalc ) then
          do kblock = 1, nblock
            call computemass(amass(kblock,:,:),transpose(coords(kblock,:,:)))
          end do	  
        else if (lflags(iOpCode) .eq. jIntForceAndDtStable) then
          do kblock = 1, nblock
            call computerhs(rhs(kblock,:),transpose(coords(kblock,:,:)),u(kblock,:),svars(kblock,:),time(iStepTime),jElem(kblock))
            dtimeStable(kblock) = time_step
          end do	
        end if
        return
      end subroutine vuel
	  
!**********************************************************************************************************
      subroutine vusdfld(nblock,nstatev,nfieldv,nprops,ndir,nshr, 
     &                   jElemUid,kIntPt,kLayer,kSecPt,stepTime,
     &                   totalTime,dt,cmname,coordMp,direct,T, 
     &                   charLength,props,stateOld,stateNew,field)

        use NumKind
        use ModelParam
        use FEM
        implicit none
	    ! variables passed in
	    integer (ikind), intent (in    ) :: nblock, nstatev, nfieldv, nprops, 
     &	          ndir, nshr, kIntPt, kLayer, kSecPt
        real (rkind), intent (in    ) :: stepTime, totalTime, dt
		real (rkind), intent (in    ) :: coordMp(nblock,*), direct(nblock,3,3),
     &            T(nblock,3,3), charLength(nblock), props(nprops), stateOld(nblock,nstatev)
	    integer (ikind), intent (in    ) :: jElemUid(nblock)
        character(len=8) :: cmname
        ! variables to be updated
		real (rkind), intent (in out) :: stateNew(nblock,nstatev), field(nblock, nfieldv)
        integer (ikind) :: k
        
		do k = 1, nblock  
		  stateNew(k,1) = allP(jElemUid(k)-NumEle)
        end do

        return
      end subroutine	  
