!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  This f90 code simulate a CF in a Perseus-like cluster
!!  This version includes the Fe enrichment and heating by stars (winds & SNIa)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! TO BE DONE: ADD FEEDBACK FROM AGN


MODULE DATA
integer :: N, i
real*8 :: pi,cmpc,cmkpc,yr,kbol,mu,mp,msun,grav
logical :: use_cooling
parameter (N=1001) 
parameter(pi=3.141592)
parameter(cmpc=3.085d18)
parameter(cmkpc=1000.*cmpc)
parameter(yr=3.156d7)
parameter(kbol=1.38d-16)
parameter(mu=0.61)
parameter(mp=1.67d-24)
parameter(msun=1.989d33)
parameter(grav=6.672e-8)
END MODULE DATA

PROGRAM ZEUS
USE DATA
!!IMPLICIT NONE
real*8 :: xa(N), xb(N), xmax, xmin, deltax, dxa(N), dxb(N)    !! more parameters than I actually use !!
real*8 :: d(N),e(N),v(N),P(N),s(N),tem(N),q(N),bri(N),ds,ent(N),tdyn(N)
real*8 :: g2a(N), g2b(N), g31a(N), g31b(N), dvl1a(N), dvl1b(N) 
real*8 :: F1(N),Flux_e(N),Flux_s(N),M(N),dstar(N),e_dstar(N),vstar(N),e_d(N)
real*8 :: divV(N),rshock(50),rsedov,rbubble,mcoldadd,mhotadd
real*8 :: Ecin,Eter,EterIN, mdot, mcold(N),mhot(N),mcoldold, mvir,&
          delta_time,timeold,mcoldtot,gasfrac,lx,lx500,nbv(N),dlogk
real*8 :: dtmin, tmax,t,C2,gam,Cv,k,t1,t2,t3,tcontR,tcontT,LumX,cfl,deltat
real*8 :: mstars,mnfw,mbh,mgal(N),mdark(N),ggal(N),gdark(N),gbh(N),conc, &
          fc,reff,ahern,rvir,rs,tem0,den0,mtot(N),gtot(N),lnd(N),yy, &
          factor,tcool(N),mgasin(N),mgas(N),vol(N),barfrac(N),deltat2 
real*8 :: dfe(N),zfe(N),dfestar(N),Ffe(N),tbv(N),twind,lwind,vwind, &
          mdotwind,temwind,rhost(N),zfest(N),zfesn,zfesol,alphast,alphasn,&
          alphatot,tst(N),tsn,t00,e00, ledd, M_dot(N)
integer :: sdr,Num,ncicli,isn,iagn,ncont,nsn,ivir,ncont2
character*7 nomefile,brifile,abundfile

real*8, EXTERNAL :: Cool       !! radiative cooling subroutine (see below) !!


!GENERATE THE TWO GRIDS xa AND xb

xmin=0.
xmax=3500.*cmkpc
deltax=50.*cmpc   !! first Delta_x
factor=1.0072     !! Delta_x(j+1) = factor*Delta_x(j)

!GRID a (the old value of r(j))
!do i=2,N   !! this is for a uniform grid !!
!   xa(i)= xmin+(xmax-xmin)*(i-2.)/(N-2.)
!end do
!xa(1)=-xa(3)

xa(1)=-deltax
xa(2)=0.
xa(3)=deltax
do i=4,N               !! use a non-uniform grid !!
   xa(i)=xa(i-1)+(xa(i-1)-xa(i-2))*factor
enddo

!GRIGLIA b (the old value of rr(j))
do i=1, N-1
	xb(i)=0.5*(xa(i)+xa(i+1))
end do
xb(N)=xb(N-1)+(xb(N-1)-xb(N-2))   !! add the last calculated Delta_xb to xb(N-1)

do i=2, N-1
	dxa(i)=xa(i+1)-xa(i)
	dxb(i)=xb(i)-xb(i-1)
end do

dxa(1)=xa(2)-xa(1) !boundary values
dxa(N)=dxa(N-1)
dxb(1)=dxb(2)
dxb(N)=xb(N)-xb(N-1)

open(20,file='grid.dat')  !! here write the grid on file !!
do i=1,N
   write(20,1001)xa(i)/cmkpc,xb(i)/cmkpc,dxa(i)/cmkpc,dxb(i)/cmkpc
enddo
close(20)
1001 format(4(1pe12.4))

!METRIC SCALE FACTORS AND CHOICE OF COORDINATES

sdr=1    !! use Cartesian coordinates, always !!

if (sdr==0) then  !! Cartesian - not used here !!

	do i=1, N
	!fattori metrici
	g2a(i)=1.
	g2b(i)=1.
	g31a(i)=1.
	g31b(i)=1.
 
	end do
	do i=1, N-1
	dvl1a(i)=xa(i+1)-xa(i)   !! !volume in cartesian coord centered in xb(i), x(b) is the old rr(j)
	end do
	dvl1a(N)=dvl1a(N-1)
	do i=2, N
	dvl1b(i)=xb(i)-xb(i-1)  !! volume in cartesian coord centered in xa(i), x(a) is the old r(j)
	end do
	dvl1b(1)=dvl1b(2)
	
else if (sdr==1) then   !! spherical coordinates
	do i=1, N
	g2a(i)=xa(i) !metrical parameter
	g31a(i)=xa(i)
	g2b(i)=xb(i)
	g31b(i)=xb(i)
	end do

	do i=1, N-1
	dvl1a(i)=(xa(i+1)**3-xa(i)**3)/3. !volume in spherical coord
	end do
	dvl1a(N)=dvl1a(N-1)
	do i=2, N
	dvl1b(i)=(xb(i)**3-xb(i-1)**3)/3.
	end do
	dvl1b(1)=dvl1b(2)

end if

! INITIAL CONDITIONS & PARAMETERS OF THE CLUSTER

 gam=5./3. !gam=gamma, gam is 5/3 bc the gas monoatomic
 cv=kbol/((gam-1.)*mu*mp)
 c2=3.            !! artificial viscosity constant
 cfl=0.01         !! initial small Delta t

 gasfrac=0.16296    !! cosmological baryon (not gas...) fraction !!
 mvir=1.4e15*msun   !! virial mass of the cluster (includes DM and gas)
 mnfw=(1.-gasfrac)*mvir   !! Total dark matter mass (remove Mgas from Mvir)
 rvir=2866.8*cmkpc  !! Virial radius
 rs=449.66*cmkpc    !! Scale radius for NFW halo (depends on mnfw!)
 conc=rvir/rs       !! concentration of NFW halo
 fc=log(1.+conc)-conc/(1.+conc) 

 mstars=1.e12*msun   !! Stellar mass of the BCG
 reff=10.*cmkpc      !! Effective radius
 ahern=reff/(1.+sqrt(2.)) !parametro di scala del profilo di Hernquist
 aml=7.5          !! stellas mass-to-light ratio
 slope=1.1        !! exponent for the power law SNIa rate
 snu= 0.3         !! a fiducial current SNIa rate - it can be changed
 sigma=2.5e7    !! 1D stellar velocity dispersion, constant for now 
 sigma2=sigma*sigma
 tst=mu*mp*sigma2/kbol  !! tst (equivalent stellar temperature) is a vector
 tsn=1.79e9       !! equivalent T for SNIa ejecta: 3/2 kT/(mu*mp) = 1e51/1.4 Msun
 zfesol=1.8e-3    !! solar Fe abundance (by mass)
 temmin=1.e4      !! minimum T allowed in the simulation 
 
 rhost=mstars/(2.*pi)*(ahern/xb)/(ahern+xb)**3  !! stellar density (Hernquist)
 if(xb(i).ge.100.*cmkpc) rhost(i)=0.   !ge is >= so stars truncated at 100 kpc
 zfest=0.8*zfesol   !! stellar Fe abundance, constant for now
 zfesn=0.744/1.4    !! this is y_Ia/M_ejecta
 
 open(15,file='rhost.dat')
 do i=2,N
    write(15,1009)xb(i)/cmkpc,rhost(i)
 enddo
 close(15)
 1009 format(2(1pe12.4))

 mbh=1.e9*msun       !! Central BH mass
 ledd=1.3d46*(mbh/(1.e8*msun))   !! Eddington luminosity - use it later (maybe)
 tnow=13.7           !! Current time in Gyr

!!  Find the grid point corresponding to rvir (it might be useful)

   do i=1,N
      if(xb(i).gt.rvir)then !!gt means >
         ivir=i-1
         goto11
      endif
   enddo
11 continue
   print*,'ivir=',ivir,real(xa(i)/cmkpc),real(rvir/cmkpc)

!! CALCULATE MASSES (AND GRAVITY)

 do i=3,N      !! notice: do loop starts with i=3
    mgal(i)=mstars*xa(i)**2/(xa(i)+ahern)**2  !! note: galactic mass 
    ggal(i)=grav*mgal(i)/xa(i)**2             !! galactic gravity
    yy=xa(i)/rs
    mdark(i)=mnfw/fc*(log(1.+yy)-yy/(1.+yy))  !! NFW halo mass
    gdark(i)=grav*mdark(i)/xa(i)**2           !! NFW gravity
    gbh(i)=grav*mbh/xa(i)**2                  !! BH gravity
 enddo
 mgal(1)=0.
 mgal(2)=0.
 ggal(1)=0.
 ggal(2)=0.
 mdark(1)=0.
 mdark(2)=0.
 gdark(1)=0.
 gdark(2)=0.
 gbh(1)=0.
 gbh(2)=0.
!g=gravity 
do i=1,N
   gtot(i)=ggal(i)+gdark(i)+gbh(i) !! total gravity includes galaxy, DM & BH
   mtot(i)=mgal(i)+mdark(i)+mbh
enddo

open(24,file='grv.dat')     !! write some interesting stuff
open(25,file='mass.dat')
do i=3,N
   write(24,1005)xa(i)/cmkpc,ggal(i),gdark(i),gbh(i),gtot(i)
   write(25,1005)xa(i)/cmkpc,mgal(i)/msun,mdark(i)/msun,mbh/msun,mtot(i)/msun
enddo
close(24)
close(25)
1005 format(5(1pe12.4))


!! Now set the initial conditions: ICM in hydrostatic equilibrium
!! The initial temperature is a fit of Perseus cluster.

    temp0=8.12e7  !! parameters for the initial T profile
    rtemp1=71.
    rtemp2=71.

!! Set the initial ICM temperature

 do i=2,N
    rkpc=xb(i)/cmkpc
    tem(i)=temp0*( &
           (1.+(rkpc/rtemp1)**3)/(2.3+(rkpc/rtemp2)**3) )
 enddo
 tem(1)=tem(2)

 den0=0.87*2.17d-24 !! Central initial gas density (for deltax=50 pc)

lnd(2)=log(den0)
lnd(1)=log(den0)
do i=3,N        !! Here solves the hydrostatic equilibrium equation
   lnd(i)=lnd(i-1) - (xb(i)-xb(i-1))*gtot(i)*   &
          (mu*mp)/(kbol*0.5*(tem(i)+tem(i-1)))  &
           - log(tem(i)) + log(tem(i-1))
enddo

!! SET THE INITIAL CONDITIONS FOR THE ICM

do i=1, N
   v(i)=0.
   d(i)=exp(lnd(i)) !density
   e(i)=cv*d(i)*tem(i) !energy
   p(i)=e(i)*(gam-1.) !pressure
   zfe(i)=0.3*zfesol  !! set the initial abundance to 0.3 solar (early SNII!)   
   dfe(i)=zfe(i)*d(i)/1.4  !! iron density !!
end do

!!!!!!!!!! Here write the initial conditions and calculate some stuff !!!!!!!!!!!

!! Calculate the initial gas mass and the baryon faction
open(20,file='mgasin.dat')
do i=2,N-1
   vol(i)=1.3333333*pi*(xa(i+1)**3-xa(i)**3)
enddo
vol(N)=vol(N-1)
mgasin(1)=0. !initial mass of gas
do i=2,N-1
   mgasin(i)=mgasin(i-1)+d(i)*vol(i)
   barfrac(i)=(mgasin(i)+mgal(i))/(mgasin(i)+mgal(i)+mdark(i))
enddo
 do i=3,N-1
   write(20,1006)xa(i)/cmkpc,mgasin(i)/msun
 enddo
1006 format(5(1pe12.4))

print*,'f_bar(rvir) = ',real(barfrac(ivir)) !! this is the baryon fraction at virial radius

 ent=kbol*tem/(d/1.937d-24)**0.66666667   !! this is the initial entropy profile

 do i=2, N-1
    tcool(i)=max(e(i)/Cool(tem(i),d(i)),3.156e7) !! cooling time: key quantity   !Cool=radiative cooling
    gmed=0.5*(gtot(i)+gtot(i+1)) !gtot media
    tdyn(i)=sqrt(2*xb(i)/gmed)  !! this is the dynamical time
 end do
 !boundary condition
 tcool(1)=tcool(2) 
 tcool(N)=tcool(N-1)
 tdyn(1)=tdyn(2) 
 tdyn(N)=tdyn(N-1)

 !! Here calculate the initial X-ray surface brightness
 
 open(20,file='briiniz.dat')   !briiniz= initial brightness ! in this file writes also tcool
 do i=2,ivir+1  !!N-1
    bri(i)=0.
    do nn=i,N-2
       epsilonx=0.
       if(tem(nn).ge.1.e6)epsilonx=Cool(tem(nn),d(nn)) 
       !tem(nn)= temperatura della cella nn
       !episilonx=emissività volumentrica-->variazione di luminosità in termini di volume, calcolata solo se T>10^6K-->it represents the cooling function
       if(nn.eq.i)then
          ds=sqrt(xa(nn+1)**2-xb(i)**2) !spessore della shell, xa(nn) e xb(nn)= coord radiali(limiti delle celle)
       else
          ds=sqrt(xa(nn+1)**2-xb(i)**2)- &
             sqrt(xa(nn)**2-xb(i)**2)
       endif
       bri(i)=bri(i)+epsilonx*ds 
    enddo
    bri(i)=max(1.d-25,2.*bri(i))   !! "2.*" in order to count the other side
    write(20,1008)xb(i)/cmkpc,bri(i),ent(i),&
                  tdyn(i)/yr,tcool(i)/yr !rr, brightness, entropy, tdyn, tcool
                  !tempo/yr= normalizzazione del tempo
 enddo
 close(20)

!! write the initial conditions in files resiniz.dat & abundiniz.dat
open(20,file='resiniz.dat')
  do i=3,N
     write(20,1000)xa(i)/cmkpc,xb(i)/cmkpc,d(i)/1.937d-24,v(i)/1.e5,&
          tem(i),p(i),tcool(i)/yr
  enddo
close(20)

open(20,file='abundiniz.dat')
  do i=2,N
     write(20,1010)xb(i)/cmkpc,dfe(i),zfe(i)/zfesol
  enddo
close(20)
1010 format(3(1pe12.4))

!!!! Define some more parameters and start the time integration !!!!

 t=0.        !!
 timeold=0.    !! this is used to calculate the cooling rate, see below
 tmax=5.e9*yr  !! final time of the simulation 
 tshift=8.7*(1.e9*yr)  !! this is for the stellar source term
 ncicli=0
 ncont=0
 ncont2=0
 deltat=tmax/5.    !! the code writes results every "deltat" (every Gyr in this case)
 deltat2=tmax/500. !! a file with some other stuff is written every "deltat2"

 nomefile='clu0000norad' !file del cluster
 brifile='bri0000norad' !file del brightness
 abundfile='abu0000norad' !file di abundance

 open(22,file='mass_time.dat')  !! writes the time evolution of masses
 open(23,file='mdot_time.dat')  !! writes the time evolution of cooling rate

!!!!!!!! HERE STARTS THE TIME INTEGRATION OF THE HYDRO EQUATIONS !!!!!!!
 
do while (t<tmax)
        ncicli=ncicli+1     !! this is just a counter of cycles
!!        if(ncicli.gt.1) goto 1111

        
!! CALCULATE DTMIN

        dtmin=1.d30   !! an arbitrary very large value
        !dtmin=passo temporale minimo consentito per rispettare condizione di CFL
        p=(gam-1.)*e !pressione, e=energia interna, gam=gamma
	do i=2, N-1
	!ridurre il passo minimo consentito=dtmin
		 dtmin=min(dtmin,(xb(i)-xb(i-1))/(abs(v(i))+sqrt(gam*p(i)/d(i))))
		 !cs=sound speed=sqrt(gam*p(i)/d(i))
		 !formula CFL: il passo temporale non può superare il tempo che l'informazione (in questo caso fluido+onde sonore) impiega ad attraversare la cella-->stabilià: il valore minimo preso su tutte le celle 
	end do
       !cfl factor is a stability factor 
       cfl=min(0.5,1.1*cfl)   !! here increase cfl, up to 0.5
        dtmin=cfl*dtmin
        t=t+dtmin
!!        print*,'ncicli, dtmin = ',ncicli, real(dtmin),real(t)


!! SOURCE STEP

!! SUBSTEP 0: STELLAR TERMS, IRON DENSITY, SN ENERGY...

   tgyr=(t+tshift)/(1.e9*yr)   !! tshift = 8.7*1.e9*yr
   ! to simulate a stellar pop of initial age "tshift"
   !alpha coeff for stellar wind
   alphast=4.7e-20*(tgyr/tnow)**(-1.26)   !! new alphast for this cycle
!!   alphast=0.         !! if you don't want stellar winds

!alpha coeff for supernova
   alphasn=4.436e-20*(snu/aml)*(tgyr/tnow)**(-slope)
!!   alphasn=0.         !! if you don't want SNIa

   alphatot=alphasn+alphast !tasso totale di perdita di massa

 do i=2,N-1
    d(i)=d(i)+dtmin*alphatot*rhost(i)  !! add stellar mass loss !d density of gas
    dfe(i)=dfe(i)+dtmin*rhost(i)*(alphast*zfest(i)/1.4+alphasn*zfesn) !! add Fe density!!
    t00=(alphast*tst(i)+alphasn*tsn)/(alphatot+1.d-30) !temperatura media del gas iniettato da parte di SNe e staalar wind in ICM
    e00=cv*t00 
!!         if(i.eq.5)print*,'e00 = ',real(e00),real(t00/1.e7)
    e(i)=e(i)+dtmin*alphatot*rhost(i)*(e00+0.5*v(i)**2) !energ interna
 enddo

!!print*,'alp ',real(alphast), real(alphasn)

 CALL BCb(d)
 CALL BCB(dfe)
 CALL BCb(e)

!! SUBSTEP I: UPDATE OF VELOCITY FOR PRESSURE GRADIENT pag22a zeus

	do i=2, N-1
		v(i)=v(i)-dtmin*2.*(P(i)-P(i-1))/((d(i)+d(i-1))*dxb(i)) &
                     - gtot(i)*dtmin	
	end do
	CALL BCa(v) !set the boundary conditions


!! CALCULATE ARTIFICIAL VISCOSITY Q pag 23a zeus 
	do i=2, N-1
		if ((v(i+1)-v(i))<0.) then !do not consider to verify the equilibrium
			q(i)=C2*d(i)*(v(i+1)-v(i))**2
		else 
			q(i)=0.
		end if
	end do

	CALL BCb(q) !set the boundary conditions

!! SUBSTEP II: UPDATE FOR THE ARTIFICIAL VISCOSITY pag 22b zeus

	do i=2, N-1
		v(i)=v(i)-dtmin*2.*(q(i)-q(i-1))/((d(i)+d(i-1))*dxb(i))
	end do

	CALL BCa(v) !set the boundary conditions

	do i=2, N-1
		e(i)=e(i)-dtmin*q(i)*(v(i+1)-v(i))/dxa(i) !viscous dissipation
	end do
	CALL BCb(e) !set the boundary conditions

!! SUBSTEP III: UPDATE FOR COMPRESSIONAL HEATING pag 23b zeus
	do i=2,N-1
		divV(i)=(g2a(i+1)*g31a(i+1)*v(i+1)-g2a(i)*g31a(i)*v(i))/dvl1a(i)
	!divV=divergence of V
	end do
	CALL BCa(divV) !set the boundary conditions

	do i=2, N-1
		e(i)=e(i)*(1.-0.5*dtmin*(gam-1.)*divV(i))/(1.+0.5*dtmin*(gam-1.)*divV(i))
	end do
	CALL BCb(e) !set the boundary conditions

!!  Here update T when it is needed


!    RADIATIVE COOLING
!    -- function Cool calculates the cooling function for T>1e4 K,
!       using a fit of the tabulated Sutherland & Dopita results.
!       Also note that we are NOT considering the cooling time as
!       a constraint for dtmin !!
use_cooling = .false.
        tem=e/(cv*d)   !! update T after the compressional source step

        do i=2, N-1
        	if (use_cooling) then
                 e(i)=e(i)-dtmin*Cool(tem(i),d(i))   !! here Cool=n_e*n_i*Lambda
                 tcool(i)=max(e(i)/Cool(tem(i),d(i)),3.156e7)
       		else 
                 tcool(i)=1.0d20  
        	end if ! tempo di raffreddamento infinito
        end do
        tcool(1)=tcool(2)
        tcool(N)=tcool(N-1)

        call BCb(e) !set the boundary conditions

        tem=e/(cv*d)

        do i=2, N-1
           if(tem(i)<1.d4) then
              tem(i)=1.d4
              e(i)=tem(i)*Cv*d(i)
           end if
        end do
        call BCb(tem) !set the boundary conditions
        call BCb(e) !set the boundary conditions
        P=(gam-1.)*e


!! TRANSPORT STEP

	do i=2, N-1
		s(i)=0.5*(d(i)+d(i-1))*v(i)  !! this is at "i" !!
		!s=rho*v is momentum
	end do	

	CALL BCa(s)

!! UPDATE OF THE DENSITY, IRON DENSITY AND ENERGY
 
     !! CALCULATE MASS AND IRON MASS FLUXES
        do i=2, N-1
		if (v(i)>0.) then
			dstar(i)=d(i-1)     !dtallar winds density ! at i !!
			dfestar(i)=dfe(i-1)  !!dfestar=stellar iron density  ! at i !!
			

		else
			dstar(i)=d(i)
			dfestar(i)=dfe(i)
		end if
	end do
	dstar(N)=dstar(N-1)
	dstar(1)=dstar(3)
	dfestar(N)=dfestar(N-1)
	dfestar(1)=dfestar(3)

	do i=2, N
		
		F1(i)=dstar(i)*v(i)*g2a(i)*g31a(i)  !gas density flux  due to stellar winds p25b zeus  !! at i !!	
		Ffe(i)=dfestar(i)*v(i)*g2a(i)*g31a(i)    !!desity ironflux  ! at i !!	
	end do

    !! CALCULATE ENERGY FLUX

	do i=2, N-1
		M_dot(i)=dstar(i)*v(i)   !! gas mass flux for consistent advection !!
	end do
	CALL BCa(M)
	
	
	do i=2, N-1
		if (v(i)>0.) then
			e_dstar(i)=e(i-1)/d(i-1)   !energy=energy/density ! at i !!
		else
			e_dstar(i)=e(i)/d(i)
		end if
	end do
	e_dstar(N)=e_dstar(N-1)
	e_dstar(1)=e_dstar(3)


     !! NOW FINALLY UPDATE DENSITY, IRON DENSITY AND ENERGY
	do i=2, N-1
		d(i)=d(i)-dtmin*(F1(i+1)-F1(i))/dvl1a(i) !densita calcolata dopo avere ricavato flusso di densità formula pag 25b zeus
		dfe(i)=dfe(i)-dtmin*(Ffe(i+1)-Ffe(i))/dvl1a(i) !new iron density
	end do 
	CALL BCb(d)
	CALL BCb(dfe)
        zfe=1.4*dfe/d       !new iron abundance ! in absolute units

	do i=2, N
		Flux_e(i)=e_dstar(i)*M_dot(i)*g2a(i)*g31a(i)	!energy gas flux			
	end do
	CALL BCa(Flux_e)

	do i=2, N-1
		e(i)=e(i)-dtmin*(Flux_e(i+1)-Flux_e(i))/dvl1a(i) !new energy of the gas 
	end do

	CALL BCb(e)
        do i=1,N
           tem(i)=e(i)/(cv*d(i)) !new temperature of the gas
        enddo


!! UPDATE OF THE MOMENTUM

	do i=2, N-1
		if ((v(i-1)+v(i))*0.5>0) then
			vstar(i)=v(i-1)       !steallar winds velocity ! at i-1/2
		else
			vstar(i)=v(i)
		end if
	end do

	CALL BCb (vstar)

	do i=1, N-1
	Flux_s(i)=vstar(i+1)*0.5*(M_dot(i)+M_dot(i+1))*g2b(i)*g31b(i)   !momentum flux ! at i+1/2
	end do
	
	do i=2, N-1
		s(i)=s(i)-dtmin/dvl1b(i)*(Flux_s(i)-Flux_s(i-1)) !new momentum
	end do

	CALL BCa(s)

	do i=2, N-1
		v(i)=2.*s(i)/(d(i)+d(i-1)) !gas velocity 
	end do

	CALL BCa(v)

!!!!  Here the hydro update ends  !!!!
       !!! NOW WRITE STUFF !!!

 if(t.ge.(ncont+1)*deltat)then  !! write results every time interval deltat
    ncont=ncont+1
    call incnam(nomefile)   !! update name of the result datafile "clu00norad**"
    open(20,file=nomefile)   !r,rr,densità,velocità,temperatura, pressione, tcool

    call incnam(brifile)    !! update name of the brightness datafile "bri00**"
    open(42,file=brifile) !

    call incnam(abundfile)  !! update name of the abund datafile "abu00**
    open(21,file=abundfile) !new iron abundance

    do i=3,N     !! starts from 3 !!
       write(20,1000)xa(i)/cmkpc,xb(i)/cmkpc,d(i)/1.937d-24,v(i)/1.e5, &
                     tem(i),p(i),tcool(i)/yr
    end do
    1000 format(7(1pe11.3))
    close(20)

    do i=2,N
        write(21,1010)xb(i)/cmkpc,dfe(i),zfe(i)/zfesol
    enddo
    close(21)

!! Calculate the X-ray surface brightness and other stuff

 ent=kbol*tem/(d/1.937d-24)**0.66666667   !! this is the "entropy"
!d/1.937d-24=per convertire la densità in numero di particelle

!brightness profile integrando lungo line-of-sight
 do i=2,ivir+1   !! X-ray brightness within rvir
    bri(i)=0.
    do nn=i,N-2
       epsilonx=0.
       if(tem(nn).ge.1.e6)epsilonx=Cool(tem(nn),d(nn)) !cool=cooling function
       if(nn.eq.i)then
          ds=sqrt(xa(nn+1)**2-xb(i)**2) !spessore della shell
       else
          ds=sqrt(xa(nn+1)**2-xb(i)**2)- & 
             sqrt(xa(nn)**2-xb(i)**2) !xa e xb sono gli estremi delle spehrical shell
       endif
       bri(i)=bri(i)+epsilonx*ds
    enddo
    bri(i)=max(1.d-25,2.*bri(i))   !! factor of 2 to count for the other side
    write(42,1008)xb(i)/cmkpc,bri(i),ent(i),&
                  tdyn(i)/yr,tcool(i)/yr !rr, brightness, entropy, tdyn, tcool
 enddo
1008 format(5(1pe12.4))
 close(42)

 endif

 if(t.ge.(ncont2+1)*deltat2)then  !! calculate and print some stuff every "deltat2"
    ncont2=ncont2+1

!! Calculate the X-ray luminosity

 lx=0. !xray luminosity
 lx500=0. !xray luminosity entro raggio R500
 do i=2, N-1   !! within rvir or R500 and for gas with X-ray temperature
    if(xa(i+1).le.rvir.and.tem(i).ge.1.e6)then
       lx=lx+Cool(tem(i),d(i))*vol(i)
    endif
    if(xa(i+1).le.0.5*rvir.and.tem(i).ge.1.e6)then  !! for r < R500
       lx500=lx500+Cool(tem(i),d(i))*vol(i)
    endif
 end do


!! Calculate the masses of hot gass, cold gas, extra cold mass

    delta_time=t-timeold !! these two command lines are used to calculate the cooling rate
    mcoldold=mcold(ivir)

    tcold=1.e5  !! this T defines the "cold gas"
    mcold(1)=0.
    mhot(1)=0.
    do i=2,N-1
       mcoldadd=0. !extra cold mass 
       mhotadd=0.
       if(tem(i).le.tcold) mcoldadd=d(i)*vol(i)
       mcold(i)=mcold(i-1)+mcoldadd
       if(tem(i).gt.tcold) mhotadd=d(i)*vol(i) !gt is >
       mhot(i)=mhot(i-1)+mhotadd
    enddo

    write(22,1100)t/(1.e9*yr),mhot(ivir)/msun,mcold(ivir)/msun,&
          lx/1.d45,lx500/1.d45 !stampati in mass_time.dat

1100 format(5(1pe12.4))

!    print*,'mcold(ivir) = ',real(mcold(ivir)/msun)
!    print*,'mcoldold = ',real(mcoldold/msun)
!    print*,'timeold,t,delta_time = ',timeold/3.156e16,t/3.156e16,&
!            delta_time/3.156e16
    mdot=(mcold(ivir)-mcoldold)/delta_time   !! this is the cooling rate
    timeold=t

    write(23,1101)t/(1.e9*yr),mdot/msun*yr
1101 format(2(1pe12.4))

 endif                   !! end calculating & printing

enddo       !!!!!! here the "do while" ends !!!!!
1111 continue

1002 format(3(1pe12.4))

close(22)
close(23)


END PROGRAM ZEUS


SUBROUTINE BCa(z1) !! reflection BC for velocity and momentum
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z1

z1(2)=0.
z1(1)=-z1(3)
z1(N)=z1(N-1)
!z1(1)=z1(2)       !! ouflow !!
!z1(N)=z1(N-1)

END SUBROUTINE BCa

SUBROUTINE BCb(z2)  !! outflow BC for density, energy...
USE DATA
IMPLICIT NONE
real*8, dimension (N) :: z2
z2(1)=z2(2) !boundary condition
z2(N)=z2(N-1) !boundary condition
END SUBROUTINE BCb


Real*8 FUNCTION Cool(tem1, d1)
USE DATA
IMPLICIT NONE
Real*8:: tem1, d1, qlamsh,temkev,kev


real*8 lambda, radfac, lambda1,qlambda0,qlambda1,met
      real*8 sigmat,maccr
      real tlog,qlog0,arg,bump1rhs,bump2lhs,tlogc,qlogc,&
           qloginfty,ppp,qlog1,tlogm,qlogm,sig

!       kev=1.16e7    !! this is a old cooling funcion - don't use !!
!       temkev=tem1/kev
!       if(tem1.ge.2.32e5)then
!          qlamsh=1.d-22*(8.6d-3*temkev**(-1.7)+0.058*temkev**0.5 &
!                        +0.063)
!       elseif(tem1.ge.2.d4.and.tem1.lt.2.32d5)then
!          qlamsh=6.72d-22*(temkev/0.02)**0.6
!       elseif(tem1.lt.2.d4.and.tem1.ge.1.d4)then
!          qlamsh=1.544d-22*(temkev/0.0017235)**6
!       else
!          qlamsh=0.
!       endif
!       Cool=0.6831*(d1/1.67d-24)**2*qlamsh

!!  Here is a new (better) fit to Lambda(T,Z):

      lambda=0.
      tlogc = 5.65 !temp critica in log_10
      qlogc = -21.566 !valore log_10 di funzione di raffreddamento metallica alla temp critica
      qloginfty = -23.1 !valore limite a cui la curva metallica può scendere alle alte temper
      ppp = .8 !pendenza serve a modellare il calo della emissività dei metalli ad alte temp
      tlogm = 5.1 !temp caratteristica
      qlogm = -20.85 !valore max della funzione di raffreddamento metallica al picco
      sig = .65 !sigma di gaussiana-->larghezza del bump metallico-->controlla quanto è ampio il picco attorno a tlogm

      tlog=log10(tem1)
      if(tlog .ge. 6.1) then
         qlog0 = -26.39 + .471*(log10(tem1 + 3.1623e6)) !emissività per H/He
      else if(tlog .ge. 4.9) then
         arg = 10.**(-(tlog-4.9)/.5) + .077302 !variabile di appoggio
         qlog0 = -22.16 + log10(arg)
      else if(tlog .ge. 4.25) then
         bump1rhs = -21.98 - ((tlog-4.25)/.55)
         bump2lhs = -22.16 - ((tlog-4.9)/.284)**2
         !bump1rhs e bump2lhs sono due curve di raccordo (right-hand side e left-hand side) usate per modellare la transizione tra regimi di raffreddamento a basse temp
         qlog0 = max(bump1rhs,bump2lhs)
      else
         qlog0 = -21.98 - ((tlog-4.25)/.2)**2
      endif
      qlambda0 = 10.**qlog0

!     emission from metals alone at solar abundance:
      if(tlog .ge. 5.65) then
         qlog1 = qlogc -ppp*(tlog - tlogc)
         qlog1 = max(qlog1,qloginfty) !emissività da metalli
      else
         qlog1 = qlogm - ((tlog - tlogm)/sig)**2
      endif
         qlambda1 = 10.**qlog1

!     final cooling coefficient:
      lambda = qlambda0 + 0.4*qlambda1   !! this is Lambda(T), for Z=0.4 solar
      Cool=lambda*0.6831*(d1/1.67d-24)**2


END FUNCTION Cool


! ================================================================
      subroutine incnam (name)
! ================================================================
! --------------------------------------------
! --- incnam increments the last 4 digits in name by +1.
!     name is type character*7 (right-justified).
! --------------------------------------------
      character*7 name

!!!      print*,'prima: ',name
      if(name(7:7) .eq. '9') go to 20 !se ultimo carattere non è 9-->incremento di 1 usando ichar e char, se è 9 passare al blocco successivo
!!!        print*,'prima name7: ',name(7:7)
      name(7:7) = char(1 + ichar(name(7:7)))
!!!        print*,'dopo name7: ',name(7:7)
      return

20    if(name(6:6) .eq. '9') go to 30
      name(6:6) = char(1 + ichar(name(6:6)))
      name(7:7) = '0'
      return
30    if(name(5:5) .eq. '9') stop
      name(5:5) = char(1 + ichar(name(5:5)))
      name(6:6) = '0'
      name(7:7) = '0'
!!!      print*,'dopo: ',name

      return
      end

!to sum up: il programma è un contatore numerico con riporto
!!Incrementa l’ultima cifra se non è '9'.
!!Se è '9', fa il “riporto” sulla cifra precedente.
!!Se anche quella è '9', fa riporto ancora più indietro.
!!Se tutte le cifre sono '9', il programma si ferma.
