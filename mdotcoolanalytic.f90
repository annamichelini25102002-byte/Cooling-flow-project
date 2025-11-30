program mdot_vs_time
  implicit none
  integer, parameter :: N=1000
  real*8 :: xa(N), xb(N), d(N), v(N), tem(N), p(N), tcool(N)
  real*8 :: mu, mp, kb, lx, mdot, Tmean, t_age, rcool, vol
  integer :: i, i_rcool, itime
  character*7 :: nomefile

  ! Costanti fisiche
  mu = 0.6d0
  mp = 1.67d-24
  kb = 1.38d-16

  ! Inizializza il nome del primo file
  nomefile = 'clu0000'

  ! Loop sugli snapshot (clu0001 ... clu0005)
  do itime=1,5 !itime=’indice del ciclo (1,2,3,4,5)Gyr
     call incnam(nomefile)   ! aggiorna il nome: clu0001, clu0002, ...
     open(10, file=nomefile, status='old')
     do i=1,N
        read(10,*) xa(i), xb(i), d(i), v(i), tem(i), p(i), tcool(i)
     end do
     close(10)

     ! Tempo in anni (1 Gyr, 2 Gyr, ...)
     t_age = itime * 1.d9

     ! Trova rcool dove tcool ~ t_age
     i_rcool = 1
     do i=1,N
        if (abs(tcool(i)-t_age) < abs(tcool(i_rcool)-t_age)) i_rcool = i !nei dati discreti non c'è mai uguaglianza perfetta!! per questo motivo cercare l’indice i che minimizza la differenza
     end do
     rcool = xa(i_rcool)

     ! Calcola Lx entro rcool usando direttamente Cool
     lx = 0.d0
     Tmean = 0.d0
     do i=2,i_rcool
        if(tem(i).ge.1.d6) then
           vol = (4.d0/3.d0) * 3.141592653589793d0 * (xb(i)**3 - xa(i)**3)
           lx = lx + Cool(tem(i), d(i)) * vol
           Tmean = Tmean + tem(i) * Cool(tem(i), d(i)) * vol
        endif
     end do
     Tmean = Tmean / lx   ! media pesata sulla luminosità

     ! Formula analitica
     mdot = 2.d0 * mu * mp * lx / (5.d0 * kb * Tmean)

     print *, 'File:', nomefile, 'rcool=', rcool, 'Mdot_analytic=', mdot
  end do

contains

  ! Funzione di raffreddamento: Cool = ne * ni * Lambda(T)
  real*8 function Cool(T, rho)
    implicit none
    real*8, intent(in) :: T, rho
    real*8 :: ne, ni, Lambda

    ! Stima ne ~ ni ~ rho/(mu*mp)
    ne = rho/(mu*mp)
    ni = rho/(mu*mp)

    !  funzione Lambda(T)
    Lambda = 1.d-23 * sqrt(T)

    Cool = ne * ni * Lambda
  end function Cool
  
    subroutine incnam(name)
    character*7 name
    if(name(7:7) .eq. '9') go to 20
    name(7:7) = char(1 + ichar(name(7:7)))
    return
20  if(name(6:6) .eq. '9') go to 30
    name(6:6) = char(1 + ichar(name(6:6)))
    name(7:7) = '0'
    return
30  if(name(5:5) .eq. '9') stop
    name(5:5) = char(1 + ichar(name(5:5)))
    name(6:6) = '0'
    name(7:7) = '0'
    return
  end subroutine incnam

end program mdot_vs_time
