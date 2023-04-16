
c----------------------------------------
c     author kitaki takaaki
c     2019/05/10:for analysis to nacho.
c----------------------------------------

      program main
      implicit none
c---------2D data-------------------------------  inpraのあたいを入れる
      integer,parameter :: nrmax=200,nthmax=240 !2DRHD-data,計算メッシュ
      integer,parameter :: n_file=1             !何個のインプットファイルを用いるか
      real(8),parameter :: r_in_array(n_file) = 2.d0
      real(8),parameter :: r_out_array(n_file) = 6.d3
      real(8),parameter :: MBH_array(n_file)= 1.d1
      character(256) :: infile_array(n_file)=(/ !読み込むファイル
     &     "./m1_mdot1d3_rk1d3_ex.dat"
     &     /)
      character(256) :: outfile_array(n_file)=(/ !出力ファイル
     &     "./rho.txt"
     &     /)

C---------Data(Read) ----------------------------- !重要 読み込むファイル(.dat)からデータを入れる
      real(8) :: dt,lapse
      real(8) :: gase(nrmax,nthmax),rade(nrmax,nthmax)
      real(8) :: rho(nrmax,nthmax),vr(nrmax,nthmax)
      real(8) :: vth(nrmax,nthmax),vphi(nrmax,nthmax)
      real(8) :: Fr(nrmax+1,nthmax),Fth(nrmax,nthmax+1)
      real(8) :: FLD_lambda_r(nrmax-1,nthmax)
      real(8) :: FLD_lambda_th(nrmax,nthmax+1)
      real(8) :: Radforce_r(nrmax,nthmax),Radforce_th(nrmax,nthmax)
      real(8) :: flux_ang(nrmax,nthmax),tgas(nrmax,nthmax)
C---------Physical Constant ----------------------
      real(8),parameter :: msun= 1.989d33    ! solar mass
      real(8),parameter :: c_vel  = 2.998d10 ! light speed
      real(8),parameter :: gc  = 6.672d-8    ! gravity
C---------Parameter (Write)-----------------------
      real(8),parameter :: pi=dacos(-1.d0)
      real(8) :: ra_s(nrmax+1),the_s(nthmax+1)
      real(8) :: ra(nrmax),the(nthmax)
      real(8) :: x(nrmax,nthmax),z(nrmax,nthmax)
      integer :: i,j,id,kd,k,kk,n_div
      character(256) :: outfile,outfile2,infile
      real(8) :: r_in,r_out,MBH,rs,L_Edd
      real(8) :: ra_ini,the_ini,R_ini,z_ini
c----------------------------------------

      do k=1,n_file
         r_in=r_in_array(k)
         r_out=r_out_array(k)
         MBH=MBH_array(k)
         rs = 2.0d0*gc*msun/c_vel**2*MBH
         L_Edd=1.258028968d38*MBH
         infile=infile_array(k)
         outfile=outfile_array(k)
         call mesh(r_in,r_out,nrmax,nthmax,ra_s,the_s,ra,the) !シミュレーションのメッシュを貼る

         open(unit=10,file=infile,status='unknown' !ファイルを読み込む,formはデータの種類
     &        ,form='unformatted')
         read(10) dt,lapse,gase,rade,rho,vr,vth,vphi,
     &        Fr,Fth,FLD_lambda_r,FLD_lambda_th,
     &        Radforce_r,Radforce_th,flux_ang
         close(10)
         call gase_to_tgas(tgas,nrmax,nthmax,gase,rho)
         write(6,*) "infile = ", infile

         write(6,*) dt,lapse
c-----output-----------------------------------
         open(unit=20,file=outfile,status='unknown'
     &        ,form='formatted')     !今回はテクストデータなので、formatted

         do i=1,nrmax
            do j=1,nthmax
               write(20,300) ra(i)*sin(the(j)),ra(i)*cos(the(j)),
     &         rho(i,j) ![rs][rad][g/cm3]
            enddo
            write(20,*) ""  !改行,gnuplotのためだけ
            !write(20,*) ""  !カラーの時はこの文はいらない
         enddo
 300     format(1p,3e21.8)
         close(20)
      enddo


c----------------------------------------
      contains   !以降サブルーチンを作る

      subroutine mesh(r_in,r_out,nrmax,nthmax,ra_s,the_s,ra,the)
      implicit none
c-----local variables----------------------------------------------------
      integer,intent(in) :: nrmax,nthmax
      real(8),intent(in) :: r_in,r_out
      real(8) :: dlogr, logr, dcos_the, cos_the
      real(8) :: mu_th,dmu_th
      integer :: i,j
      real(8),intent(out) :: ra_s(nrmax+1),the_s(nthmax+1)
      real(8),intent(out) :: ra(nrmax),the(nthmax)
C=== Mesh-Construction(Spherical coordinate) ============================
C--- r-direction -------------------
      dlogr = dlog10(r_out) - dlog10(r_in)
      dlogr = dlogr/dble(nrmax) ! grid-spacing for r-direction

      do i=1,nrmax+1
         logr = dlog10(r_in) + dlogr*dble(i-1)
         ra_s(i) = 10.d0**(logr) ! cell surface[rs]
      enddo

      do i=1,nrmax
         ra(i) = 0.5*(ra_s(i)+ra_s(i+1))
      enddo

C---  theta-direction ---------------
C--   theta= 0 - pi/2 -------
      dmu_th = 1.d0/(dble(nthmax)-0.5d0)
      do j = 1,nthmax
         mu_th = 1.d0 - dmu_th*dble(j-1)
         the_s(j) = dacos(mu_th) ! cell surface
      enddo
c      the_s(nthmax+1) = pi/2.0
      the_s(nthmax+1) = pi - the_s(nthmax)
      do j = 1, nthmax
         the(j) = 0.5d0*(the_s(j) + the_s(j+1))! cell center
      enddo

      end


      subroutine gase_to_tgas(tgas,nrmax,nthmax,gase,rho)
      implicit none                                             !型宣言を忘れていた時に教えてくる
      integer,intent(in) :: nrmax,nthmax
      real(8),intent(in) :: gase(nrmax,nthmax),rho(nrmax,nthmax)
      real(8),intent(out) :: tgas(nrmax,nthmax)
      real(8),parameter :: gamma = 5.d0/3.d0 ! spec. heting ratio
      real(8),parameter :: mu = 0.5d0 ! mean moleculer weight
      real(8),parameter :: mp = 1.672649d-24 ! proton mass
      real(8),parameter :: kb  = 1.38066d-16 ! boltzman const.
      integer :: i,j
      do i = 1, nrmax
         do j = 1, nthmax
            tgas(i,j) = gase(i,j)*(gamma-1.d0)*mu*mp/kb/rho(i,j)
         enddo
      enddo
      end                       !calc_tgas

      end                       !main
