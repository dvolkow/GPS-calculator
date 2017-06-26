!---программа тестирования процедур определения координат спутников
!---и координат места наблюдения---выполнил Волков Д.В.--АО СПбГУ--

program main_gps
	implicit none
	!---бортовые эфемериды:
	real(8) :: Toe, WN, ec, A0_5, Omega_0, i_0, w, M_0, delta_n, Omega_vel, IDOT
	real(8) :: C_uc, C_us, C_rc, C_rs, C_ic, C_is
	!---наблюдательные данные:
	real(8) :: P, t_obs
	!---результирующие переменные:
	real(8) :: X, Y, Z  !---положение в земной системе координат
	real(8), dimension(4) :: mX, mY, mZ, mP_0, mt 
	real(8) :: lon, lat, dt, k, ln_SEC, lt_SEC
	integer ln_DEG, ln_M,  lt_DEG, lt_M 

	!---чтение RINEX-файла (здесь должна быть процедура его чтения,
	!---которая пока что симулируется чтением обычного файла с дан-
	!---ными из пособия)-------------------------------------------
	open(unit=10, file='data.dat')
	read(10,*) P, t_obs
	read(10,*) Toe
	read(10,*) WN
	read(10,*) ec
	read(10,*) A0_5
	read(10,*) Omega_0
	read(10,*) i_0
	read(10,*) w
	read(10,*) M_0
	read(10,*) delta_n
	read(10,*) Omega_vel
	read(10,*) IDOT
	read(10,*) C_uc, C_us
	read(10,*) C_rc, C_rs
	read(10,*) C_ic, C_is
	
	!---вызов процедуры, опр. коорд. спутника в земной СК
	call coord_sp(Toe,WN,ec,A0_5,Omega_0,i_0,w,M_0,delta_n,Omega_vel,IDOT,C_uc,C_us,C_rc,C_rs,C_ic,C_is,P,t_obs,X,Y,Z)

	!---определение координат места и поправки часов
	read(10,*) mX(:)
	read(10,*) mY(:)
	read(10,*) mZ(:)
	read(10,*) mP_0(:)
	read(10,*) mt(:)

	call geo_coord(mX, mY, mZ, mP_0, mt, lon, lat, dt)

	k=180d0/(4d0*atan(1.d0))

!	write(*,*) lon*k, lat*k, dt
	write(*,*) "----------------------------------------"
	write(*,*) "Поправка часов = ", dt
	call rad_to_deg(lon, ln_DEG, ln_M, ln_SEC)
	call rad_to_deg(lat, lt_DEG, lt_M, lt_SEC)
	write(*,*) "Широта ",lt_DEG," градусов ",lt_M," минут ",lt_SEC," секунд"
	write(*,*) "Долгота ",ln_DEG," градусов ",ln_M," минут ",ln_SEC," секунд"

end

subroutine rad_to_deg(rad, deg, mint, sec)
	implicit none
	integer deg, mint
	real(8) rad, sec

	sec = 206264.8d0*rad
	deg = INT(sec / 3600)
	sec = sec - 3600*deg
	mint = INT(sec / 60)
	sec = sec - 60*mint
endsubroutine rad_to_deg	
