!---Определениe координат спутников--------------------------------
!---и координат места наблюдения---
!-------------Author: Волков Д.В.-------АО мат-меха СПбГУ----------


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
    integer ln_DEG, ln_M,  lt_DEG, lt_M, i

    !---чтение RINEX-файла (здесь должна быть процедура его чтения,
    !---которая пока что симулируется чтением обычного файла с дан-
    !---ными из пособия)-------------------------------------------

    !--1стр---P, t_obs, delta_t------------------------------------
    !--2стр---***----С_rs----delta_n----M_0------------------------
    !--3cтр---С_uc---ec----C_us-------A0_5-------------------------
    !--4стр---T_oe---C_ic----Omega_0----C_is-----------------------
    !--5стр---i_0----C_rc----w--------Omega_vel--------------------
    !--6стр---IDOT---****----WN---------****-----------------------
    
    open(unit=21, file='6_ef.txt')
    open(unit=22, file='17_ef.txt')
    open(unit=23, file='19_ef.txt')
    open(unit=24, file='24_ef.txt')

    do i = 1, 4
        write(*,*) "Результаты № ", i
        read(i+20,*) P, t_obs, mt(i)
        read(i+20,*) C_rs, delta_n, M_0
        read(i+20,*) C_uc, ec, C_us, A0_5
        read(i+20,*) Toe, C_ic, Omega_0, C_is
        read(i+20,*) i_0, C_rc, w, Omega_vel
        read(i+20,*) IDOT, lt_SEC, WN, lt_SEC
        mP_0(i) = P
        !---вызов процедуры, опр. коорд. спутника в земной СК
        call coord_sp(Toe,WN,ec,A0_5,Omega_0,i_0,w,M_0,delta_n,Omega_vel,IDOT,C_uc,C_us,C_rc,C_rs,C_ic,C_is,P,t_obs,X,Y,Z)
        mX(i) = X
        mY(i) = Y
        mZ(i) = Z
        write(*,*) "------------------------------"
    enddo

    write(*,*) "Полная система спутников даёт решение:"
    call geo_coord(mX, mY, mZ, mP_0, mt, lon, lat, dt)

    k=180d0/(4d0*atan(1.d0))

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
