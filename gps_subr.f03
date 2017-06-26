!-------Определение координат спутника по бортовым эфемеридам--------
!---------------------исп. переменные:-------------------------------
!-------Бортовые Эфемериды: Toe - эпоха
!-------------------------- WN  - неделя
!-------------------------- ec   - эксцентриситет орбиты
!-------------------------- A0_5 - квадр. корень из большой полуоси 
!-------------------------- Omega_0 - долгота восходящего узла
!-------------------------- i_0 - наклон орбиты
!-------------------------- w   - аргумент перигея
!-------------------------- M_0 - средняя аномалия на начало эпохи
!-------------------------- delta_n - отклонение значения среднего движения
!-------------------------- Omega_vel - скорость изменения долготы в. узла
!-------------------------- IDOT - скорость измен. наклона орбиты
!-------------------------- C_uc - амплитуда квадр. поправки арг. широты
!-------------------------- C_us - амплитуда синфазн. поправки арг. широты
!-------------------------- C_rc - амплитуда квадрю поправки радиуса орбиты
!-------------------------- C_rs - амплитуда синфазн. поправки радиуса орбиты
!-------------------------- C_ic - амплитуда квадрю. поправки наклона орбиты
!-------------------------- C_is - амплитуда синфазн. поправки наклона орбиты
!-------Наблюдательные дан: P    - псеводальность до спутника
!-------------------------- t_obs - момент наблюдения
!----Фундаментальные пост.: mu   - геоцентрическая грав. постоянная
!-------------------------- Omega_vel_e - угловая скорость вр. Земли
!-------------------------- c    - скорость света
!-------------------------- a    - большая полуось эллипсоида WGS-84



subroutine coord_sp(Toe,WN,ec,A0_5,Omega_0,i_0,w,M_0,delta_n,Omega_vel,IDOT,C_uc,C_us,C_rc,C_rs,C_ic,C_is,P,t_obs,X,Y,Z)
	implicit none
	!---бортовые эфемериды
	real(8), intent(in) :: Toe, WN, ec, A0_5, Omega_0, i_0, w, M_0, delta_n, Omega_vel, IDOT
	real(8), intent(in) :: C_uc, C_us, C_rc, C_rs, C_ic, C_is
	!---наблюдательные данные
	real(8), intent(in) :: P, t_obs
	!---фундаментальные постоянные
	real(8), parameter :: mu=3.986004418e14, Omega_vel_e=7.2921151467e-5, c=299792458.0, a=6378137.0 
	!---вспомогательные переменные
	real(8) :: n_0 !-------среднее движение
	real(8) :: t   !-------время от начала эпохи
	real(8) :: n   !-------исправленное среднее движение
	real(8) :: M   !-------средняя аномалия
	real(8) :: E   !-------эксцентрическая аномалия
	real(8) :: nu  !-------истинная аномалия
	real(8) :: fi  !-------аргумент широты
	real(8) :: delta_u, delta_r, delta_i !--поправки широты, радиуса и наклона
	real(8) :: u, r, i !---испр. широта, радиус и наклон
	real(8) :: X_orb, Y_orb  !----координаты в орбитальной плоскости
	real(8) :: Omega !---испр. долгота восх. узла
	real(8), intent(out) :: X, Y, Z  !---положение в земной системе координат


	n_0 = sqrt(mu/(A0_5)**6)
    write(*,*) "1. n_0 = ", n_0

	t = t_obs - P/c - Toe
	if ( t > 302400 ) then 
		t = t - 604800
	else if ( t < -302400 ) then
		t = t + 604800
	endif
    write(*,*) "3. t = ", t

	n = n_0 + delta_n
    write(*,*) "4. n = ", t

	M = M_0 + n*t
    write(*,*) "5. M = ", t

	call E_calc(M, ec, E)
    write(*,*) "6. E = ", E

	call nu_calc(ec, E, nu)
    write(*,*) "7. nu = ", nu

	fi = w + nu
    write(*,*) "8. fi = ", fi


	delta_u = C_us*sin(2*fi) + C_uc*cos(2*fi)
	delta_r = C_rs*sin(2*fi) + C_rc*cos(2*fi)
	delta_i = C_is*sin(2*fi) + C_ic*cos(2*fi)
    write(*,*) "9. delta_u  = ", delta_u
    write(*,*) "9. delta_r  = ", delta_r
    write(*,*) "9. delta_i  = ", delta_i

	u = fi + delta_u
	r = A0_5**2*(1 - ec*cos(E)) + delta_r
	i = i_0 + delta_i + IDOT*t
    write(*,*) "10. u  = ", u
    write(*,*) "10. r  = ", r
    write(*,*) "10. i  = ", i

	X_orb = r*cos(u)
	Y_orb = r*sin(u)
    write(*,*) "11. X_orb  = ", X_orb
    write(*,*) "11. Y_orb  = ", Y_orb

	Omega = Omega_0 + t*(Omega_vel - Omega_vel_e) - Toe*Omega_vel_e
    write(*,*) "12. Omega  = ", Omega

	X = X_orb*cos(Omega) - Y_orb*cos(i)*sin(Omega)
	Y = X_orb*sin(Omega) + Y_orb*cos(i)*cos(Omega)
	Z = Y_orb*sin(i)
    write(*,*) "13. X  = ", X
    write(*,*) "13. Y  = ", Y
    write(*,*) "13. Z  = ", Z

endsubroutine coord_sp


!---процедура решения уравнения Кеплера:

subroutine E_calc(M, ec, E)
	implicit none
	real(8), intent(in) :: M, ec
	real(8), intent(out) :: E
	real(8) :: E_old, E_new, funE

	E_old = 0
	E_new = 1

	do while (abs(funE(M,ec,E_new))>0.0000000000000001)  
		E_old = E_new
		E_new = E_old - funE(M,ec,E_old)/(ec*cos(E_old) - 1)	
	enddo

	E = E_new
endsubroutine E_calc

function funE(M, ec, E) result(f)
	implicit none
	real(8) :: M, ec, E, f
	f = M - E + ec*sin(E)
endfunction funE	 

!---процедура нахождения истинной аномалии

subroutine nu_calc(ec, E, nu)
	implicit none
	real(8), intent(in) :: ec, E	
	real(8), intent(out) :: nu

	nu = 2*atan(tan(E/2)*sqrt((1+ec)/(1-ec)) )

endsubroutine nu_calc
