!-------Определение координат места и поправки часов приемника-------
!-------------Author: Волков Д.В.-------АО мат-меха СПбГУ----------

subroutine geo_coord(X, Y, Z, P_0, delta_ts, lon, lat, dt)
    implicit none
	real(8), dimension(4), intent(in) :: X, Y, Z, P_0, delta_ts
	real(8), intent(out) :: lon, lat, dt
    integer(2) i
	real(8), dimension(4) :: P
	real(8), parameter :: c = 299792458.d0, a_WGS = 6378137.d0, b_WGS = 6356752.d0, level_spbu=0.0005d0
    !------вспомогательные суммы, разности и пр.:
	real(8), dimension(3) :: delta_P, delta_X, delta_Y, delta_Z, A
	real(8), dimension(3) :: sigma_P, sigma_X, sigma_Y, sigma_Z
    !------рабочая матрица
	real(8), dimension(3,3) :: matrix
	real(8) :: Dbig, Db_1, Db_2, Db_3, d1, d2, d3
	real(8) :: Xa, Ya, Za, ax, ay, az, a4, b4, c4
	real(8) :: eps_1, eps_2, e_t, e_st, X_fin, Y_fin, Z_fin, psi

    forall (i=1:4)
        P(i) = P_0(i) + c*delta_ts(i)
    end forall
    
    do i=1,3
        delta_P(i) = P(i+1) - P(i)
        delta_Z(i) = Z(i+1) - Z(i)
        delta_X(i) = X(i+1) - X(i)
        delta_Y(i) = Y(i+1) - Y(i)

        sigma_P(i) = P(i) + P(i+1)
        sigma_Z(i) = Z(i) + Z(i+1)
        sigma_Y(i) = Y(i) + Y(i+1)
        sigma_X(i) = X(i) + X(i+1)
    
        A(i) = (delta_X(i)*sigma_X(i)+delta_Y(i)*sigma_Y(i)+delta_Z(i)*sigma_Z(i)-delta_P(i)*sigma_P(i))/2
    enddo

    call create_mx(delta_X(1),delta_Y(1),delta_Z(1),delta_X(2),delta_Y(2),delta_Z(2),delta_X(3),delta_Y(3),delta_Z(3),matrix)
    call det3x3(matrix, Dbig)

    call create_mx(A(1),delta_Y(1),delta_Z(1),A(2),delta_Y(2),delta_Z(2),A(3),delta_Y(3),delta_Z(3),matrix)
    call det3x3(matrix, Db_1)

    call create_mx(delta_X(1),A(1),delta_Z(1),delta_X(2),A(2),delta_Z(2),delta_X(3),A(3),delta_Z(3),matrix)
    call det3x3(matrix, Db_2)

    call create_mx(delta_X(1),delta_Y(1),A(1),delta_X(2),delta_Y(2),A(2),delta_X(3),delta_Y(3),A(3),matrix)
    call det3x3(matrix, Db_3)

    call create_mx(-1*delta_P(1),delta_Y(1),delta_Z(1),-1*delta_P(2),delta_Y(2),delta_Z(2),-1*delta_P(3),delta_Y(3),delta_Z(3),matrix)
    call det3x3(matrix, d1)

    call create_mx(delta_X(1),-1*delta_P(1),delta_Z(1),delta_X(2),-1*delta_P(2),delta_Z(2),delta_X(3),-1*delta_P(3),delta_Z(3),matrix)
    call det3x3(matrix, d2)

    call create_mx(delta_X(1),delta_Y(1),-1*delta_P(1),delta_X(2),delta_Y(2),-1*delta_P(2),delta_X(3),delta_Y(3),-1*delta_P(3),matrix)
    call det3x3(matrix, d3)

    Xa = Db_1/Dbig
    Ya = Db_2/Dbig
    Za = Db_3/Dbig

    ax = d1/Dbig
    ay = d2/Dbig
    az = d3/Dbig
    !----коэфф. кв. уравнения:
    a4 = 1 - ax*ax - ay*ay - az*az
    b4 = P(4) + (X(4) - Xa)*ax + (Y(4) - Ya)*ay + (Z(4) - Za)*az
    c4 = P(4)*P(4) - (X(4) - Xa)*(X(4) - Xa) - (Y(4) - Ya)*(Y(4) - Ya) - (Z(4) - Za)*(Z(4) - Za)

    eps_1 = (-1*b4 + sqrt(b4*b4 - a4*c4))/a4
    eps_2 = (-1*b4 - sqrt(b4*b4 - a4*c4))/a4
    if(abs(eps_1/c)>level_spbu) then
        dt = eps_2/c
    else if(abs(eps_2/c)>level_spbu) then
        dt = eps_1/c
    endif

    X_fin = Xa + dt*c*ax
    Y_fin = Ya + dt*c*ay
    Z_fin = Za + dt*c*az

    !----параметры трансформации
    e_t = (a_WGS*a_WGS - b_WGS*b_WGS)/(a_WGS*a_WGS)
    e_st = (a_WGS*a_WGS - b_WGS*b_WGS)/(b_WGS*b_WGS)

    !----долгота
    lon = atan(Y_fin/X_fin)

    !----широта
    psi = atan((a_WGS*Z_fin)/(b_WGS*sqrt(X_fin*X_fin + Y_fin*Y_fin)))
    lat = atan((Z_fin + e_st*b_WGS*(sin(psi)**3))/(sqrt(X_fin*X_fin + Y_fin*Y_fin) - e_t*a_WGS*(cos(psi)**3)))
    
endsubroutine geo_coord

!----процедурка для подсчета определителя 3х3

subroutine det3x3(matr, det_out)
    implicit none
	real(8), dimension(3,3), intent(in) :: matr
	real(8), intent(out) :: det_out

    det_out = matr(1,1)*(matr(2,2)*matr(3,3) - matr(2,3)*matr(3,2)) - matr(1,2)*(matr(2,1)*matr(3,3)-matr(2,3)*matr(3,1))
    det_out = det_out + matr(1,3)*(matr(2,1)*matr(3,2)-matr(2,2)*matr(3,1))

endsubroutine det3x3

!---процедурка заполнения рабочей матрицы

subroutine create_mx(m11,m12,m13,m21,m22,m23,m31,m32,m33,matrix)
    implicit none
	real(8), intent(in) :: m11,m12,m13,m21,m22,m23,m31,m32,m33
	real(8), dimension(3,3), intent(out) :: matrix

    matrix(1,1) = m11
    matrix(1,2) = m12
    matrix(1,3) = m13

    matrix(2,1) = m21
    matrix(2,2) = m22
    matrix(2,3) = m23

    matrix(3,1) = m31
    matrix(3,2) = m32
    matrix(3,3) = m33

endsubroutine create_mx
