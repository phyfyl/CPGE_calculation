  subroutine cpge
     !> Calculate CPGE tensor

     use wmpi
     use para
     implicit none
    
     integer :: iR, ik, ikx, iky, ikz
     integer :: m, n, i, j, ie,jf
     integer :: ierr, knv3

     real(dp) :: kdotr, mu,freq,Beta_fake
     real(dp) :: k(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_bulk(:, :)
     complex(dp), allocatable :: velocity_wann(:, :, :), velocity_Ham(:, :, :)

     complex(dp) :: ratio

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:),freqs(:)
     real(dp), allocatable :: beta_op(:,:, :,:)
     real(dp), allocatable :: beta_op_mpi(:,:, :,:)
     real(dp) :: beta_op(3,3,OmegaNum,freqnum)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi

     real(dp),allocatable :: beta_temp( :, :)

     allocate(beta_temp(3,3))
     allocate( W (Num_wann))
     allocate( energy(OmegaNum))
     allocate( freqs(freqnum))
     allocate( beta_op(3,3,omeganum,freqnum))
     allocate( beta_op_mpi(3,3,omeganum,freqnum))
     allocate(velocity_wann(Num_wann, Num_wann, 3), velocity_Ham(Num_wann,Num_wann, 3))
     allocate(UU(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann),Hamk_bulk(Num_wann, Num_wann))
     allocate(Amat(Num_wann, Num_wann))
     W=0d0; velocity_wann= 0d0; UU= 0d0; UU_dag= 0d0; velocity_Ham= 0d0
     beta_op    = 0d0
     beta_op_mpi= 0d0
     Hamk_bulk=0d0
     
     !> energy
     do ie=1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo ! ie
     
     do ie=1, freqnum
        if (freqnum>1) then
           freqs(ie)= freqmin+ (freqmax-freqmin)* (ie-1d0)/dble(freqnum-1)
        else
           freqs= freqmin
        endif
     enddo ! ie

     knv3= Nk1*Nk2*Nk3
     beta_k_mpi=0d0
     beta_k=0d0
     call now(time_start) 
     do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
           call now(time_end) 
           write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
           ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
           time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, Hamk_bulk)
   
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
  
        !> get velocity operator in Hamiltonian basis
        UU_dag= conjg(transpose(UU))

        call dHdk_atomicgauge_new(k, velocity_wann)

        !> unitility rotate velocity
        do i=1, 3
           call mat_mul(Num_wann, velocity_wann(:, :, i), UU, Amat)
           call mat_mul(Num_wann, UU_dag, Amat, velocity_Ham(:, :, i))
        enddo

        do ie=1, OmegaNum
           mu = energy(ie)
           do jf=1, freqnum
           beta_temp=0d0
           call beta_k_cal(mu,freqs(jf),W, velocity_Ham,beta_temp)
           beta_op_mpi(:,:, ie, jf)=beta_op_mpi(:,:, ie,jf)+beta_temp(:,:)
           enddo ! jf
        enddo ! ie
     enddo ! ik
#if defined (MPI)
     call mpi_allreduce(beta_op_mpi,beta_op,size(beta_op),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     beta_op= beta_op_mpi
#endif

     beta_op=beta_op/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*(pi*Echarge**3)/(4*hbar**2)
 
    outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        open(unit=outfileindex, file='beta.txt')
        write(outfileindex, '("#",a)')' Circular Photogalvanic Effect Tensor'
        write(outfileindex, '("#",a13, 20a16)')'freqs(eV)', '\beta'
        do ie=1, omeganum
        do jf=1,freqnum
           write(outfileindex,'(200E30.8)')energy(ie)/eV2Hartree,freqs(jf)/eV2Hartree,beta_op(1,1,ie,jf),&
            beta_op(1,2,ie,jf), beta_op(1,3,ie,jf),&
             beta_op(2,1,ie,jf),beta_op(2,2,ie,jf),beta_op(2,3,ie,jf),&
             beta_op(3,1,ie,jf),beta_op(3,2,ie,jf),beta_op(3,3,ie,jf)
        enddo
        enddo
        close(outfileindex)
     endif

     deallocate( W, Hamk_bulk, UU, energy)
     deallocate( beta_temp, beta_op, beta_op_mpi)
 
     return
  end subroutine cpge


  subroutine beta_k_cal(mu,omega1,W,Vmn_Ham,beta_op)
     use wmpi
     use para
     implicit none

     !> D_mn^H=V_mn/(En-Em) for m!=n
     !> D_nn^H=0 
     complex(dp), intent(in) :: Vmn_Ham(Num_wann, Num_wann,3)
     real(dp), intent(in) :: W(Num_wann)
     real(dp),intent(in) :: mu,omega1

     !> Berry curvature vectors for all bands
     real(dp), intent(out) :: beta_op(3,3)

     real(dp) :: Beta_fake,delta_nm,ff(Num_wann),fe,fnm,lor_f,Vnn(3)
     real(dp), external :: fermi,fermi_E
     integer :: m, n, l,i,j,k
     complex(dp) :: Vnm(3)

     beta_op=0d0
     Beta_fake=1d0/Eta_Arc

     do n=1,Num_wann
         ff(n)=fermi(W(n)-mu, Beta_fake)
     enddo

       do m=1,Num_wann-1
          if (ff(m)<eps6) exit
          do n=m+1,Num_wann
            delta_nm=W(n)-W(m)
            if (delta_nm<eps6 .or. abs(ff(n)-1d0)<eps6) cycle
            fnm=(ff(n)-ff(m))/delta_nm**2 
            lor_f=delta_op/pi/((delta_nm-omega1)**2+delta_op**2)
            Vnn(:)=Vmn_Ham(m,m,:)-Vmn_Ham(n,n,:)       

            Vnm(1)=Vmn_Ham(n,m,3)*Vmn_Ham(m,n,2)
            Vnm(2)=Vmn_Ham(n,m,1)*Vmn_Ham(m,n,3)
            Vnm(3)=Vmn_Ham(n,m,2)*Vmn_Ham(m,n,1)

            do l=1,3 
               beta_op(l,:)=beta_op(l,:)+aimag(fnm*Vnn(l)*Vnm*lor_f)
            enddo
          enddo
        enddo

     return
  end subroutine beta_k_cal

  subroutine dHdk_atomicgauge_new(k, velocity_Wannier)
   !> Velocity operator in Wannier basis using atomic gauge
   !> First derivate of H(k); dH(k)/dk
   use para, only : Nrpts, irvec, crvec,Origin_cell, HmnR, ndegen, &
       Num_wann, dp, Rcut, pi,  &
      zi, soc, zzero
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> velocity operator in Wannier basis using atomic gauge 
   complex(dp), intent(out) :: velocity_Wannier(Num_wann, Num_wann, 3)
   complex(dp) :: vnm_R(Num_wann, Num_wann, 3),vnm_o(Num_wann,Num_wann),vnm_L(Num_wann,Num_wann,3)

   integer :: iR, N0,i1, i2, i

   real(dp) :: pos1_direct(3), pos2_direct(3)
   real(dp) :: pos1_cart(3),pos2_cart(3),kdotrm(Num_wann),ratiom(Num_wann)
   real(dp) :: kdotr, dis_direct(3),dis_cart(3),kpi(3)
   complex(dp) :: ratio,factor(Num_wann,Num_wann)

   vnm_R=zzero
   vnm_o=zzero
   vnm_L=zzero
   kpi=2*pi*k

   do iR=1,(Nrpts-1)/2

      kdotr = kpi(1)*irvec(1,iR)+kpi(2)*irvec(2,iR)+kpi(3)*irvec(3,iR)
      ratio = (-sin(kdotr)+zi*cos(kdotr))/ndegen(iR)
      factor=HmnR(:, :, iR)*ratio

      do i=1,3
      vnm_R(:, :, i)=vnm_R(:, :, i)+crvec(i,iR)*factor
      enddo
      vnm_o=vnm_o+factor

    enddo ! iR

   do i=1,3
   vnm_R(:,:,i)=vnm_R(:,:,i)+conjg(transpose(vnm_R(:,:,i)))
   enddo
   vnm_o=vnm_o-conjg(transpose(vnm_o))

   N0=(Nrpts+1)/2
   vnm_o=vnm_o+zi*HmnR(:,:,N0)/ndegen(N0)

   vnm_L=0d0
   do i1=1, Num_wann-1
      pos1_direct= Origin_cell%wannier_centers_direct(:, i1)
      pos1_cart= Origin_cell%wannier_centers_cart(:, i1)
      do i2=i1+1,Num_wann
         pos2_direct= Origin_cell%wannier_centers_direct(:, i2)
         pos2_cart= Origin_cell%wannier_centers_cart(:, i2)
         dis_direct=pos2_direct-pos1_direct
         dis_cart = pos2_cart-pos1_cart
         kdotr = kpi(1)*dis_direct(1) + kpi(2)*dis_direct(2) + kpi(3)*dis_direct(3)
         ratio = cos(kdotr)+zi*sin(kdotr)
         vnm_L(i1, i2, :)=(vnm_R(i1,i2,:)+dis_cart*vnm_o(i1,i2))*ratio
      enddo ! i2
   enddo ! i1

  do i=1,3
     velocity_Wannier(:,:,i) = vnm_L(:,:,i)+conjg(transpose(vnm_L(:,:,i)))
  enddo
  do i1=1,Num_wann
     velocity_Wannier(i1,i1,:)=vnm_R(i1,i1,:)
  enddo

end subroutine dHdk_atomicgauge_new
