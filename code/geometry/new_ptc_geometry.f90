program new_ptc_geometry
use madx_ptc_module
use pointer_lattice
implicit none

logical(lp) :: doit
integer :: i, j, mf, pos, example,poss
real(dp) :: b0
real(dp), dimension(3) :: a, d
real(dp), dimension(6) :: fix1, fix2, mis, x
type(real_8), dimension(6) :: y1, y2
type(layout), pointer :: L1, L2, L3, L4, L5, L6
type(layout), pointer :: PSR1, PSR2, Fig8, Col1, Col2
type(fibre), pointer :: p1, p2, b, f
type(internal_state) :: state
type(fibre_appearance), pointer :: next_doko
! some simple fitting objects for the fake collider complex
type(pol_block) :: qf(2), qd(2)
type(normalform) :: n1, n2
type(damap) :: id
type(taylor) :: eq(4)
type(gmap) :: g
!-----------------------------------
!-----------------------------------
    interface
       subroutine build_PSR(PSR)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: PSR
       end subroutine build_PSR
       subroutine build_Quad_for_Bend(PSR)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: PSR
       end subroutine build_Quad_for_Bend
       subroutine build_PSR_minus(PSR)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: PSR
       end subroutine build_PSR_minus
    end interface

call ptc_ini_no_append
write(6,*) "Do you want to read an accelerator complex ? yes=1 "
read(5,*) poss
if(poss==1) then
 call read_universe_database(m_u,'m_u.txt',arpent=my_false)
 call read_universe_pointed(M_u,M_t,'m_t.txt')
 call create_dna(M_u,m_t); col2=>m_t%end; ! col2 pointed needed
 goto 999
end if


!==========================!
!== set up DNA sequences ==!
!==========================!
call append_empty_layout(m_u) ! DNA sequence 1
call set_up(m_u%end)
L1 => m_u%end
call build_PSR(L1); L1%name="L1";

call append_empty_layout(m_u) ! DNA sequence 2
call set_up(m_u%end)
L2 => m_u%end
call  build_Quad_for_Bend(L2); L2%name="L2";

call append_empty_layout(m_u) ! DNA sequence 3
call set_up(m_u%end)
L3 => m_u%end
call build_PSR_minus(L3); L3%name="L3";

call append_empty_layout(m_u) ! DNA sequence 4
call set_up(m_u%end)
L4 => m_u%end
call build_PSR(L4); L4%name="L4";

call append_empty_layout(m_u) ! DNA sequence 5
call set_up(m_u%end)
L5 => m_u%end
call build_PSR_minus(L5); L5%name="L5";

call append_empty_layout(m_u) ! DNA sequence 6
call set_up(m_u%end)
L6 => m_u%end
call build_PSR(L6); L6%name="L6";

!================================!
!== create "trackable" layouts ==!
!================================!

!== PSR1 : forward ring  (layout 1 of m_t)
call append_empty_layout(m_t)
PSR1 => m_t%end

p1 => L1%start
p2 => L2%start
do i = 1, L1%n
  if(p1%mag%name == "B") then
    ! read bends from L2
    call append_point(PSR1, p2)
    f => PSR1%end
    d = p1%chart%f%o - f%chart%f%o
    call translate(f, d)
    call compute_entrance_angle(f%chart%f%mid, p1%chart%f%mid, a)
    call rotate(f, f%chart%f%o, a, basis = f%chart%f%mid)
    p2 => p2%next
  else
    call append_point(PSR1, p1)
  end if
  p1 => p1%next
end do ! elements in PSR1 now in correct locations

f => PSR1%start
do i = 1, PSR1%n
  if(f%mag%name == "B_QUAD") then
    call find_patch(f%previous, f, next = .true.)
    call find_patch(f, f%next, next = .false.)
  end if
  f => f%next
end do ! PSR1 now patched

PSR1%name = "PSR 1"
PSR1%closed = .true.
call ring_L(PSR1, .true.) ! make it a ring topologically

!== PSR2 : backward ring (layout 2 of m_t)
call append_empty_layout(m_t)
PSR2 => m_t%end

p1 => L1%end
p2 => L2%end
do i = 1, L1%n
  if(p1%mag%name == "B") then
    call append_point(PSR2, p2)
    p2 => p2%previous
  else
    call append_point(PSR2, p1)
  end if
  f => PSR2%end
  f%dir = -1
  f%charge = -1
  p1 => p1%previous
end do

f => PSR2%start
do i = 1, PSR2%n
  if(f%mag%name == "B_QUAD") then
    call find_patch(f%previous, f, next = .true.)
    call find_patch(f, f%next, next = .false.)
  end if
  f => f%next
end do

PSR2%name = "PSR 2"
PSR2%closed = .true.
call ring_l(PSR2, .true.) ! make it a ring topologically


!== Fig8 : figure-eight lattice (layout 3 of m_t)
d = zero
d(3) = -40.d0
call translate(L4, d)
a = zero
a(2) = pi
call rotate(L3, L3%start%chart%f%a, a)
call move_to(L4, p1, "B", pos)
d = p1%chart%f%a - L3%end%chart%f%b
call translate(L3, d)

call append_empty_layout(m_t)
Fig8 => m_t%end
p1 => L4%start
do i = 1, L4%n
  call append_point(Fig8, p1)
  p1 => p1%next
end do

call append_point(Fig8, p1)
p1 => p1%next
call append_point(Fig8, p1)
p1 => p1%next
call append_point(Fig8, p1)

p1 => L3%end
do i = 1, L3%n
  call append_point(Fig8, p1)
  Fig8%end%dir = -1
  if(p1%mag%name == "B") call add(p1,1,0, -p1%mag%bn(1))
  p1 => p1%previous
end do

p1 => L4%end%previous%previous
call append_point(Fig8, p1)
p1 => p1%next
call append_point(Fig8, p1)
p1 => p1%next
call append_point(Fig8, p1)

Fig8%name = "Figure-Eight"
Fig8%closed = .true.
call ring_l(Fig8, .true.) ! make it topologically closed

p1 => Fig8%start
do i = 1, Fig8%n
  call check_need_patch(p1, p1%next, 1.d-10, pos)
  if(pos /= 0) call find_patch(p1, p1%next, next = .false.)
  p1 => p1%next
end do

!== Col1 : lower collider ring (m_t layout 4)
!== Col2 : upper collider ring (m_t layout 5)
d = zero
d(3) = 40.d0
call translate(L6, d)
a = zero
a(2) = pi
call rotate(L5, L5%start%chart%f%a, a)
call move_to(L6, p1, "B", pos)
d = p1%chart%f%a - L5%end%chart%f%b
call translate(L5, d)

call append_empty_layout(m_t)
Col1 => m_t%end
p1 => L6%start
do i = 1, L6%n
  call append_point(Col1, p1)
  p1 => p1%next
end do

Col1%name = "Collider 1"
Col1%closed = .true.
call ring_l(Col1, .true.) ! make it a ring topologically

call append_empty_layout(m_t)
Col2 => m_t%end
p1 => L6%start%next%next
do i = 1, 6
  call append_point(Col2, p1)
  Col2%end%dir = -1
  p1 => p1%previous
end do
p1 => L5%start
do i = 1, L3%n
  call append_point(Col2, p1)
  p1 => p1%next
end do

Col2%name = "Collider 2"
Col2%closed = .true.
call ring_l(Col2, .true.) ! make it a ring topologically

p1 => Col2%start
do i = 1, Col2%n
  call check_need_patch(p1, p1%next, 1.d-10, pos)
  if(pos /= 0) call find_patch(p1, p1%next, next = .false.)
  p1 => p1%next
end do


!=======================!
!== set up DNA arrays ==!
!=======================!

allocate(PSR1%DNA(2))
PSR1%DNA(1)%L => L1
do i = 2, 2
  PSR1%DNA(i)%L => PSR1%DNA(i-1)%L%next ! L2
end do

allocate(PSR2%DNA(2))
PSR2%DNA(1)%L => L1
do i = 2, 2
  PSR2%DNA(i)%L => PSR2%DNA(i-1)%L%next ! L2
end do

allocate(Fig8%DNA(2))
Fig8%DNA(1)%L => L3
do i = 2, 2
  Fig8%DNA(i)%L => Fig8%DNA(i-1)%L%next ! L4
end do

allocate(Col1%DNA(2))
Col1%DNA(1)%L => L5
do i = 2, 2
  Col1%DNA(i)%L => Col1%DNA(i-1)%L%next ! L6
end do

allocate(Col2%DNA(2))
Col2%DNA(1)%L => L5
do i = 2, 2
  Col2%DNA(i)%L => Col2%DNA(i-1)%L%next ! L6
end do


!== maps with polymorphs
!== simple fit

qf(1) = 0
qf(1)%name = "qf"
qf(1)%ibn(2) = 1
qf(2) = 0
qf(2)%name = "qf"
qf(2)%ibn(2) = 3
qd(1) = 0
qd(1)%name = "qd"
qd(1)%ibn(2) = 2
qd(2) = 0
qd(2)%name = "qd"
qd(2)%ibn(2) = 4
Col1%dna(1)%L = qf(1)
Col1%dna(1)%L = qd(1)
CoL1%dna(2)%L = qf(2)
CoL1%dna(2)%L = qd(2)

101 continue
state = default0 + only_4d0

fix1 = 0.d0
fix2 = 0.d0;
call init(state, 2, c_%np_pol) ! c_%np_pol is automatically computed
call find_orbit(CoL1, fix1, 1, state, 1.d-6)
call find_orbit(Col2, fix2, 1, state, 1.d-6)
call alloc(y1)
call alloc(y2)
call alloc(id)
call alloc(n1)
call alloc(n2)
call alloc(eq);
id=1 ! identity damap
y1 = id + fix1 ! this is permitted in ptc only (not fpp)
y2 = id + fix2 ! closed orbit added to map
call track(Col1, y1, 1, +state) ! unary + activates knobs
call track(Col2, y2, 1, +state)
n1 = y1 ! normal forms: abused of language permitted by ptc
n2 = y2 ! normally one should do => damap=y; normalform=damap
write(6,*) " tunes 1 "
write(6,*) n1%tune(1:2)
write(6,*) " tunes 2 "
write(6,*) n2%tune(1:2)
eq(1) = n1%dhdj%v(1) - 0.254d0
eq(2) = n1%dhdj%v(2) - 0.255d0
eq(3) = n2%dhdj%v(1) - 0.130d0
eq(4) = n2%dhdj%v(2) - 0.360d0
do i = 1, 4
  eq(i) = eq(i) <= c_%npara
end do

call kanalnummer(mf,"eq.txt")
do i=1,4
  call daprint(eq(i), mf)
end do
close(mf)

call kill(y1)
call kill(y2)
call kill(id)
call kill(n1)
call kill(n2)
call kill(eq)
call init(1,4)
call alloc(g,4)
call kanalnummer(mf,"eq.txt")
do i = 1, 4
  call read(g%v(i), mf)
end do
close(mf)

g = g.oo.(-1)
tpsafit(1:4) = g
set_tpsafit = .true.
set_element = .true.
Col1%dna(1)%L = qf(1)
Col1%dna(1)%L = qd(1)
Col1%dna(2)%L = qf(2)
Col1%dna(2)%L = qd(2)
set_tpsafit = .false.
set_element = .false.
call kill(g)
write(6,*) " more "
read(5,*) i
if(i == 1) goto 101
call kill_para(Col1%dna(1)%l)
call kill_para(Col1%dna(2)%l)

!===============================!
!== set up Siamese and Girder ==!
!===============================!

call move_to(Col1, p1, 67)
call move_to(Col2, p2, 7)
p1%mag%siamese => p2%mag
p2%mag%siamese => p1%mag
call move_to(Col1, p1, 4)
call move_to(Col2, p2, 70)
p1%mag%siamese => p2%mag
p2%mag%siamese => p1%mag

call move_to(Col1, p1, 68)
f => p1  ! remember start of girder linked-list
do i = 2, 7
  p2 => p1%next
  p1%mag%girders => p2%mag
  p1 => p1%next
end do
call move_to(Col2, p2,  7)
p1%mag%girders => p1%mag%siamese
p1%mag%siamese%girders => p2%mag
p2%mag%girders => p2%mag%siamese
call move_to(Col1, p1, 67)
call move_to(Col2, p2, 14)
p1%mag%girders => p2%mag
p2%mag%girders => f%mag

call move_to(Col1, p1, 1)
call alloc_af(p1%mag%girder_frame, girder = .true.)
p1%mag%girder_frame%ent = p1%mag%parent_fibre%chart%f%ent
p1%mag%girder_frame%a   = p1%mag%parent_fibre%chart%f%a
p1%mag%girder_frame%exi = p1%mag%parent_fibre%chart%f%ent
p1%mag%girder_frame%b   = p1%mag%parent_fibre%chart%f%a

a = p1%mag%girder_frame%a
a(3) = a(3) - 5.d0
call move_to(Col1, b, 67)
call alloc_af(b%mag%siamese_frame)
call find_patch(b%mag%p%f%a, b%mag%p%f%ent, &
                a, p1%mag%girder_frame%ent, &
                b%mag%siamese_frame%d, b%mag%siamese_frame%angle)
a = p1%mag%girder_frame%a
a(3) = a(3) + 5.d0
call move_to(Col1, b, 4)
call alloc_af(b%mag%siamese_frame)
call find_patch(b%mag%p%f%a, b%mag%p%f%ent, &
                a, p1%mag%girder_frame%ent, &
                b%mag%siamese_frame%d, b%mag%siamese_frame%angle)


!===========================!
!== example misalignments ==!
!===========================!

999 continue
write(6,*) "Example # (from the manual) 1--11 ?"
write(6,*) "Input 0 and then nothing is done!"
read(5,*) example

call move_to(Col2, p2, 7)
if(example == 1) then
  mis = 0.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis)
elseif(example == 2) then
  mis = 0.d0
  mis(1) = 2.0d0
  call misalign_girder(p2, mis)
elseif(example == 3) then
  mis = 0.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis)
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_girder(p2, mis, add = .false.)
elseif(example == 4) then
  mis = 0.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis)
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_girder(p2, mis, add = .true.)
elseif(example == 5) then
  mis = 0.d0
  mis(1) = 2.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis)
elseif(example == 6) then
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_siamese(p2, mis)
elseif(example == 7) then
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_siamese(p2, mis)
  mis = 0.d0
  mis(1) = 2.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis, add = .false.)
elseif(example == 8) then
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_siamese(p2, mis)
  mis = 0.d0
  mis(1) = 2.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis, add = .true.)
elseif(example == 9) then
  mis = 0.d0
  mis(1) = 2.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis)
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_siamese(p2, mis, add = .false.)
elseif(example == 10) then
  mis = 0.d0
  mis(1) = 2.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis)
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_siamese(p2, mis, add = .true.)
elseif(example == 11) then
  mis = 0.d0
  mis(1) = 2.d0
  mis(5) = pi / 8.d0
  call misalign_girder(p2, mis)
  mis = 0.d0
  mis(1) = 2.d0
  call misalign_siamese(p2, mis, add = .false., &
                        preserve_girder = .true.)
end if

call kanalnummer(mf,"doko.txt")
call TIE_MAD_UNIVERSE(m_u)
p1=>m_u%start%start
do i=1,m_u%nf
 write(mf,'(i4,1x,a16)') i,p1%mag%name
 next_doko=>p1%mag%doko
  do while(associated(next_doko))
   write(mf,'(a120)') next_doko%parent_fibre%parent_layout%name 
   p2=>next_doko%parent_fibre
   next_doko=>next_doko%next
 enddo
 p1=>p1%n
end do
close(mf)
 
l1=>m_t%start
do i=1, m_t%n
 write(6,*) l1%name
 do j=1,size(l1%dna)
  write(6,*) i,l1%dna(j)%l%name
 enddo
 l1=>l1%next
enddo

call ptc_end(graphics_maybe=2,flat_file=.true.)

end program new_ptc_geometry


!=================================================================
subroutine  build_PSR(PSR)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: PSR

real(dp) :: ang, brho, kd, kf, Larc
type(fibre) :: b, d1, d2, qd, qf
type(layout) :: cell
!-----------------------------------

call make_states(.false.)
exact_model = .true.
madlength = .false.

ang = (twopi * 36.d0 / 360.d0)
Larc = 2.54948d0
brho = 1.2d0 * (Larc / ang)
call set_mad(brho = brho, method = 6, step = 10)
madkind2 = drift_kick_drift

kf =  2.72d0 / brho
kd = -1.92d0 / brho

d1 = drift("D1", 2.28646d0)
d2 = drift("D2", 0.45d0)
qf = quadrupole("QF", 0.5d0, kf)
qd = quadrupole("QD", 0.5d0, kd)
b  = rbend("B", Larc, ang)
cell = d1 + qd + d2 + b + d2 + qf + d1

PSR = 10 * cell
PSR = .ring.PSR

call survey(PSR)
end subroutine build_PSR


!=================================================================
subroutine  build_PSR_minus(PSR)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: PSR

real(dp) :: ang, brho, kd, kf, Larc
type(fibre) :: b, d1, d2, qd, qf
type(layout) :: cell
!-----------------------------------

call make_states(.false.)
exact_model = .true.
madlength = .false.

ang = (twopi * 36.d0 / 360.d0)
Larc = 2.54948d0
brho = 1.2d0 * (Larc / ang)
call set_mad(brho = brho, method = 6, step = 10)
madkind2 = drift_kick_drift

kf =  2.72d0 / brho
kd = -1.92d0 / brho

d1 = drift("D1", 2.28646d0)
d2 = drift("D2", 0.45d0)
qf = quadrupole("QF", 0.5d0, kf)
qd = quadrupole("QD", 0.5d0, kd)
b  = rbend("B", Larc, ang)
cell = d1 + qd + d2 + b + d2 + qf + d1

PSR = b + d2 + qf + d1 + 8 * cell + d1 + qd + d2 + b
PSR = .ring.PSR

call survey(PSR)
end subroutine build_PSR_minus

!=================================================================

subroutine  build_Quad_for_Bend(PSR)
use madx_ptc_module
use pointer_lattice
implicit none
    type(keywords) key

type(layout),target :: PSR
integer code,EXCEPTION,ik
logical(lp) doneit
real(dp) :: ang, ang2, brho, b1, Larc, Lq
type(fibre) :: b
!-----------------------------------

call make_states(.false.)
exact_model = .true.
default = default + nocavity  
call update_states
madlength = .false.

ang = (twopi * 36.d0 / 360.d0)
Larc = 2.54948d0
brho = 1.2d0 * (Larc / ang)
madkind2 = drift_kick_drift
ang2 = ang / two
b1 = ang / Larc
Lq = Larc * sin(ang2) / ang2

    CALL SET_MADx(brho=brho,METHOD=6,STEP=10)

code=2
    do ik=1,10
    call zero_key(key)
    call append_empty(psr)
    
    key%list%name="B_QUAD"
    key%list%l=Lq

    select case(code)
    case(2)
       key%magnet="quadrupole"
       key%list%k(2)=0.d0
       key%list%ks(2)=0.d0
       key%list%k(1)=b1
       key%list%ks(1)=0.d0
       key%list%BEND_FRINGE=.true.
       key%list%PERMFRINGE=.true.
       key%model= "DRIFT_KICK       "

    end select
    
        call create_fibre(psr%end,key,EXCEPTION)

    enddo
   
    psr%closed=.true.

    doneit=.true.
    call ring_l(psr,doneit)

call survey(PSR)
end subroutine build_Quad_for_Bend
