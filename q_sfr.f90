program q_sfr_cal
implicit none
real,dimension(17)::ind,d,ed,b,eb,a,ea,ni,en,ext,f_o,sfr,beta,e_sfr,e_beta,e_ext,e_fo,m,e_m
real,parameter:: pi = 4.0*atan(1.0)
integer::i,n
n=17
open(unit=21,file="q_sfr_cal.dat",action="write",status="new")
open(unit=11,file="q.dat",action="read",status="old")
do i=1,n
read(11,*)ind(i), d(i), ed(i), b(i), eb(i), a(i), ea(i), ni(i), en(i)
end do
do i=1,n
if (b(i) .ne. 0) then
beta(i)=a(i)/b(i)
else 
beta(i)=0
end if
if (beta(i) .le. 2.86) then
sfr(i) = 299.7702*pi*a(i)*d(i)**2/1E11
ext(i)=0
f_o(i)=0
e_sfr(i) = sqrt((sfr(i)*ea(i)/a(i))**2 + (2*sfr(i)*ed(i)/d(i))**2)
e_beta(i)=sqrt((beta(i)*ea(i)/a(i))**2 + (beta(i)*eb(i)/b(i))**2)
e_ext(i)=0.0
e_fo(i)=0.0
else if (beta(i) > 2.86) then
ext(i)=2.45*1.97*log10(beta(i)/2.86)
f_o(i) = a(i)*10**(ext(i)/2.5)
sfr(i) = 299.7702*pi*f_o(i)*d(i)**2/1E11

e_beta(i) = sqrt((beta(i)*ea(i)/a(i))**2 + (beta(i)*eb(i)/b(i))**2)
e_ext(i) = 2.45*1.97*e_beta(i)/beta(i)
e_fo(i) = sqrt((f_o(i)*ea(i)/a(i))**2 + (f_o(i)*log(10.0)*e_ext(i)/2.5)**2)
e_sfr(i) = sqrt((sfr(i)*e_fo(i)/f_o(i))**2 + (2*sfr(i)*ed(i)/d(i))**2)
end if

if (ni(i) >0) then
m(i) = 8.743 + 0.462*log10(ni(i)/a(i))
e_m(i) = 0.462*sqrt((en(i)/ni(i))**2+(ea(i)/a(i))**2)
else if (ni(i).eq.0) then
m(i)=0
e_m(i)=0
end if 
write(21,*)ind(i), beta(i),e_beta(i),ext(i),e_ext(i),f_o(i),e_fo(i),sfr(i),e_sfr(i),m(i),e_m(i)
end do
 end program
